version 1.0

workflow OutlierExclusion {
  input {
    # File with list of contigs, one per line, in the JoinRawCalls VCF.
    File contigs_file

    # Output VCF from JoinRawCalls in GATK-SV.
    File join_raw_calls_vcf
    File join_raw_calls_vcf_index

    # See MakeSvCountsDb
    Array[Array[String]] svtypes_to_filter = [["DEL", 5000, 25000], ["DUP", 5000, 25000]]

    # See DetermineOutlierSamples
    File? wgd_scores
    Float? min_wgd_score
    Float? max_wgd_score
    Float iqr_multiplier = 8.0

    # See FormatOutlierSamples and DetermineOutlierVariants
    File? outlier_samples

    # Output VCFs from ClusterBatch in GATK-SV
    Array[File]? clustered_depth_vcfs
    Array[File]? clustered_manta_vcfs
    Array[File]? clustered_wham_vcfs
    Array[File]? clustered_melt_vcfs
    Array[File]? clustered_scramble_vcfs
    # See DetermineOutlierVariants
    Float min_outlier_sample_prop = 1.0

    # Output VCF from FilterGenotypes in GATK-SV
    File filter_genotypes_vcf
    String output_prefix

    # Has all the external dependencies
    String base_docker
    # Adds the OutlierExclusion related scripts
    String pipeline_docker
  }

  Array[String] contigs = read_lines(contigs_file)

  scatter (contig in contigs) {
    call GetContigFromVcf {
      input:
        vcf = join_raw_calls_vcf,
        vcf_index = join_raw_calls_vcf_index,
        contig = contig,
        base_docker = base_docker
    }

    call GetJoinRawCallsClusters {
      input:
        vcf_or_bcf = GetContigFromVcf.contig_bcf,
        base_docker = base_docker
    }

    if (!defined(outlier_samples)) {
      call ConvertVcfOrBcfToTsv {
        input:
          vcf_or_bcf = GetContigFromVcf.contig_bcf,
          base_docker = base_docker
      }
    }
  }

  call GatherSVs {
      input:
          sv_tsvs=select_all(ConvertVcfOrBcfToTsv.tsv)
  }
  
  call GatherClusters {
      input:
          cluster_tsvs=GetJoinRawCallsClusters.clusters
  }

  if (!defined(outlier_samples)) {
    call CountSVs {
      input:
        svs_tsv = GatherSVs.out,
        filters = svtypes_to_filter,
        pipeline_docker = pipeline_docker
    }

    call DetermineOutlierSamples {
      input:
        sv_counts_tsv = CountSVs.sv_counts_tsv,
        sv_filters_tsv = CountSVs.sv_filters_tsv,
        wgd_scores = wgd_scores,
        min_wgd_score = min_wgd_score,
        max_wgd_score = max_wgd_score,
        iqr_multiplier = iqr_multiplier,
        pipeline_docker = pipeline_docker
    }
  }

  File outlier_samples_list = select_first([outlier_samples, DetermineOutlierSamples.sv_count_outlier_samples])
  
  Array[Array[File]] clustered_vcfs = select_all([clustered_depth_vcfs,
    clustered_manta_vcfs,
    clustered_wham_vcfs,
    clustered_melt_vcfs,
    clustered_scramble_vcfs])

  scatter (batch_vcfs in transpose(clustered_vcfs)) {
      call ConvertBatchToTsv {
          input:
              vcfs=batch_vcfs
      }
      
      call DetermineOutlierVariants {
          input:
              variants_tsv = ConvertBatchToTsv.tsv,
              outlier_samples_tsv = outlier_samples_list,
              jrc_clusters_tsv = GatherClusters.out,
              min_outlier_sample_prop = min_outlier_sample_prop,
              pipeline_docker = pipeline_docker
      }
  }
  

  call FlagOutlierVariants {
    input:
      output_prefix = output_prefix,
      filter_genotypes_vcf = filter_genotypes_vcf,
      outlier_variants = DetermineOutlierVariants.outlier_variants,
      pipeline_docker = pipeline_docker
  }

  output {
    File outlier_annotated_vcf = FlagOutlierVariants.outlier_annotated_vcf
    File outlier_annotated_vcf_index = FlagOutlierVariants.outlier_annotated_vcf_index
    File? sv_count_outlier_samples = DetermineOutlierSamples.sv_count_outlier_samples
    File? wgd_outlier_samples = DetermineOutlierSamples.wgd_outlier_samples
  }
}

#=======================================================================
# Extract a single contig from a VCF.
#
# Inputs
# ------
# vcf: VCF file.
# vcf_index: The index file for `vcf`.
# contig: Contig to extract.
# base_docker: Path to Docker image.
#
# Outputs
# -------
# contig_bcf: BCF file only containing records on `contig`.
# contig_bcf_index: Index file for `contig_bcf`.
#=======================================================================
task GetContigFromVcf {
  input {
    File vcf
    File vcf_index
    String contig
    String base_docker
  }

  Int disk_size_gb = ceil(size(vcf, "GB") * 1.2) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  String output_bcf = "${contig}.bcf"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools view --output-type b --regions '~{contig}' --write-index=csi \
      --output '~{output_bcf}' '~{vcf}'
  >>>

  output {
    File contig_bcf = output_bcf
    File contig_bcf_index = "${output_bcf}.csi"
  }
}

#=======================================================================
# Reformat an SV VCF or BCF into a TSV.
#
# Each record is expanded into multiple rows, one for each carrier.
# This is done to facilitate counting SVs per sample, which is more
# complicated when there are multiple carriers per row.
#
# The reformatted TSV will have four columns:
# 1. Variant ID
# 2. SV type
# 3. SV length
# 4. Sample ID
#
# Inputs
# ------
# vcf_or_bcf: File to convert, in VCF or BCF format.
# base_docker: Path to Docker image.
#
# Outputs
# -------
# tsv: Reformatted file compressed with Zstandard.
#=======================================================================
task ConvertVcfOrBcfToTsv {
  input {
    File vcf_or_bcf
    String base_docker
  }

  Int disk_size_gb = ceil(size(vcf_or_bcf, "GB") * 2.0) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  String output_prefix = sub(basename(vcf_or_bcf), "\\.(bcf|vcf\\.gz)$", "")
  String output_tsv = "${output_prefix}-tidy.tsv.zst"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools view --output-type u --include '(FILTER = "." || FILTER = "PASS") && INFO/SVLEN != "." && INFO/SVLEN > 0' '~{vcf_or_bcf}' \
      | bcftools view --output-type u --exclude 'INFO/SVTYPE = "BND"' \
      | bcftools query --include 'GT ~ "1"' --format '[%ID\t%ALT{0}\t%INFO/SVLEN\t%SAMPLE\n]' \
      | awk -F'\t' '{sub(/^</, "", $2); sub(/>$/, "", $2); print}' OFS='\t' \
      | zstd -q -c > '~{output_tsv}'
  >>>

  output {
    File tsv = output_tsv
  }
}

#=======================================================================
# Extract the SV clusters from a JoinRawCalls VCF or BCF.
#
# In GATK-SV, JoinRawCalls clusters the variants from ClusterBatch and
# reports one VCF record per cluster. The members of each cluster are
# listed in the INFO/MEMBERS field as a comma-separated list of IDs.
# Given a set of member IDs, an efficient way to find all the clusters
# containing these variants is to do a join operation. However, that
# requires each cluster member to be placed in its own row, which is
# what this task does.
#
# The clusters from the output of JoinRawCalls is reformatted into a
# TSV file of two columns:
# 1. Cluster ID (ID field of VCF)
# 2. Member ID (ID from INFO/MEMBERS), one row per cluster member
#
# Note that JoinRawCalls does not produce a BCF, but this task accepts
# one because the task GetContigFromVcf in this workflow outputs a BCF.
#
# Inputs
# ------
# vcf_or_bcf: Output of JoinRawCalls from GATK-SV, in VCF or BCF format.
#   The extension of this file is expected to be ".bcf" or ".vcf.gz".
#   Any other extension will result in a oddly named output file.
# base_docker: Path to Docker image.
#
# Outputs
# -------
# clusters: TSV file of clusters compressed with Zstandard.
#=======================================================================
task GetJoinRawCallsClusters {
  input {
    File vcf_or_bcf
    String base_docker
  }

  Int disk_size_gb = ceil(size(vcf_or_bcf, "GB") * 1.2) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  String output_prefix = sub(basename(vcf_or_bcf), "\\.(bcf|vcf\\.gz)$", "")
  String output_tsv = "${output_prefix}-sv_clusters.tsv.zst"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools query --format '%ID\t%INFO/MEMBERS\n' '~{vcf_or_bcf}' \
      | awk -F'\t' '$2 {split($2, a, /,/); for (i in a) print $1"\t"a[i]}' \
      | zstd -q -c > '~{output_tsv}'
  >>>

  output {
    File clusters = output_tsv
  }
}

task GatherSVs {
    input {
        Array[File] sv_tsvs
    }
    command <<<
        zstdcat ~{sep=' ' sv_tsvs} > all_svs.tsv
    >>>
    output {
        File out = "all_svs.tsv"
    }
}

task GatherClusters {
    input {
        Array[File] cluster_tsvs
    }
    command <<<
        zstdcat ~{sep=' ' cluster_tsvs} > all_clusters.tsv
    >>>
    output {
        File out = "all_clusters.tsv"
    }
}


task CountSVs {
  input {
    File svs_tsv
    Array[Array[String]] filters
    String pipeline_docker
  }

  Int disk_size_gb = ceil(size(svs_tsv, "GB") * 2) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 4
    disks: "local-disk ${disk_size_gb} SSD"
    docker: pipeline_docker
    maxRetries: 1
    memory: "8 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail
    
    python3 '/opt/outlier-exclusion/scripts/count_svs.py' \
        '~{svs_tsv}' \
        '~{write_tsv(filters)}' \
        sv_counts.tsv
  >>>

  output {
    File sv_counts_tsv = "sv_counts.tsv"
    File sv_filters_tsv = "~{write_tsv(filters)}"
  }
}

task DetermineOutlierSamples {
  input {
    File sv_counts_tsv
    File sv_filters_tsv
    File? wgd_scores
    Float min_wgd_score = -0.2
    Float max_wgd_score = 0.2
    Float iqr_multiplier

    String pipeline_docker
  }

  Float input_size = size([sv_counts_tsv, wgd_scores], "GB")
  Int disk_size_gb = ceil(input_size * 1.2) + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: pipeline_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 1
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    python3 '/opt/outlier-exclusion/scripts/determine_outlier_samples.py' \
      '~{sv_counts_tsv}' \
      '~{sv_filters_tsv}' \
      sv_count_outlier_samples.tsv \
      wgd_outlier_samples.tsv \
      '~{iqr_multiplier}' \
      ~{if defined(wgd_scores) then "--wgd-scores '" + wgd_scores + "'" else ""} \
      ~{if defined(wgd_scores) then "--min-wgd-score " + min_wgd_score else ""} \
      ~{if defined(wgd_scores) then "--max-wgd-score " + max_wgd_score else ""}
  >>>

  output {
    File sv_count_outlier_samples = "sv_count_outlier_samples.tsv"
    File wgd_outlier_samples = "wgd_outlier_samples.tsv"
  }
}

task ConvertBatchToTsv {
    input {
        Array[File] vcfs
    }
    command <<<
    bcftools query --include 'GT ~ "1" & INFO/SVTYPE != "BND"' \
      --format '[%ID\t%ALT{0}\t%INFO/SVLEN\t%SAMPLE\n]' \
      --vcf-list '~{write_lines(vcfs)}' \
        | awk -F'\t' '{sub(/^</, "", $2); sub(/>$/, "", $2); print}' OFS='\t' \
        | zstd -q -c > variants.tsv.zst
    >>>
    output {
        File tsv = "variants.tsv.zst"
    }
}


task DetermineOutlierVariants {
  input {
    File variants_tsv
    File outlier_samples_tsv
    File jrc_clusters_tsv
    Float min_outlier_sample_prop

    String pipeline_docker
  }


  Float input_size = size(variants_tsv, "GB") + size(jrc_clusters_tsv, "GB") + size(outlier_samples_tsv, "GB")
  Int disk_size_gb = ceil(input_size) + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: pipeline_docker
    maxRetries: 1
    memory: "2 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    python3 '/opt/outlier-exclusion/scripts/determine_outlier_variants.py' \
      '~{variants_tsv}' \
      '~{jrc_clusters_tsv}' \
      '~{outlier_samples_tsv}' \
      outlier_variants.list \
      '~{min_outlier_sample_prop}'
  >>>

  output {
    File outlier_variants = "outlier_variants.list"
  }
}


#=======================================================================
# Flag the outlier variants in the FilterGenotypes VCF.
#
# Inputs
# ------
# output_prefix: Prefix given to the output VCF.
# filter_genotypes_vcf: Output VCF from FilterGenotypes. A downstream
#   VCF is also permitted. The only requirement is that there must be
#   a "TRUTH_VID" key in the INFO field matching the ID in the VCF
#   from JoinRawCalls.
# outlier_variants: Files with lists of outlier variant IDs, one per
#   line.
# pipeline_docker: Path to Docker image.
#
# Outputs
# -------
# outlier_annotated_vcf: VCF with outlier variants flagged. The FILTER
#   field will have the "OUTLIER" flag if it is an outlier.
# outlier_annotated_vcf: VCF index for `outlier_annotated_vcf`.
#=======================================================================
task FlagOutlierVariants {
  input {
    String output_prefix
    File filter_genotypes_vcf
    Array[File] outlier_variants

    String pipeline_docker
  }

  Int disk_size_gb = ceil(size(filter_genotypes_vcf, "GB") * 2.0 + size(outlier_variants, "GB")) + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: pipeline_docker
    maxRetries: 1
    memory: "2 GiB"
    preemptible: 3
  }

  String output_vcf = "${output_prefix}-outlier_flagged.vcf.gz"
  String output_vcf_index = "${output_vcf}.tbi"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    cat '~{write_lines(outlier_variants)}' \
      | xargs cat \
      | LC_ALL=C sort -u > outlier_variants.list

    bgzip -cd '~{filter_genotypes_vcf}' \
      | gawk -f /opt/outlier-exclusion/scripts/flag_outliers.awk outlier_variants.list - \
      | bgzip -c > "~{output_vcf}"
    bcftools index --tbi "~{output_vcf}"
  >>>

  output {
    File outlier_annotated_vcf = output_vcf
    File outlier_annotated_vcf_index = output_vcf_index
  }
}
