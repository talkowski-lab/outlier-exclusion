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

  call MakeJoinRawCallsClustersDb {
    input:
      clusters = GetJoinRawCallsClusters.clusters,
      base_docker = base_docker
  }

  if (!defined(outlier_samples)) {
    call MakeSvCountsDb {
      input:
        tsvs = select_all(ConvertVcfOrBcfToTsv.tsv),
        filters = svtypes_to_filter,
        pipeline_docker = pipeline_docker
    }

    call DetermineOutlierSamples {
      input:
        sv_counts_db = MakeSvCountsDb.sv_counts_db,
        wgd_scores = wgd_scores,
        min_wgd_score = min_wgd_score,
        max_wgd_score = max_wgd_score,
        iqr_multiplier = iqr_multiplier,
        pipeline_docker = pipeline_docker
    }
  }

  if (defined(outlier_samples)) {
    call FormatOutlierSamples {
      input:
        outlier_samples = select_first([outlier_samples]),
        base_docker = base_docker
    }
  }

  File outlier_samples_db = select_first([FormatOutlierSamples.db, DetermineOutlierSamples.sv_counts_db_with_outliers])
  Array[Array[File]] clustered_vcfs = select_all([clustered_depth_vcfs,
    clustered_manta_vcfs,
    clustered_wham_vcfs,
    clustered_melt_vcfs,
    clustered_scramble_vcfs])

  scatter (batch_vcfs in transpose(clustered_vcfs)) {
    call DetermineOutlierVariants {
      input:
        clustered_vcfs = batch_vcfs,
        outlier_samples_db = outlier_samples_db,
        min_outlier_sample_prop = min_outlier_sample_prop,
        jrc_clusters_db = MakeJoinRawCallsClustersDb.jrc_clusters_db,
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
    File? sv_count_outlier_samples = DetermineOutlierSamples.outlier_samples
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

    bcftools view --output-type b --regions '~{contig}' --write-index=csi
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

#=======================================================================
# Load SV clusters from JoinRawCalls into a DuckDB database.
#
# Inputs
# ------
# clusters: TSV files as produced by GetJoinRawCallsClusters.
# base_docker: Path to Docker image.
#
# Outputs
# -------
# jrc_clusters_db: DuckDB database with clusters loaded into a table
#   named 'jrc_clusters'.
#=======================================================================
task MakeJoinRawCallsClustersDb {
  input {
    Array[File] clusters
    String base_docker
  }

  Int disk_size_gb = ceil(size(clusters, "GB")) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    duckdb jrc_clusters.duckdb 'CREATE TABLE jrc_clusters (vid VARCHAR, member VARCHAR);'
    gawk '$0 {print "COPY jrc_clusters FROM \047"$0"\047 (FORMAT CSV, DELIMITER \047\\t\047, HEADER false);"}' \
      '~{write_lines(clusters)}' \
      | duckdb jrc_clusters.duckdb
  >>>

  output {
    File jrc_clusters_db = "jrc_clusters.duckdb"
  }
}

#=======================================================================
# Count the number of SVs per sample for the given filters.
#
# Inputs
# ------
# tsvs: TSV files as output by ConvertVcfOrBcfToTsv.
# filters: Filters for counting SVs. Each filter is a three-element
#   Array of the form `[SVTYPE, MIN_SVLEN, MAX_SVLEN]` and the SV counts
#   per sample will be computed for each filter. Note that the filters
#   are not checked for overlap so if there is overlap between two
#   filters, a variant could be counted twice. On Terra, all elements of
#   the filters may have to given as Strings because type coercion may
#   not work. Negative and positive infinity may be expressed using
#   "-Inf" and "Inf", respectively.
# pipeline_docker: Path to Docker image.
#
#
# Outputs
# -------
# sv_counts_db: DuckDB database with counts of SVs per sample. There
#   will be one table for the filters and one table of counts per
#   filter.
#=======================================================================
task MakeSvCountsDb {
  input {
    Array[File] tsvs
    Array[Array[String]] filters
    String pipeline_docker
  }

  Int disk_size_gb = ceil(size(tsvs, "GB") * 2) + 16

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

    duckdb svs.duckdb << 'EOF'
    CREATE TABLE svs (
        vid VARCHAR,
        svtype VARCHAR,
        svlen INTEGER,
        sample VARCHAR
    );
    EOF
    gawk '{print "COPY svs FROM \047"$0"\047 (FORMAT CSV, DELIMITER \047\\t\047, HEADER false);"}' \
      '~{write_lines(tsvs)}' \
      | duckdb svs.duckdb

    duckdb sv_counts.duckdb << 'EOF'
    CREATE TABLE sv_filters (
        svtype VARCHAR,
        min_svlen DOUBLE,
        max_svlen DOUBLE
    );
    COPY sv_filters
    FROM '~{write_tsv(filters)}' (
        FORMAT CSV,
        DELIMITER '\t',
        HEADER false
    );
    CREATE SEQUENCE id_sequence START 1;
    ALTER TABLE sv_filters ADD COLUMN id INTEGER DEFAULT nextval('id_sequence');

    UPDATE sv_filters
    SET min_svlen = trunc(min_svlen)
    WHERE isfinite(min_svlen);
    UPDATE sv_filters
    SET max_svlen = trunc(max_svlen)
    WHERE isfinite(max_svlen);
    EOF

    python3 '/opt/outlier-exclusion/scripts/count_svs.py' sv_counts.duckdb svs.duckdb
  >>>

  output {
    File sv_counts_db = "sv_counts.duckdb"
  }
}

#=======================================================================
# Find outlier samples based on SV counts and optionally WGD scores.
#
# Finding outlier samples based on SV counts is done by first computing
# the interquartile range (IQR) of SVs per sample. Then the IQR is
# multiplied by a scaling factor to determine the upper limit of SVs per
# sample. Any sample with more SVs than the limit is marked as an
# outlier. Thresholds and outliers are determined independently for each
# counts filter.
#
# Finding outlier samples based on WGD scores is done with simple
# thresholding. Any sample with a score less than the specified
# minimum or greater than the specified maximum is an outlier.
#
# Inputs
# ------
# sv_counts_db: DuckDB database with counts of SVs per sample.
# wgd_scores: Two-column TSV file of WGD scores.
#   1. Sample ID
#   2. WGS score
# min_wgd_score: Minimum WGD score. Samples with scores less than this
#   are outliers.
# max_wgd_score: Maximum WGD score. Samples with scores greater than
#   this are outliers.
# iqr_multiplier: Scaling factor for IQR threshold. Samples with more
#   then `iqr_multiplier * IQR(SVs per sample)` SVs are outliers.
# pipeline_docker: Path to Docker image.
#
# Outputs
# -------
# sv_counts_db_with_outliers: DuckDB database with tables for outlier
#   samples. The tables are added to `sv_counts_db`.
# outlier_samples: TSV file listing all the outlier samples found by
#   counting SVs per sample and their associated filters.
# wgd_outlier_samples: TSV file listing all the outlier samples found
#   by thresholding on WGD scores.
#=======================================================================
task DetermineOutlierSamples {
  input {
    File sv_counts_db
    File? wgd_scores
    Float min_wgd_score = -0.2
    Float max_wgd_score = 0.2
    Float iqr_multiplier

    String pipeline_docker
  }

  Float input_size = size([sv_counts_db, wgd_scores], "GB")
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

    mv '~{sv_counts_db}' sv_counts_with_outliers.duckdb
    python3 '/opt/outlier-exclusion/scripts/determine_outlier_samples.py' \
      sv_counts_with_outliers.duckdb \
      '~{iqr_multiplier}' \
      ~{if defined(wgd_scores) then "--wgd-scores '" + wgd_scores + "'" else ""} \
      ~{if defined(wgd_scores) then "--min-wgd-score " + min_wgd_score else ""} \
      ~{if defined(wgd_scores) then "--max-wgd-score " + max_wgd_score else ""}

    python3 '/opt/outlier-exclusion/scripts/dump_outlier_samples.py' \
      sv_counts_with_outliers.duckdb \
      wgd_outlier_samples.tsv \
      dump

    printf 'sample_id\tcount\tsvtype\tmin_svlen\tmax_svlen\n' > sv_count_outlier_samples.tsv
    find dump -type f -name '*.tsv' -exec cat '{}' \; >> sv_count_outlier_samples.tsv
  >>>

  output {
    File sv_counts_db_with_outliers = "sv_counts_with_outliers.duckdb"
    File outlier_samples = "sv_count_outlier_samples.tsv"
    File wgd_outlier_samples = "wgd_outlier_samples.tsv"
  }
}

#=======================================================================
# Load a user-provided TSV file of outlier samples into a database.
#
# Inputs
# ------
# outlier_samples: A two-column TSV file of outlier samples.
#   1. Sample ID
#   2. SV type
# base_docker: Path to Docker image.
#
# Outputs
# -------
# db: DuckDB database.
#=======================================================================
task FormatOutlierSamples {
  input {
    File outlier_samples
    String base_docker
  }

  Int disk_size_gb = ceil(size(outlier_samples, "GB")) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  command <<<
    duckdb outlier_samples.duckdb << 'EOF'
    CREATE TABLE outlier_samples (
        sample VARCHAR,
        svtype VARCHAR
    );
    COPY outlier_samples
    FROM '~{outlier_samples}' (
        FORMAT CSV,
        DELIMITER '\t',
        HEADER false
    );
    EOF
  >>>

  output {
    File db = 'outlier_samples.duckdb'
  }
}

#=======================================================================
# Find outlier enriched variants in the ClusterBatch VCFs.
#
# For each variant, the carrier samples (those with a "1" in GT) are
# collected, the proprotion of carriers that are outlier samples is
# computed, and then variants with a proportion of outlier samples
# greater than the given threshold are labeled as outlier variants.
#
# Outlier samples are defined in different ways depending on the
# workflow inputs. Outlier samples found by counting the number of SVs
# per sample are only considered to be outliers for variants matching
# their filter SV type. For example, if a sample was determined to be an
# outlier from counting DELs in the range 5-25Kb, it would only count
# as an outlier in DEL ClusterBatch variants.
#
# Outlier samples found from WGD scores are considered outlier samples
# in all SV types. Outlier variants can only be found in this way if
# WGD scores were provided in DetermineOutlierSamples.
#
# Outlier sample proportions using SV count outliers and WGD outliers
# are computed independently, but the lists of outlier variants found
# from each are merged into one.
#
# Outlier samples given as a custom input must have an accompanying SV
# type and are matched in that way. If these are given as inputs, the
# outlier samples from SV counts and WGD scores will not be used.
#
# All the outlier variants in the ClusterBatch VCFs are matched to
# their JoinRawCalls IDs using the information from the JoinRawCalls
# clusters.
#
# Inputs
# ------
# clustered_vcfs: The VCFs from GATK-SV ClusterBatch. When running the
#   OutlierExclusion workflow, these VCFs will be all the raw algorithm
#   VCFs for a single batch, but the task will accept VCFs from mixed
#   batches.
# outlier_samples_db: DuckDB database containing the outlier samples
#   (output from DetermineOutlierSamples).
# jrc_clusters_db: DuckDB database containing the JoinRawCalls clusters
#   (output from MakeJoinRawCallsClustersDb).
# min_outlier_sample_prop: Minimum proportion of outlier samples a
#   variant must have to be considered an outlier.
# pipeline_docker: Path to Docker image.
#
# Outputs
# -------
# outlier_variants: List of IDs of the outlier variants. These IDs are
#   the ones in the JoinRawCalls VCF.
#=======================================================================
task DetermineOutlierVariants {
  input {
    Array[File] clustered_vcfs
    File outlier_samples_db
    File jrc_clusters_db
    Float min_outlier_sample_prop

    String pipeline_docker
  }


  Float input_size = size(clustered_vcfs, "GB") + size(jrc_clusters_db, "GB") + size(outlier_samples_db, "GB")
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

    bcftools query --include 'GT ~ "1" & INFO/SVTYPE != "BND"' \
      --format '[%ID\t%ALT{0}\t%INFO/SVLEN\t%SAMPLE\n]' \
      --vcf-list '~{write_lines(clustered_vcfs)}' \
        | awk -F'\t' '{sub(/^</, "", $2); sub(/>$/, "", $2); print}' OFS='\t' \
        | zstd -q -c > variants.tsv.zst
    duckdb variants.duckdb 'CREATE TABLE variants (vid VARCHAR, svtype VARCHAR, svlen INTEGER, sample VARCHAR);'
    duckdb variants.duckdb "COPY variants FROM 'variants.tsv.zst' (HEADER false, DELIMITER '\t');"

    python3 '/opt/outlier-exclusion/scripts/determine_outlier_variants.py' \
      '~{outlier_samples_db}' \
      '~{jrc_clusters_db}' \
      variants.duckdb \
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
