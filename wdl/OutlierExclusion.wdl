version 1.0

workflow OutlierExclusion {
  input {

    File? optional

    # GetContigFromVcf --------------------------------------------------------
    File contigs_file

    File join_raw_calls_vcf
    File join_raw_calls_vcf_index
    String base_docker
    String pipeline_docker

    # MakeSvCountsDb ----------------------------------------------------------
    Array[Array[String]] svtypes_to_filter = [["DEL", 5000, 250000], ["DUP", 5000, 25000]]

    # DetermineOutlierSamples -------------------------------------------------
    # See DetermineOutlierSamples task for description of inputs.
    File? wgd_scores
    Float? min_wgd_score
    Float? max_wgd_score
    Float iqr_multiplier = 8.0

    # See FormatOutlierSamples
    # This takes precedence over outliers from WGD and SV counts.
    File? outlier_samples

    # DetermineOutlierVariants ------------------------------------------------
    Array[File]? clustered_depth_vcfs
    Array[File]? clustered_manta_vcfs
    Array[File]? clustered_wham_vcfs
    Array[File]? clustered_melt_vcfs
    Array[File]? clustered_scramble_vcfs
    Float min_outlier_sample_prop = 1.0

    # FlagOutlierVariants -----------------------------------------------------
    File filter_genotypes_vcf
    String output_prefix
  }

  Array[String] contigs = read_lines(contigs_file)

  scatter (contig in contigs) {
    call GetContigFromVcf {
      input:
        vcf = join_raw_calls_vcf,
        vcf_index = join_raw_calls_vcf_index,
        contig = contig,
        runtime_docker = base_docker
    }

    call GetJoinRawCallsClusters {
      input:
        vcf_or_bcf = GetContigFromVcf.contig_bcf,
        runtime_docker = base_docker
    }

    if (!defined(outlier_samples)) {
      call ConvertVcfOrBcfToTsv {
        input:
          vcf_or_bcf = GetContigFromVcf.contig_bcf,
          runtime_docker = base_docker
      }
    }
  }

  call MakeJoinRawCallsClustersDb {
    input:
      clusters = GetJoinRawCallsClusters.clusters,
      runtime_docker = base_docker
  }

  if (defined(ConvertVcfOrBcfToTsv.tsv)) {
    call MakeSvCountsDb {
      input:
        tsvs = select_all(ConvertVcfOrBcfToTsv.tsv),
        filters = svtypes_to_filter,
        runtime_docker = pipeline_docker
    }

    call DetermineOutlierSamples {
      input:
        sv_counts_db = MakeSvCountsDb.sv_counts_db,
        wgd_scores = wgd_scores,
        min_wgd_score = min_wgd_score,
        max_wgd_score = max_wgd_score,
        iqr_multiplier = iqr_multiplier,
        runtime_docker = pipeline_docker
    }
  }

  if (defined(outlier_samples)) {
    call FormatOutlierSamples {
      input:
        outlier_samples = select_first([outlier_samples]),
        runtime_docker = base_docker
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
        runtime_docker = pipeline_docker
    }
  }

  call FlagOutlierVariants {
    input:
      output_prefix = output_prefix,
      filter_genotypes_vcf = filter_genotypes_vcf,
      outlier_variants = DetermineOutlierVariants.outlier_variants,
      runtime_docker = pipeline_docker
  }

  output {
    File outlier_annotated_vcf = FlagOutlierVariants.outlier_annotated_vcf
    File outlier_annotated_vcf_index = FlagOutlierVariants.outlier_annotated_vcf_index
    File? sv_count_outlier_samples = DetermineOutlierSamples.outlier_samples
    File? wgd_outlier_samples = DetermineOutlierSamples.wgd_outlier_samples
  }
}

# Extract a single contig from a VCF and output as a BCF.
task GetContigFromVcf {
  input {
    File vcf
    File vcf_index
    String contig
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(vcf, "GB") * 1.2) + 16

  runtime {
    memory: "1 GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 8
  }

  String output_bcf = "${contig}.bcf"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools view --output-type u --regions '~{contig}' --write-index=csi
  >>>

  output {
    File contig_bcf = output_bcf
    File contig_bcf_index = "${output_bcf}.csi"
  }
}

# Reformat a VCF or BCF into a TSV that can be ingested by the SV counting task.
# The columns are:
# "variant ID" "SV type" "SV length" "sample ID"
# one row per carrier.
task ConvertVcfOrBcfToTsv {
  input {
    File vcf_or_bcf
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(vcf_or_bcf, "GB") * 2.0) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: runtime_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  String output_prefix = sub(basename(vcf_or_bcf), "\\.(bcf|vcf\\.gz)$", "")
  String output_tsv = "${output_prefix}-tidy.tsv.gz"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools view --output-type u --include '(FILTER = "." || FILTER = "PASS") && INFO/SVLEN != "." && INFO/SVLEN > 0' '~{vcf_or_bcf}' \
      | bcftools view --output-type u --exclude 'INFO/SVTYPE = "BND"' \
      | bcftools query --include 'GT ~ "1"' --format '[%ID\t%ALT{0}\t%INFO/SVLEN\t%SAMPLE\n]' \
      | awk -F'\t' '{sub(/^</, "", $2); sub(/>$/, "", $2); print}' OFS='\t' \
      | gzip -c > '~{output_tsv}'
  >>>

  output {
    File tsv = output_tsv
  }
}

# Extract the SV clusters from a JoinRawCalls VCF or BCF.
# The output is a tab-delimited file with the JoinedRawCalls variant ID in the
# first column and member variant ID (from ClusterBatch VCFs) in the second,
# one row per cluster member.
task GetJoinRawCallsClusters {
  input {
    File vcf_or_bcf
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(vcf_or_bcf, "GB") * 1.2) + 16

  runtime {
    memory: "1 GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 8
  }

  String output_prefix = sub(basename(vcf_or_bcf), "\\.(bcf|vcf\\.gz)$", "")
  String output_tsv = "${output_prefix}-sv_clusters.tsv.gz"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools query --format '%ID\t%INFO/MEMBERS\n' '~{vcf_or_bcf}' \
      | awk -F'\t' '$2 {split($2, a, /,/); for (i in a) print $1"\t"a[i]}' \
      | gzip -c > '~{output_tsv}'
  >>>

  output {
    File clusters = output_tsv
  }
}

# Make the database of JoinRawCalls variant clusters.
task MakeJoinRawCallsClustersDb {
  input {
    Array[File] clusters
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(clusters, "GB")) + 16

  runtime {
    memory: "1 GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 8
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    duckdb jrc_clusters.duckdb 'CREATE TABLE jrc_clusters (vid VARCHAR, member VARCHAR);'

    while read -r f; do
      duckdb jrc_clusters.duckdb "COPY jrc_clusters FROM '${f}' (FORMAT CSV, DELIMITER '\t', HEADER false);"
    done < '~{write_lines(clusters)}'
  >>>

  output {
    File jrc_clusters_db = "jrc_clusters.duckdb"
  }
}

# Make a database of SVs from a set of TSVs in the format output by ConvertVcfOrBcfToTsv
# and then make a database of SV counts per sample.
task MakeSvCountsDb {
  input {
    Array[File] tsvs
    Array[Array[String]] filters
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(tsvs, "GB") * 2) + 16

  runtime {
    memory: "8 GiB"
    disks: "local-disk ${disk_size_gb} SSD"
    cpus: 4
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 8
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

# Find outlier samples based on SV counts and optionally WGS score.
task DetermineOutlierSamples {
  input {
    File sv_counts_db
    # TSV of WGD scores of all samples. Sample ID in column 1, score in column 2.
    # Outliers found by WGD score will be used to call outlier variants across all SV types.
    File? wgd_scores
    # Samples with WGD less than this are outliers
    Float min_wgd_score = -0.2
    # Samples with WGD greater than this are outliers
    Float max_wgd_score = 0.2
    # Scaling factor for the SV counts per sample IQR. Samples with greater than
    # IQR * iqr_multiplier number of SVs are outliers.
    Float iqr_multiplier

    String runtime_docker
  }

  Float input_size = size([sv_counts_db, wgd_scores], "GB")
  Int disk_size_gb = ceil(input_size * 1.2) + 16

  runtime {
    memory: "1 GiB"
    cpu: 1
    bootDiskSizeGb: 8
    disks: "local-disk ${disk_size_gb} HDD"
    preemptible: 1
    maxRetries: 1
    docker: runtime_docker
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

# Load a user-provided TSV of outlier samples into a database.
# The format of the TSV is two columns, sample ID and SV type, no header.
# Each sample will only be used to define outlier variants in its matching
# variant class.
task FormatOutlierSamples {
  input {
    File outlier_samples
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(outlier_samples, "GB")) + 16

  runtime {
    memory: "1 GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 8
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

# Find variants in the ClusterBatch VCFs that are enriched for outlier samples.
task DetermineOutlierVariants {
  input {
    Array[File] clustered_vcfs
    File outlier_samples_db
    File jrc_clusters_db
    Float min_outlier_sample_prop

    String runtime_docker
  }


  Float input_size = size(clustered_vcfs, "GB") + size(jrc_clusters_db, "GB") + size(outlier_samples_db, "GB")
  Int disk_size_gb = ceil(input_size) + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: runtime_docker
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
        | gzip -c > variants.tsv.gz
    duckdb variants.duckdb "COPY variants FROM 'variants.tsv.gz' (HEADER false, DELIMITER '\t');"

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

# Add an OUTLIER filter flag to outlier enriched variants in FilterGenotypes VCF.
task FlagOutlierVariants {
  input {
    String output_prefix
    File filter_genotypes_vcf
    Array[File] outlier_variants

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(filter_genotypes_vcf, "GB") * 2.0 + size(outlier_variants, "GB")) + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${disk_size_gb} HDD"
    docker: runtime_docker
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

    while read -r f; do
      cat "${f}"
    done < '~{write_lines(outlier_variants)}' \
      | LC_ALL=C sort -u > outlier_variants.list

    bgzip -cd '~{filter_genotypes_vcf}' \
      | gawk -f /opt/outlier-exclusion/scripts/flag_outliers.awk outlier_variants.list - \
      | bgzip -c > "${output_vcf}"
    bcftools index --tbi "~{output_vcf}"
  >>>

  output {
    File outlier_annotated_vcf = output_vcf
    File outlier_annotated_vcf_index = output_vcf_index
  }
}
