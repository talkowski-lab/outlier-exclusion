version 1.0

workflow OutlierExclusion {
  input {

    File? optional

    # GetContigFromVcf --------------------------------------------------------
    File contigs_file

    File join_raw_calls_vcf
    File join_raw_calls_vcf_index
    String docker

    # MakeSvCountsDb ----------------------------------------------------------
    Array[Array[String]] svtypes_to_filter = [["DEL", 5000, 250000], ["DUP", 5000, 25000]]

    # DetermineOutlierSamples -------------------------------------------------
    # See DetermineOutlierSamples task for description of inputs.
    File? wgd_scores
    Float? min_wgd_score
    Float? max_wgd_score
    Float iqr_multiplier = 8.0
    File? outlier_samples

    # DetermineOutlierVariants ------------------------------------------------
    Array[File]? clustered_depth_vcfs
    Array[File]? clustered_manta_vcfs
    Array[File]? clustered_wham_vcfs
    Array[File]? clustered_melt_vcfs
    Array[File]? clustered_scramble_vcfs
    Array[File]? clustered_depth_vcf_indicies
    Array[File]? clustered_manta_vcf_indicies
    Array[File]? clustered_wham_vcf_indicies
    Array[File]? clustered_melt_vcf_indicies
    Array[File]? clustered_scramble_vcf_indicies
    Float min_outlier_sample_prop = 1.0

    # FlagOutlierVariants -----------------------------------------------------
    File filtered_vcf
    File filtered_vcf_index
    String cohort_prefix
  }

  Array[String] contigs = read_lines(contigs_file)

  scatter (contig in contigs) {
    call GetContigFromVcf {
      input:
        vcf = join_raw_calls_vcf,
        vcf_index = join_raw_calls_vcf_index,
        contig = contig,
        runtime_docker = docker
    }

    call GetJoinRawCallsClusters {
      input:
        bcf = GetContigFromVcf.contig_bcf,
        runtime_docker = docker
    }

    if (!defined(outlier_samples)) {
      call ConvertBcfToTsv {
        input:
          bcf = GetContigFromVcf.contig_bcf,
          bcf_index = GetContigFromVcf.contig_bcf_index,
          runtime_docker = docker
      }
    }
  }

  call MakeJoinRawCallsClustersDb {
    input:
      clusters = GetJoinRawCallsClusters.clusters,
      runtime_docker = docker
  }

  if (defined(ConvertBcfToTsv.tsv)) {
    call MakeSvDb {
      input:
        tsvs = select_all(ConvertBcfToTsv.tsv),
        runtime_docker = docker
    }

    call MakeSvCountsDb {
      input:
        sv_db = MakeSvDb.sv_db,
        filters = svtypes_to_filter,
        runtime_docker = docker
    }

    call DetermineOutlierSamples {
      input:
        sv_counts_db = MakeSvCountsDb.sv_counts_db,
        wgd_scores = wgd_scores,
        min_wgd_score = min_wgd_score,
        max_wgd_score = max_wgd_score,
        iqr_multiplier = iqr_multiplier,
        runtime_docker = docker
    }
  }

  if (defined(outlier_samples)) {
    call FormatOutlierSamples {
      input:
        outlier_samples = select_first([outlier_samples]),
        runtime_docker = docker
    }
  }

  File outlier_samples_db = select_first([FormatOutlierSamples.db, DetermineOutlierSamples.sv_counts_db_with_outliers])
  Array[File] looper = select_first([clustered_depth_vcfs, clustered_manta_vcfs,
    clustered_wham_vcfs, clustered_melt_vcfs, clustered_scramble_vcfs])

  Array[File] depths = select_first([clustered_depth_vcfs, []])
  Array[File] mantas = select_first([clustered_manta_vcfs, []])
  Array[File] whams = select_first([clustered_wham_vcfs, []])
  Array[File] melts = select_first([clustered_melt_vcfs, []])
  Array[File] scrambles = select_first([clustered_scramble_vcfs, []])
  Array[File] depth_idxs = select_first([clustered_depth_vcf_indicies, []])
  Array[File] manta_idxs = select_first([clustered_manta_vcf_indicies, []])
  Array[File] wham_idxs = select_first([clustered_wham_vcf_indicies, []])
  Array[File] melt_idxs = select_first([clustered_melt_vcf_indicies, []])
  Array[File] scramble_idxs = select_first([clustered_scramble_vcf_indicies, []])

  scatter (i in range(length(looper))) {
    call DetermineOutlierVariants {
      input:
        clustered_depth_vcf = if defined(clustered_depth_vcfs) then depths[i] else optional,
        clustered_manta_vcf = if defined(clustered_manta_vcfs) then mantas[i] else optional,
        clustered_wham_vcf = if defined(clustered_wham_vcfs) then whams[i] else optional,
        clustered_melt_vcf = if defined(clustered_melt_vcfs) then melts[i] else optional,
        clustered_scramble_vcf = if defined(clustered_scramble_vcfs) then scrambles[i] else optional,
        clustered_depth_vcf_index = if defined(clustered_depth_vcf_indicies) then depth_idxs[i] else optional,
        clustered_manta_vcf_index = if defined(clustered_manta_vcf_indicies) then manta_idxs[i] else optional,
        clustered_wham_vcf_index = if defined(clustered_wham_vcf_indicies) then wham_idxs[i] else optional,
        clustered_melt_vcf_index = if defined(clustered_melt_vcf_indicies) then melt_idxs[i] else optional,
        clustered_scramble_vcf_index = if defined(clustered_scramble_vcf_indicies) then scramble_idxs[i] else optional,
        outlier_samples_db = outlier_samples_db,
        min_outlier_sample_prop = min_outlier_sample_prop,
        jrc_clusters_db = MakeJoinRawCallsClustersDb.jrc_clusters_db,
        runtime_docker = docker
    }
  }

  call FlagOutlierVariants {
    input:
      cohort_prefix = cohort_prefix,
      filtered_vcf = filtered_vcf,
      filtered_vcf_index = filtered_vcf_index,
      outlier_variants = DetermineOutlierVariants.outlier_variants,
      runtime_docker = docker
  }

  output {
    File outlier_annotated_vcf = FlagOutlierVariants.outlier_annotated_vcf
    File outlier_annotated_vcf_index = FlagOutlierVariants.outlier_annotated_vcf_index
    File? determined_outlier_samples = DetermineOutlierSamples.outlier_samples
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
    bootDiskSizeGb: 16
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

# Reformat a BCF into a TSV that can be ingested by the SV counting task.
# The columns are:
# "variant ID" "SV type" "SV length" "sample ID"
# one row per carrier.
task ConvertBcfToTsv {
  input {
    File bcf
    File bcf_index
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(bcf, 'GB') * 2.0) + 16

  runtime {
    memory: "1 GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  String output_prefix = basename(bcf, ".bcf")
  String output_tsv = "${output_prefix}-tidy.tsv.gz"

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools view --output-type u --include '(FILTER = "." || FILTER = "PASS") && INFO/SVLEN != "." && INFO/SVLEN > 0' '~{bcf}' \
      | bcftools view --output-type u --exclude 'INFO/SVTYPE = "BND"' \
      | bcftools query --include 'GT ~ "1"' --format '[%ID\t%ALT{0}\t%INFO/SVLEN\t%SAMPLE\n]' \
      | awk -F'\t' '{sub(/^</, "", $2); sub(/>$/, "", $2); print}' OFS='\t' \
      | gzip -c > '~{output_tsv}'
  >>>

  output {
    File tsv = output_tsv
  }
}

# Extract the SV clusters from a JoinRawCalls BCF.
# The output is a tab-delimited file with the JoinedRawCalls variant ID in the
# first column and member variant ID (from ClusterBatch VCFs) in the second,
# one row per cluster member.
task GetJoinRawCallsClusters {
  input {
    File bcf
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(bcf, "GB") * 1.2) + 16

  runtime {
    memory: "1 GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  String bcf_prefix = basename(bcf, ".bcf")

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools query --format '%ID\t%INFO/MEMBERS\n' '~{bcf}' \
      | awk -F'\t' '$2 {split($2, a, /,/); for (i in a) print $1"\t"a[i]}' \
      | gzip -c > '~{bcf_prefix}-sv_clusters.tsv.gz'
  >>>

  output {
    File clusters = "${bcf_prefix}-sv_clusters.tsv.gz"
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
    bootDiskSizeGb: 16
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

# Make a database of SVs from a set of TSVs in the format output by ConvertBcfToTsv.
task MakeSvDb {
  input {
    Array[File] tsvs
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(tsvs, "GB") * 1.2) + 16

  runtime {
    memory: "1 GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
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

    while read -r f; do
      duckdb svs.duckdb "COPY svs FROM '${f}' (FORMAT CSV, DELIMITER '\t', HEADER false);"
    done < '~{write_lines(tsvs)}'
  >>>

  output {
    File sv_db = "svs.duckdb"
  }
}

# Make a database of SV counts given the database from MakeSvDb and a set of filters
# which determine the SV types and size ranges to summarize.
task MakeSvCountsDb {
  input {
    File sv_db
    Array[Array[String]] filters
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(sv_db, "GB") * 1.2) + 16

  runtime {
    memory: "2 GiB"
    disks: "local-disk ${disk_size_gb} HDD"
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

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

    python3 '/opt/outlier-exclusion/scripts/count_svs.py' sv_counts.duckdb '~{sv_db}'
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

  Float input_size = size([sv_counts_db, wgd_scores], 'GB')
  Int disk_size_gb = ceil(input_size * 1.2) + 16

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
      wgd_outliers.tsv \
      dump

    printf 'sample_id\tcount\tsvtype\tmin_svlen\tmax_svlen\n' > sv_count_outlier_samples.tsv
    find dump -type f -name '*.tsv' -exec cat '{}' \; >> sv_count_outlier_samples.tsv
  >>>

  runtime {
    memory: "1 GiB"
    cpu: 1
    bootDiskSizeGb: 16
    disks: "local-disk ${disk_size_gb} HDD"
    preemptible: 1
    maxRetries: 1
    docker: runtime_docker
  }

  output {
    File sv_counts_db_with_outliers = "sv_counts_with_outliers.duckdb"
    File outlier_samples = "sv_count_outlier_samples.tsv"
    File wgd_outlier_samples = "wgd_outlier_samples.tsv"
  }
}

task FormatOutlierSamples {
  input {
    File outlier_samples
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(outlier_samples, 'GB')) + 16

  runtime {
    memory: '1GiB'
    disks: 'local-disk ${disk_size_gb} HDD'
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
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

task DetermineOutlierVariants {
  input {
    File? clustered_depth_vcf
    File? clustered_manta_vcf
    File? clustered_wham_vcf
    File? clustered_melt_vcf
    File? clustered_scramble_vcf
    File? clustered_depth_vcf_index
    File? clustered_manta_vcf_index
    File? clustered_wham_vcf_index
    File? clustered_melt_vcf_index
    File? clustered_scramble_vcf_index
    File outlier_samples_db
    File jrc_clusters_db
    Float min_outlier_sample_prop

    String runtime_docker
  }

  Array[File] clusterbatch_vcfs = select_all([
    clustered_depth_vcf, clustered_manta_vcf, clustered_wham_vcf,
    clustered_melt_vcf
  ])
  Array[File] clusterbatch_vcf_indicies = select_all([
    clustered_depth_vcf_index, clustered_manta_vcf_index, clustered_wham_vcf_index,
    clustered_melt_vcf_index
  ])

  Int disk_size_gb = ceil(
    size(clusterbatch_vcfs, 'GB') * 2.0
    + size(jrc_clusters_db, 'GB')
    + size(outlier_samples_db, 'GB')
  ) + 16

  runtime {
    memory: '2GiB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    mkdir clusterbatch_vcfs
    cat '~{write_lines(clusterbatch_vcfs)}' '~{write_lines(clusterbatch_vcf_indicies)}' | while read -r vcf; do
      mv -t clusterbatch_vcfs "${vcf}" 
    done

    find clusterbatch_vcfs -name '*.vcf.gz' -print \
      | awk -F'/' '{bn=$NF; sub(/\.cluster_batch\.(depth|wham|manta|melt|scramble)\.vcf.gz$/, "", bn); print bn}' \
      | sort -u > clusterbatch_ids.list

    mkdir clusterbatch_dbs
    cat > reformat.bash << 'EOF'
    set -o errexit
    set -o nounset
    set -o pipefail
    query_vcf() {
      bcftools query \
        --include 'GT ~ "1" & INFO/SVTYPE != "BND"' \
        --format '[%ID\t%ALT{0}\t%INFO/SVLEN\t%SAMPLE\n]' "$1" \
        | awk -F'\t' '{sub(/^</, "", $2); sub(/>$/, "", $2); print}' OFS='\t' >> "$2"
    }
    tsv="clusterbatch_dbs/${1}_variants.tsv"
    db="clusterbatch_dbs/${1}_variants.duckdb"
    : > "${tsv}"
    rm -f -- "${db}"
    duckdb "${db}" 'CREATE TABLE variants (vid VARCHAR, svtype VARCHAR, svlen INTEGER, sample VARCHAR);'
    depth="clusterbatch_vcfs/${1}.cluster_batch.depth.vcf.gz"
    wham="clusterbatch_vcfs/${1}.cluster_batch.wham.vcf.gz"
    manta="clusterbatch_vcfs/${1}.cluster_batch.manta.vcf.gz"
    melt="clusterbatch_vcfs/${1}.cluster_batch.melt.vcf.gz"
    scramble="clusterbatch_vcfs/${1}.cluster_batch.scramble.vcf.gz"
    test -r "${depth}" && query_vcf "${depth}" "${tsv}"
    test -r "${wham}" && query_vcf "${wham}" "${tsv}"
    test -r "${manta}" && query_vcf "${manta}" "${tsv}"
    test -r "${melt}" && query_vcf "${melt}" "${tsv}"
    test -r "${scramble}" && query_vcf "${scramble}" "${tsv}"
    duckdb "${db}" "COPY variants FROM '${tsv}' (FORMAT CSV, HEADER false, DELIMITER '\t');"
    EOF
    xargs -L 1 -P 0 -- bash reformat.bash < clusterbatch_ids.list

    python3 '/opt/outlier-exclusion/scripts/determine_outlier_variants.py' \
      '~{outlier_samples_db}' \
      '~{jrc_clusters_db}' \
      clusterbatch_dbs \
      '~{min_outlier_sample_prop}' \
      | LC_ALL=C sort -u > outlier_variants.list
  >>>

  output {
    File outlier_variants = 'outlier_variants.list'
  }
}

task FlagOutlierVariants {
  input {
    String cohort_prefix
    File filtered_vcf
    File filtered_vcf_index
    Array[File] outlier_variants

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(filtered_vcf, 'GB') * 2.0) + 16

  runtime {
    memory: '2GiB'
    cpu: 4
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    while read -r f; do
      cat "${f}"
    done < '~{write_lines(outlier_variants)}' \
      | LC_ALL=C sort -u > outlier_variants.list

    bcftools query --include 'INFO/TRUTH_VID != ""' \
      --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/TRUTH_VID\n' \
      '~{filtered_vcf}' \
      | LC_ALL=C sort -k6,6 > filtered_vcf_variants.tsv
    LC_ALL=C join -1 6 -2 1 -o 1.1,1.2,1.3,1.4,1.5 -t $'\t' \
      filtered_vcf_variants.tsv \
      outlier_variants.list > filtered_calls_outliers.tsv

    awk -F'\t' '{print $0"\toutlier"}' 'filtered_calls_outliers.tsv' \
      | sort -k1,1 -k2,2n > annotations.tsv
    bgzip annotations.tsv
    tabix --begin 2 --end 2 --sequence 1 annotations.tsv.gz

    printf '##FILTER=<ID=outlier,Description="Variant enriched by outlier samples">' \
      > header.txt

    bcftools annotate \
      --annotations annotations.tsv.gz \
      --columns 'CHROM,POS,REF,ALT,~ID,.FILTER' \
      --header-lines header.txt \
      --output '~{cohort_prefix}-outliers_annotated.vcf.gz' \
      --output-type z \
      --threads 4 \
      --write-index=tbi \
      '~{filtered_vcf}'
  >>>

  output {
    File outlier_annotated_vcf = '~{cohort_prefix}-outliers_annotated.vcf.gz'
    File outlier_annotated_vcf_index = '~{cohort_prefix}-outliers_annotated.vcf.gz.tbi'
  }
}
