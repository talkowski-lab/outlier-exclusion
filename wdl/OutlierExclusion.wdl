version 1.0

workflow OutlierExclusion {
  input {

    # GetContigsArray ---------------------------------------------------------
    File joined_raw_calls_vcf
    File joined_raw_calls_vcf_index
    String docker

    # MakeTidyVCF -------------------------------------------------------------
    Array[String] svtypes_to_filter = ['DEL;0;inf', 'DUP;0;inf']

    # DetermineOutlierSamples -------------------------------------------------
    File wgd_scores
    Float min_wgd_score = -0.2
    Float max_wgd_score = 0.2
    Float iqr_multiplier = 8.0

    # DetermineOutlierVariants ------------------------------------------------
    Array[File] clustered_depth_vcfs
    Array[File] clustered_manta_vcfs
    Array[File] clustered_wham_vcfs
    Array[File] clustered_melt_vcfs
    Array[File] clustered_depth_vcf_indicies
    Array[File] clustered_manta_vcf_indicies
    Array[File] clustered_wham_vcf_indicies
    Array[File] clustered_melt_vcf_indicies
    File filtered_vcf
    File filtered_vcf_index
    Float min_outlier_sample_prop = 1.0

    # FlagOutlierVariants -----------------------------------------------------
    String cohort_prefix
  }

  call GetContigsArray  {
    input:
      vcf = joined_raw_calls_vcf,
      vcf_index = joined_raw_calls_vcf_index,
      runtime_docker = docker
  }

  scatter (contig in GetContigsArray.contigs) {
    call MakeTidyVCF {
      input:
        vcf = joined_raw_calls_vcf,
        vcf_index = joined_raw_calls_vcf_index,
        contig = contig,
        filters = svtypes_to_filter,
        runtime_docker = docker
    }

    call GetJoinedRawCallsClusters {
      input:
        vcf = joined_raw_calls_vcf,
        vcf_index = joined_raw_calls_vcf_index,
        contig = contig,
        runtime_docker = docker
    }
  }

  call MakeJoinedRawCallsDB {
    input:
      tidy_vcfs = MakeTidyVCF.tidy_vcf,
      clusters = GetJoinedRawCallsClusters.clusters,
      runtime_docker = docker
  }

  call CountSVsPerGenome {
    input:
      joined_raw_calls_db = MakeJoinedRawCallsDB.joined_raw_calls_db,
      svtypes_to_filter = svtypes_to_filter,
      runtime_docker = docker
  }

  call DetermineOutlierSamples {
    input:
      sv_counts_db = CountSVsPerGenome.sv_counts_db,
      wgd_scores = wgd_scores,
      min_wgd_score = min_wgd_score,
      max_wgd_score = max_wgd_score,
      iqr_multiplier = iqr_multiplier,
      runtime_docker = docker
  }

  scatter (i in range(length(clustered_depth_vcfs))) {
    call DetermineOutlierVariants {
      input:
        clustered_depth_vcfs = [clustered_depth_vcfs[i]],
        clustered_manta_vcfs = [clustered_manta_vcfs[i]],
        clustered_wham_vcfs = [clustered_wham_vcfs[i]],
        clustered_melt_vcfs = [clustered_melt_vcfs[i]],
        clustered_depth_vcf_indicies = [clustered_depth_vcf_indicies[i]],
        clustered_manta_vcf_indicies = [clustered_manta_vcf_indicies[i]],
        clustered_wham_vcf_indicies = [clustered_wham_vcf_indicies[i]],
        clustered_melt_vcf_indicies = [clustered_melt_vcf_indicies[i]],
        sv_counts_db = DetermineOutlierSamples.sv_counts_db_with_outliers,
        min_outlier_sample_prop = min_outlier_sample_prop,
        joined_raw_calls_db = MakeJoinedRawCallsDB.joined_raw_calls_db,
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
    File outlier_samples = DetermineOutlierSamples.outlier_samples
  }
}

task GetContigsArray {
  input {
    File vcf
    File vcf_index

    String runtime_docker
  }

  Int disk_size_gib = ceil(size(vcf, 'GiB')) + 8

  runtime {
    memory: '256MiB'
    disks: 'local-disk ${disk_size_gib} HDD'
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }
  
  command <<<
    bcftools index --stats '~{vcf}' \
      | cut -f 1 > contigs.list
  >>>

  output {
    Array[String] contigs = read_lines('contigs.list')
  }
}

task MakeTidyVCF {
  input {
    File vcf
    File vcf_index
    String contig
    Array[String] filters
    String runtime_docker
  }

  Int disk_size_gib = ceil(size(vcf, 'GiB') * 2.0) + 8

  runtime {
    memory: '1GiB'
    disks: 'local-disk ${disk_size_gib} HDD'
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

    bcftools view --output-type u --regions '~{contig}' '~{vcf}' \
      | bcftools query --include 'INFO/SVTYPE != "BND" & GT ~ "1"' \
          --format '[%ID\t%ALT{0}\t%INFO/SVLEN\t%SAMPLE\n]' \
      | awk -F'\t' '{sub(/^</, "", $2); sub(/>$/, "", $2); print}' OFS='\t' \
      | gzip -c > '~{contig}_tidy_vcf.tsv.gz'
  >>>

  output {
    File tidy_vcf = '${contig}_tidy_vcf.tsv.gz'
  }
}

task GetJoinedRawCallsClusters {
  input {
    File vcf
    File vcf_index 
    String contig
    String runtime_docker
  }

  Int disk_size_gib = ceil(size(vcf, 'GiB') * 2.0) + 8

  runtime {
    memory: '512MiB'
    disks: 'local-disk ${disk_size_gib} HDD'
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

    bcftools view --output-type u --regions '~{contig}' '~{vcf}' \
      | bcftools query --format '%ID\t%INFO/MEMBERS\n' \
      | awk -F'\t' '$2 {split($2, a, /,/); for (i in a) print $1"\t"a[i]}' \
      | gzip -c > '~{contig}_sv_clusters.tsv.gz'
  >>>

  output {
    File clusters = '${contig}_sv_clusters.tsv.gz'
  }
}

task MakeJoinedRawCallsDB {
  input {
    Array[File] tidy_vcfs
    Array[File] clusters

    String runtime_docker
  }

  Int disk_size_gib = ceil(size(flatten([tidy_vcfs, clusters]), 'GiB') * 3.0) + 8

  runtime {
    memory: '4GiB'
    disks: 'local-disk ${disk_size_gib} HDD'
    cpus: 4
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    duckdb joined_raw_calls.duckdb << 'EOF'
    CREATE TABLE joined_raw_calls_svlens (
        vid VARCHAR,
        svtype VARCHAR,
        svlen INTEGER,
        sample VARCHAR
    );
    CREATE TABLE joined_raw_calls_clusters (
        vid VARCHAR,
        member VARCHAR
    );
    EOF

    while read -r f; do
      duckdb joined_raw_calls.duckdb "COPY joined_raw_calls_svlens FROM '${f}' (FORMAT CSV, DELIMITER '\t', HEADER false);"
    done < '~{write_lines(tidy_vcfs)}'

    while read -r f; do
      duckdb joined_raw_calls.duckdb "COPY joined_raw_calls_clusters FROM '${f}' (FORMAT CSV, DELIMITER '\t', HEADER false);"
    done < '~{write_lines(clusters)}'
  >>>

  output {
    File joined_raw_calls_db = 'joined_raw_calls.duckdb'
  }
}

task CountSVsPerGenome {
  input {
    File joined_raw_calls_db
    Array[String] svtypes_to_filter

    String runtime_docker
  }

  Int disk_size_gib = ceil(size(joined_raw_calls_db, 'GiB') * 2.0) + 8

  runtime {
    memory: '4GiB'
    disks: 'local-disk ${disk_size_gib} HDD'
    cpus: 1
    preemptible: 3
    maxRetries: 0
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
    FROM '~{write_lines(svtypes_to_filter)}' (
        FORMAT CSV,
        DELIMITER ';',
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

    python3 '/opt/outlier-exclusion/scripts/count_svs.py' \
      sv_counts.duckdb '~{joined_raw_calls_db}'
  >>>

  output {
    File sv_counts_db = 'sv_counts.duckdb'
  }
}

task DetermineOutlierSamples {
  input {
    File sv_counts_db
    File wgd_scores
    Float min_wgd_score
    Float max_wgd_score
    Float iqr_multiplier

    String runtime_docker
  }

  Float input_size = size([sv_counts_db, wgd_scores], 'GiB')
  Int disk_size_gib = ceil(input_size * 1.5) + 8

  command <<<
    set -euo pipefail

    mv '~{sv_counts_db}' sv_counts_with_outliers.duckdb
    python3 '/opt/outlier-exclusion/scripts/determine_outlier_samples.py' \
      sv_counts_with_outliers.duckdb \
      '~{iqr_multiplier}' \
      '~{wgd_scores}' \
      '~{min_wgd_score}' \
      '~{max_wgd_score}'

    python3 '/opt/outlier-exclusion/scripts/dump_outlier_samples.py' \
      sv_counts_with_outliers.duckdb \
      dump

    printf 'sample_id\tcount\tsvtype\tmin_svlen\tmax_svlen\n' > outlier_samples.tsv
    find dump -type f -name '*.tsv' -exec cat '{}' \; >> outlier_samples.tsv
  >>>

  runtime {
    memory: '4GiB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gib} HDD'
    preemptible: 1
    maxRetries: 1
    docker: runtime_docker
  }

  output {
    File sv_counts_db_with_outliers = 'sv_counts_with_outliers.duckdb'
    File outlier_samples = 'outlier_samples.tsv'
  }
}

task DetermineOutlierVariants {
  input {
    Array[File] clustered_depth_vcfs
    Array[File] clustered_manta_vcfs
    Array[File] clustered_wham_vcfs
    Array[File] clustered_melt_vcfs
    Array[File] clustered_depth_vcf_indicies
    Array[File] clustered_manta_vcf_indicies
    Array[File] clustered_wham_vcf_indicies
    Array[File] clustered_melt_vcf_indicies
    File sv_counts_db
    File joined_raw_calls_db
    Float min_outlier_sample_prop

    String runtime_docker
  }

  Array[File] clusterbatch_vcfs = flatten([
    clustered_depth_vcfs, clustered_manta_vcfs, clustered_wham_vcfs,
    clustered_melt_vcfs
  ])
  Array[File] clusterbatch_vcf_indicies = flatten([
    clustered_depth_vcf_indicies, clustered_manta_vcf_indicies, clustered_wham_vcf_indicies,
    clustered_melt_vcf_indicies
  ])

  Int disk_size_gib = ceil(
    size(clusterbatch_vcfs, 'GiB') * 2.0
    + size(joined_raw_calls_db, 'GiB')
  ) + 16

  runtime {
    memory: '4GiB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gib} HDD'
    preemptible: 1
    maxRetries: 0
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
      | awk -F'/' '{bn=$NF; sub(/\.cluster_batch\.(depth|wham|manta|melt)\.vcf.gz$/, "", bn); print bn}' \
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
    test -r "${depth}" && query_vcf "${depth}" "${tsv}"
    test -r "${wham}" && query_vcf "${wham}" "${tsv}"
    test -r "${manta}" && query_vcf "${manta}" "${tsv}"
    test -r "${melt}" && query_vcf "${melt}" "${tsv}"
    duckdb "${db}" "COPY variants FROM '${tsv}' (FORMAT CSV, HEADER false, DELIMITER '\t');"
    EOF
    xargs -L 1 -P 0 -- bash reformat.bash < clusterbatch_ids.list

    python3 '/opt/outlier-exclusion/scripts/determine_outlier_variants.py' \
      '~{sv_counts_db}' \
      '~{joined_raw_calls_db}' \
      clusterbatch_dbs \
      '~{min_outlier_sample_prop}' \
      | LC_ALL=C sort -u > joined_raw_calls_outlier_variants.list
  >>>

  output {
    File outlier_variants = 'joined_raw_calls_outlier_variants.list'
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

  Int disk_size_gib = ceil(size(filtered_vcf, 'GiB') * 2.0) + 8

  runtime {
    memory: '2GiB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gib} HDD'
    preemptible: 1
    maxRetries: 0
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
      --output '~{cohort_prefix}_outliers_annotated_calls.vcf.gz' \
      --output-type z \
      --write-index=tbi \
      '~{filtered_vcf}'
  >>>

  output {
    File outlier_annotated_vcf = '~{cohort_prefix}_outliers_annotated_calls.vcf.gz'
    File outlier_annotated_vcf_index = '~{cohort_prefix}_outliers_annotated_calls.vcf.gz.tbi'
  }
}
