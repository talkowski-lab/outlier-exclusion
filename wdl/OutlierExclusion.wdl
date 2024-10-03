version 1.0

workflow OutlierExclusion {
  input {

    # MakeJoinedRawCallsDB ------------------------------------------------------
    File joined_raw_calls_vcf
    File joined_raw_calls_vcf_index
    String docker

    # CountSVsPerGenome -------------------------------------------------------
    File count_svs_script
    Array[String] svtypes_to_filter = ['DEL;-Inf;Inf', 'DUP;-Inf;Inf']

    # DetermineOutlierSamples -------------------------------------------------
    File wgd_scores
    File determine_outlier_samples_script
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
    File determine_outlier_variants_script
    File concordance_vcf
    File concordance_vcf_index
    Float min_outlier_sample_prop = 1.0

    # FlagOutlierVariants -----------------------------------------------------
    String cohort_prefix

    # ApplyManualFilter -------------------------------------------------------
    String filter_name = 'manual_filter'
    String? bcftools_filter
  }

  call MakeJoinedRawCallsDB {
    input:
      joined_raw_calls_vcf = joined_raw_calls_vcf,
      joined_raw_calls_vcf_index = joined_raw_calls_vcf_index,
      runtime_docker = docker
  }

  call CountSVsPerGenome {
    input:
      joined_raw_calls_db = MakeJoinedRawCallsDB.joined_raw_calls_db,
      svtypes_to_filter = svtypes_to_filter,
      count_svs_script = count_svs_script,
      runtime_docker = docker
  }

  call DetermineOutlierSamples {
    input:
      sv_counts_db = CountSVsPerGenome.sv_counts_db,
      wgd_scores = wgd_scores,
      min_wgd_score = min_wgd_score,
      max_wgd_score = max_wgd_score,
      iqr_multiplier = iqr_multiplier,
      determine_outlier_samples_script = determine_outlier_samples_script,
      runtime_docker = docker
  }

  call DetermineOutlierVariants {
    input:
      clustered_depth_vcfs = clustered_depth_vcfs,
      clustered_manta_vcfs = clustered_manta_vcfs,
      clustered_wham_vcfs = clustered_wham_vcfs,
      clustered_melt_vcfs = clustered_melt_vcfs,
      clustered_depth_vcf_indicies = clustered_depth_vcf_indicies,
      clustered_manta_vcf_indicies = clustered_manta_vcf_indicies,
      clustered_wham_vcf_indicies = clustered_wham_vcf_indicies,
      clustered_melt_vcf_indicies = clustered_melt_vcf_indicies,
      sv_counts_db = DetermineOutlierSamples.sv_counts_db_with_outliers,
      joined_raw_calls_db = MakeJoinedRawCallsDB.joined_raw_calls_db,
      runtime_docker = docker
  }

  call FlagOutlierVariants {
    input:
      cohort_prefix = cohort_prefix,
      concordance_vcf = concordance_vcf,
      concordance_vcf_index = concordance_vcf_index,
      outlier_variants = DetermineOutlierVariants.outlier_variants,
      runtime_docker = docker
  }

  if (defined(bcftools_filter)) {
    call ApplyManualFilter {
      input:
        cohort_prefix = cohort_prefix,
        vcf = FlagOutlierVariants.outlier_annotated_vcf,
        vcf_index = FlagOutlierVariants.outlier_annotated_vcf_index,
        filter_name = filter_name,
        bcftools_filter = select_first([bcftools_filter, '']),
        runtime_docker = docker
    }
    File manual_filtered_vcf = ApplyManualFilter.hard_filtered_vcf
    File manual_filtered_vcf_index = ApplyManualFilter.hard_filtered_vcf_index
  }

  output {
    File manual_filtered_and_flagged_vcf = select_first([manual_filtered_vcf, FlagOutlierVariants.outlier_annotated_vcf])
    File manual_filtered_and_flagged_vcf_index = select_first([manual_filtered_vcf_index, FlagOutlierVariants.outlier_annotated_vcf_index])
  }
}

task MakeJoinedRawCallsDB {
  input {
    File joined_raw_calls_vcf
    File joined_raw_calls_vcf_index

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(joined_raw_calls_vcf, 'GB') * 10.0 + 8.0)

  runtime {
    memory: '2 GB'
    disks: 'local-disk ${disk_size_gb} HDD'
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

    bcftools query \
      --include 'GT ~ "1" & INFO/SVLEN > 0' \
      --format '[%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%SAMPLE\n]' \
      '~{joined_raw_calls_vcf}' \
      | awk -F'\t' '$2 ~ /^INS/{$2 = "INS"; print}' OFS='\t' > joined_raw_calls_svlens.tsv

    bcftools query \
      --format '%ID\t%INFO/MEMBERS\n' \
      filtered.bcf \
      | awk -F'\t' '{split($2, a, /,/); for (i in a) print $1"\t"a[i]}' \
      > joined_raw_calls_clusters.tsv

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
    COPY joined_raw_calls_svlens
    FROM 'joined_raw_calls_svlens.tsv' (
        FORMAT CSV,
        DELIMITER '\t',
        HEADER false
    );
    COPY joined_raw_calls_clusters
    FROM 'joined_raw_calls_clusters.tsv' (
        FORMAT CSV,
        DELIMITER '\t',
        HEADER false
    );
    EOF
  >>>

  output {
    File joined_raw_calls_db = 'joined_raw_calls.duckdb'
  }
}

task CountSVsPerGenome {
  input {
    File joined_raw_calls_db
    Array[String] svtypes_to_filter
    File count_svs_script

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(joined_raw_calls_db, 'GB') * 2.0 + 8.0)

  runtime {
    memory: '4 GB'
    disks: 'local-disk ${disk_size_gb} HDD'
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

    python3 '~{count_svs_script}' sv_counts.duckdb '~{joined_raw_calls_db}'
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

    File determine_outlier_samples_script

    String runtime_docker
  }

  Float input_size = size([sv_counts_db, wgd_scores], 'GB')
  Int disk_size_gb = ceil(input_size * 1.5) + 8

  command <<<
    set -euo pipefail

    python3 '~{determine_outlier_samples_script}' \
      '~{sv_counts_db}' \
      '~{iqr_multiplier}' \
      '~{wgd_scores}' \
      '~{min_wgd_score}' \
      '~{max_wgd_score}'
  >>>

  runtime {
    memory: '4 GB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 1
    maxRetries: 1
    docker: runtime_docker
  }

  output {
    File sv_counts_db_with_outliers = sv_counts_db
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

    File determine_outlier_variants_script

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

  Int disk_size_gb = ceil(
    size(clusterbatch_vcfs, 'GB') * 2.0
    + size(joined_raw_calls_db, 'GB')
    + 16.0
  )

  runtime {
    memory: '4 GB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 1
    maxRetries: 0
    docker: runtime_docker
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    mkdir clusterbatch_vcfs
    cat '~{clusterbatch_vcfs}' '~{clusterbatch_vcf_indicies}' | while read -r vcf; do
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
      bcftools query --include 'GT ~ "1"' \
        --format '[%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%SAMPLE\n]' "$1" >> "$2"
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

    python3 '~{determine_outlier_variants_script}' \
      '~{sv_counts_db}' \
      '~{joined_raw_calls_db}' \
      clusterbatch_dbs \
      '~{min_outlier_sample_prop}' \
      | LC_ALL=C sort -u > joined_raw_calls_outlier_variants.list
  >>>

  output {
    File outlier_variants = 'joined_raw_calls_outliers_variants.list'
  }
}

task FlagOutlierVariants {
  input {
    String cohort_prefix
    File concordance_vcf
    File concordance_vcf_index
    File outlier_variants

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(concordance_vcf) * 2.0 + 8.0)

  runtime {
    memory: '2 GB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 1
    maxRetries: 0
    docker: runtime_docker
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    bcftools query --include 'INFO/TRUTH_VID != ""' \
      --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/TRUTH_VID\n' \
      '~{concordance_vcf}' \
      | LC_ALL=C sort -k6,6 > concordance_variants.tsv
    LC_ALL=C join -1 6 -2 1 -o 1.1,1.2,1.3,1.4,1.5 -t $'\t' \
      concordance_variants.tsv \
      '~{outlier_variants}' > concordance_calls_outliers.tsv

    awk -F'\t' '{print $0"\toutlier"}' 'concordance_calls_outliers.tsv' \
      | sort -k1,1 -k2,2n > annotations.tsv
    bgzip annotations.tsv
    tabix --begin 2 --end 2 --sequence 1 annotations.tsv.gz

    printf '##FILTER=<ID=outlier,Description="Variant enriched by outlier samples">' \
      > header.txt

    bcftools annotate \
      --annotations annotations.tsv.gz \
      --columns 'CHROM,POS,REF,ALT,~ID,.FILTER' \
      --header-lines header.txt \
      --output '~{cohort_prefix}_outliers_annotated_concordance_calls.vcf.gz' \
      --output-type z \
      --write-index=tbi \
      '~{concordance_vcf}'
  >>>

  output {
    File outlier_annotated_vcf = '~{cohort_prefix}_outliers_annotated_concordance_calls.vcf.gz'
    File outlier_annotated_vcf_index = '~{cohort_prefix}_outliers_annotated_concordance_calls.vcf.gz.tbi'
  }
}

task ApplyManualFilter {
  input {
    String cohort_prefix
    File vcf
    File vcf_index
    String filter_name
    String bcftools_filter

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(vcf, 'GB') * 2.0 + 8.0)

  runtime {
    memory: '1 GB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 1
    maxRetries: 1
    docker: runtime_docker
  }

  command <<<
    set -euo pipefail

    bcftools view -e '~{bcftools_filter}' ~{vcf} -Oz \
      -o '~{cohort_prefix}.~{filter_name}.vcf.gz'
    tabix '~{cohort_prefix}.~{filter_name}.vcf.gz'
  >>>

  output {
    File hard_filtered_vcf = '~{cohort_prefix}.~{filter_name}.vcf.gz'
    File hard_filtered_vcf_index = '~{cohort_prefix}.~{filter_name}.vcf.gz.tbi'
  }
}
