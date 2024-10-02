version 1.0

workflow OutlierExclusion {
  input {

    # MakeJoinedRawCallsDB ------------------------------------------------------
    File joined_raw_calls_vcf
    File joined_raw_calls_vcf_index
    File duckdb_zip = 'https://github.com/duckdb/duckdb/releases/download/v1.1.1/duckdb_cli-linux-amd64.zip'
    String docker

    # CountSVsPerGenome -------------------------------------------------------
    File count_svs_script
    Array[String] svtypes_to_filter = ['DEL;-Inf;Inf', 'DUP;-Inf;Inf']

    # DetermineOutlierSamples -------------------------------------------------
    String cohort_prefix
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

    # ReformatConcordanceVCF --------------------------------------------------
    File reformat_vcf_header_script

    # FlagOutlierVariantsStep0 ------------------------------------------------
    File update_outlier_discovery_flag_script
    File outlier_size_range_variants_script
    File final_update_outlier_discovery_flag_script

    # FlagOutlierVariantsStep1 ------------------------------------------------
    Float fraction_of_outlier_samples
    File flag_outlier_variants_based_on_size_range_counts_script
    File proportion_of_outlier_samples_associated_with_variant_script

    # ApplyManualFilter -------------------------------------------------------
    String filter_name
    String bcftools_filter
  }

  call MakeJoinedRawCallsDB {
    input:
      joined_raw_calls_vcf = joined_raw_calls_vcf,
      joined_raw_calls_vcf_index = joined_raw_calls_vcf_index,
      duckdb_zip = duckdb_zip,
      runtime_docker = bcftools_docker
  }

  call CountSVsPerGenome {
    input:
      joined_raw_calls_db = MakeJoinedRawCallsDB.joined_raw_calls_db,
      svtypes_to_filter = svtypes_to_filter,
      count_svs_script = count_svs_script,
      duckdb_zip = duckdb_zip,
      runtime_docker = linux_docker
  }

  call DetermineOutlierSamples {
    input:
      sv_counts_db = CountSVsPerGenome.sv_counts_db,
      wgd_scores = wgd_scores,
      min_wgd_score = min_wgd_score,
      max_wgd_score = max_wgd_score,
      iqr_multiplier = iqr_multiplier,
      determine_outlier_samples_script = determine_outlier_samples_script,
      runtime_docker = svtk_docker
  }

  call DetermineOutlierVariants {
    input:
      outlier_samples = DetermineOutlierSamples.sv_counts_db,
      clustered_depth_vcfs = clustered_depth_vcfs,
      clustered_manta_vcfs = clustered_manta_vcfs,
      clustered_wham_vcfs = clustered_wham_vcfs,
      clustered_melt_vcfs = clustered_melt_vcfs,
      clustered_depth_vcf_indicies = clustered_depth_vcf_indicies,
      clustered_manta_vcf_indicies = clustered_manta_vcf_indicies,
      clustered_wham_vcf_indicies = clustered_wham_vcf_indicies,
      clustered_melt_vcf_indicies = clustered_melt_vcf_indicies,
      concordance_vcf = concordance_vcf,
      concordance_vcf_index = concordance_vcf_index,
      joined_raw_calls_db = MakeJoinedRawCallsDB.joined_raw_calls_db,
      duckdb_zip = duckdb_zip,
      runtime_docker = bcftools_docker
  }

  call ReformatConcordanceVCF {
    input:
      cohort_prefix = cohort_prefix,
      concordance_vcf = concordance_vcf,
      concordance_vcf_index = concordance_vcf_index,
      reformat_vcf_header_script = reformat_vcf_header_script,
      runtime_docker = svtk_docker
  }

  call FlagOutlierVariantsStep0 {
    input:
      reformatted_concordance_vcf = ReformatConcordanceVCF.reformatted_concordance_vcf,
      reformatted_concordance_vcf_index = ReformatConcordanceVCF.reformatted_concordance_vcf_index,
      concordance_outlier_variants = DetermineOutlierVariants.outlier_variants,
      cohort_prefix = cohort_prefix,
      outlier_samples = DetermineOutlierSamples.outlier_samples,
      min_svlen = countsvs_min_svlen,
      max_svlen = countsvs_max_svlen,
      svtype = countsvs_svtype,
      update_outlier_discovery_flag_script = update_outlier_discovery_flag_script,
      outlier_size_range_variants_script = outlier_size_range_variants_script,
      final_update_outlier_discovery_flag_script = final_update_outlier_discovery_flag_script,
      runtime_docker = svtk_docker
  }

  call FlagOutlierVariantsStep1 {
    input:
      step0_outlier_flagged_vcf = FlagOutlierVariantsStep0.outlier_flagged_vcf,
      step0_outlier_flagged_vcf_index = FlagOutlierVariantsStep0.outlier_flagged_vcf_index,
      cohort_prefix = cohort_prefix,
      outlier_samples = DetermineOutlierSamples.outlier_samples,
      svtype = countsvs_svtype,
      min_svlen = countsvs_min_svlen,
      max_svlen = countsvs_max_svlen,
      fraction_of_outlier_samples = fraction_of_outlier_samples,
      flag_outlier_variants_based_on_size_range_counts_script = flag_outlier_variants_based_on_size_range_counts_script,
      proportion_of_outlier_samples_associated_with_variant_script = proportion_of_outlier_samples_associated_with_variant_script,
      runtime_docker = svtk_docker
  }

  call RemoveOutlierSamples {
    input:
      final_flagged_vcf = FlagOutlierVariantsStep1.outlier_flagged_vcf,
      final_flagged_vcf_index = FlagOutlierVariantsStep1.outlier_flagged_vcf_index,
      outlier_samples = DetermineOutlierSamples.outlier_samples,
      cohort_prefix = cohort_prefix,
      runtime_docker = bcftools_docker
  }

  call ApplyManualFilter {
    input:
      cohort_prefix = cohort_prefix,
      vcf = RemoveOutlierSamples.outliers_removed_vcf,
      vcf_index = RemoveOutlierSamples.outliers_removed_vcf_index,
      filter_name = filter_name,
      bcftools_filter = bcftools_filter,
      runtime_docker = bcftools_docker
  }

  output {
    File manual_filtered_and_flagged_vcf = ApplyManualFilter.hard_filtered_vcf
    File manual_filtered_and_flagged_vcf_index = ApplyManualFilter.hard_filtered_vcf_index
  }
}

task MakeJoinedRawCallsDB {
  input {
    File joined_raw_calls_vcf
    File joined_raw_calls_vcf_index

    File duckdb_zip
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

    unzip '~{duckdb_zip}'
    chmod u+x ./duckdb
     
    ./duckdb joined_raw_calls.duckdb << 'EOF'
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

    File duckdb_zip
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

    unzip '~{duckdb_zip}'
    chmod u+x ./duckdb

    ./duckdb sv_counts.duckdb << 'EOF'
    CREATE TABLE sv_filters (
        svtype VARCHAR,
        min_svlen DOUBLE,
        max_svlen DOUBLE
    );
    COPY sv_filters
    FROM '~{write_line(svtypes_to_filter)}' (
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
    File sv_counts_db = 'sv_counts_db'
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
    File concordance_vcf
    File concordance_vcf_index
    File sv_counts_db
    File joined_raw_calls_db
    Float min_outlier_sample_prop

    File determine_outlier_variants_script

    File duckdb_zip
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
    + size(concordance_vcf, 'GB')
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

    unzip '~{duckdb_zip}'
    chmod u+x ./duckdb

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
    ./duckdb "${db}" 'CREATE TABLE variants (vid VARCHAR, svtype VARCHAR, svlen INTEGER, sample VARCHAR);'
    depth="clusterbatch_vcfs/${1}.cluster_batch.depth.vcf.gz"
    wham="clusterbatch_vcfs/${1}.cluster_batch.wham.vcf.gz"
    manta="clusterbatch_vcfs/${1}.cluster_batch.manta.vcf.gz"
    melt="clusterbatch_vcfs/${1}.cluster_batch.melt.vcf.gz"
    test -r "${depth}" && query_vcf "${depth}" "${tsv}"
    test -r "${wham}" && query_vcf "${wham}" "${tsv}"
    test -r "${manta}" && query_vcf "${manta}" "${tsv}"
    test -r "${melt}" && query_vcf "${melt}" "${tsv}"
    ./duckdb "${db}" "COPY variants FROM '${tsv}' (FORMAT CSV, HEADER false, DELIMITER '\t');"
    EOF
    xargs -L 1 -P 0 -- bash reformat.bash < clusterbatch_ids.list

    python3 '~{determine_outlier_variants_script}' \
      '~{sv_counts_db}' \
      '~{joined_raw_calls_db}' \
      clusterbatch_dbs \
      '~{min_outlier_sample_prop}' \
      LC_ALL=C sort -u > joined_raw_calls_outlier_variants.list

    bcftools query --include 'INFO/TRUTH_VID != ""' \
      --format '%CHROM\t%POS%ID\t%INFO/TRUTH_VID\n' \
      '~{concordance_vcf}' \
      | LC_ALL=C sort -k4,4 > concordance_vids.tsv
    LC_ALL=C join -1 4 -2 1 -o 1.1,1.2,1.3 -t $'\t' \
      concordance_vids.tsv \
      joined_raw_calls_outliers_variants.list > concordance_calls_outlier_vids.tsv
  >>>

  output {
    File outlier_variants = 'concordance_calls_outlier_vids.tsv'
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

    awk -F'\t' '{print $0"\toutlier"}' '~{outlier_variants}' > annotations.tsv

    bcftools annotate \
      --annotations annotations.tsv \
      --columns 'CHROM,POS,~ID,.FILTER' \
      --header-lines '##FILTER=<ID=outlier,Description="Variant enriched by outlier samples">' \
      --output '~{cohort_prefix}_outliers_annotated_concordance.vcf.gz' \
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
