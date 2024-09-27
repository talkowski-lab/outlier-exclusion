version 1.0

workflow OutlierExclusion {
  input {

    # MakeJoinedRawCallsDB ------------------------------------------------------
    File joined_raw_calls_vcf
    File joined_raw_calls_vcf_index
    File duckdb_zip
    String bcftools_docker

    # CountSVsPerGenome -------------------------------------------------------
    String countsvs_svtype
    Int countsvs_max_svlen
    Int countsvs_min_svlen
    String linux_docker

    # DetermineOutlierSamples -------------------------------------------------
    String cohort_prefix
    File wgd_scores
    File determine_outlier_samples_script
    Float min_wgd_score
    Float max_wgd_score
    Int iqr_multiplier
    String svtk_docker

    # DetermineOutlierVariants ------------------------------------------------
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
      make_tables_sql = make_tables_sql,
      duckdb_zip = duckdb_zip,
      runtime_docker = bcftools_docker
  }

  call CountSVsPerGenome {
    input:
      svtype = countsvs_svtype,
      min_svlen = countsvs_min_svlen,
      max_svlen = countsvs_max_svlen
      MakeJoinedRawCallsDB.joined_raw_calls_db,
      templater_awk = templater_awk,
      count_svs_sql_tmpl = count_svs_sql_tmpl,
      duckdb_zip = duckdb_zip,
      runtime_docker = linux_docker
  }

  call DetermineOutlierSamples {
    input:
      cohort_prefix = cohort_prefix,
      sv_counts_per_genome_all = CountSVsPerGenome.sv_counts_per_genome_all,
      sv_counts_per_genome_filtered = CountSVsPerGenome.sv_counts_per_genome_filtered,
      wgd_scores = wgd_scores,
      min_wgd_score = min_wgd_score,
      max_wgd_score = max_wgd_score,
      iqr_multiplier = iqr_multiplier,
      determine_outlier_samples_script = determine_outlier_samples_script,
      runtime_docker = svtk_docker
  }

  call DetermineOutlierVariants {
    input:
      cohort_prefix = cohort_prefix,
      outlier_samples = DetermineOutlierSamples.outlier_samples_list,
      clustered_depth_vcfs = clustered_depth_vcfs,
      clustered_manta_vcfs = clustered_manta_vcfs,
      clustered_wham_vcfs = clustered_wham_vcfs,
      clustered_melt_vcfs = clustered_melt_vcfs,
      clustered_depth_vcf_indicies = clustered_depth_vcf_indicies,
      clustered_manta_vcf_indicies = clustered_manta_vcf_indicies,
      clustered_wham_vcf_indicies = clustered_wham_vcf_indicies,
      clustered_melt_vcf_indicies = clustered_melt_vcf_indicies
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
      reformat_vcf_header_script = reformat_vcf_header_script
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
      vcf = Remove_outlier_samples.outlier_samples_removed_vcf,
      vcf_index = Remove_outlier_samples.outlier_samples_removed_vcf_index,
      filter_name = filter_name,
      bcftools_filter = bcftools_filter,
      runtime_docker = bcftools_docker,
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

  Int disk_size_gb = ceil(size(join_raw_calls_vcf, 'GB') * 10.0)

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

    bcftools filter \
      --include 'INFO/SVTYPE ~ "^DEL" || INFO/SVTYPE ~ "^DUP" || INFO/SVTYPE ~ "^INV" || INFO/SVTYPE ~ "^INS"' \
      --output-type u
      --output filtered.bcf \
      '~{joined_raw_calls_vcf}'

    bcftools query \
      --include 'GT ~ "1"' \
      --format '[%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%SAMPLE\n]' \
      filtered.bcf > joined_raw_calls_svlens.tsv

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
    String svtype
    Int max_svlen
    Int min_svlen
    File joined_raw_calls_db

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

    ./duckdb '~{joined_raw_calls_db}' << 'EOF'
    COPY
    (SELECT sample, svtype, COUNT(*) AS count
    FROM joined_raw_calls_svlens
    GROUP BY sample, svtype)
    TO 'joined_raw_calls_svcounts_ALL.tsv'
    (DELIMITER '\t', HEADER true);

    PREPARE query_svtype AS
    SELECT sample, COUNT(*) AS count
    FROM joined_raw_calls_svlens
    WHERE svtype = ? AND svlen >= ? AND svlen <= ?
    GROUP BY sample;

    .mode tabs
    .once 'joined_raw_calls_svcount_~{svtype}.tsv'
    EXECUTE query_svtype('~{svtype}', ~{min_svlen}, ~{max_svlen});
    EOF
  >>>

  output {
    File sv_counts_per_genome_all = 'joined_raw_calls_svcounts_ALL.tsv'
    File sv_counts_per_genome_filtered = 'joined_raw_calls_svcount_~{svtype}.tsv'
  }
}

task DetermineOutlierSamples {
  input {
    String cohort_prefix
    File sv_counts_per_genome_all
    File sv_counts_per_genome_filtered
    File wgd_scores
    Float min_wgd_score
    Float max_wgd_score
    Int iqr_multiplier

    File determine_outlier_samples_script

    String runtime_docker
  }

  Float input_size = size([sv_counts_per_genome_all, sv_counts_per_genome_filtered, wgd_scores], 'GB')
  Int mem_gb = ceil(input_size * 1.2)
  Int disk_size_gb = ceil(input_size * 1.3) + 8

  command <<<
    set -euo pipefail

    python3 '~{determine_outlier_samples_script}' \
      -s '~{sv_counts_per_genome}' \
      -r '~{sv_counts_per_genome_filtered}'
      -i ~{iqr_multiplier} \
      -w '~{wgd_scores}' \
      -l '~{wgd_lower_cutoff}' \
      -hi '~{wgd_higher_cutoff}' \
      -o '~{cohort_prefix}_outlier_sample_list.txt'
  >>>

  runtime {
    memory: '~{mem_gb} GB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 1
    maxRetries: 1
    docker: runtime_docker
  }

  output {
    File outlier_samples_list = '~{cohort_prefix}_outlier_sample_list.txt'
  }
}

task DetermineOutlierVariants {
  input {
    String cohort_prefix
    File outlier_samples
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
    File joined_raw_calls_db

    File duckdb_zip
    String runtime_docker
  }

  Array[File] clusterbatch_vcfs = flatten([
    clustered_depth_vcfs, clustered_manta_vcfs, clustered_wham_vcfs,
    clustered_melt_vcfs
  ])

  Int disk_size_gb = ceil(
    size(clusterbatch_vcfs, 'GB')
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

    ./duckdb outlier_samples.duckdb << 'EOF'
    CREATE TABLE outlier_samples (
        sample VARCHAR
    );

    COPY outlier_samples
    FROM '~{outlier_samples}' (
        FORMAT CSV,
        DELIMITER '\t',
        HEADER false
    );

    CREATE TABLE variants (
        vid VARCHAR,
        sample VARCHAR
    );
    EOF

    while read -r vcf; do
      bcftools query --include 'GT ~ "1"' \
        --format '[%ID\t%SAMPLE\n]' \
        "${vcf}"
    done < '~{write_lines(clusterbatch_vcfs)}' \
      | ./duckdb outliers_samples.duckdb "COPY variants FROM '/dev/stdin' (FORMAT CSV, DELIMITER '\t', HEADER false);"

    ./duckdb outliers_samples.duckdb > clusterbatch_outliers.list << 'EOF'
    .mode tabs
    SELECT vid
    FROM (
        SELECT l.vid AS vid, COUNT(*) FILTER(r.sample IS NULL) AS nils
        FROM variants l
        LEFT JOIN outlier_samples r
        USING (sample)
        GROUP BY vid
    )
    WHERE nils = 0;
    EOF

    ./duckdb '~{joined_raw_calls_db}' > 'joined_raw_calls_outliers.list' << 'EOF'
    CREATE TEMP TABLE t1 (
        vid VARCHAR
    );

    COPY t1
    FROM 'clusterbatch_outliers.list' (
        FORMAT CSV,
        DELIMITER '\t',
        HEADER false
    );

    .mode tabs
    .header off
    SELECT jrc.vid
    FROM joined_raw_calls_clusters jrc
    JOIN t1 ON (jrc.member = t1.vid)
    ORDER BY jrc.vid ASC;
    EOF

    bcftools query --include 'INFO/TRUTH_VID != ""' \
      --format '%ID\t%INFO/TRUTH_VID\n' \
      '~{concordance_vcf}' \
      | LC_ALL=C sort -k2,2 > concordance_vids.tsv
    LC_ALL=C join -1 2 -2 1 -o 1.1 \
      concordance_vids.tsv \
      joined_raw_calls_outliers.list > '~{cohort_prefix}_concordance_calls_outlier_vids.list'
  >>>

  output {
    File outlier_variants = '${cohort_prefix}_concordance_calls_outlier_vids.list'
  }
}

task ReformatConcordanceVCF {
  input {
    String cohort_prefix
    File concordance_vcf
    File concordance_vcf_index
    File reformat_vcf_header_script

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(concordance_vcf) * 2.0 + 8.0)
  command <<<
    set -euo pipefail

    bcftools view -h '~{concordance_vcf}' \
      | grep '##FILTER' > concordance_header.txt

    echo '##FILTER=<ID=outlier_discovered,Description="Variant with only outlier samples associated">' >> concordance_header.txt

    python3 '~{reformat_vcf_header_script}' \
      -i '~{concordance_vcf}' \
      -hf concordance_header.txt \
      -o '~{cohort_prefix}_reformatted_concordance.vcf.gz'

    tabix -f '~{cohort_prefix}_reformatted_concordance.vcf.gz'
  >>>

  runtime {
    memory: '2 GB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 1
    maxRetries: 0
    docker: runtime_docker
  }

  output {
    File reformatted_concordance_vcf = '~{cohort_prefix}_reformatted_concordance.vcf.gz'
    File reformatted_concordance_vcf_index = '~{cohort_prefix}_reformatted_concordance.vcf.gz.tbi'
  }
}

task FlagOutlierVariantsStep0 {
  input {
    File reformatted_concordance_vcf
    File reformatted_concordance_vcf_index
    File concordance_outlier_variants
    String cohort_prefix
    File outlier_samples
    Int min_svlen
    Int max_svlen
    String svtype

    File update_outlier_discovery_flag_script
    File outlier_size_range_variants_script
    File final_update_outlier_discovery_flag_script

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(reformatted_concordance_vcf, 'GB') * 5.0 + 8.0)

  runtime {
    memory: '4 GB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb} HDD'
    preemptible: 1
    docker: docker
  }

  command <<<
    set -euo pipefail

    python3 '~{update_outlier_discovery_flag_script}' \
      -i '~{reformatted_concordance_vcf}' \
      -v '~{concordance_outlier_variants}' \
      -o '~{cohort_prefix}_outlier_flagged.vcf.gz'
    tabix '~{cohort_prefix}_outlier_flagged.vcf.gz'

    svtk vcf2bed --include-filters -i ALL \
      '~{cohort_prefix}_outlier_flagged.vcf.gz' \
      '~{cohort_prefix}_outlier_flagged.bed'

    python3 '~{outlier_size_range_variants_script}' \
      -o '~{outlier_samples_file}' \
      -i '~{cohort_prefix}_outlier_flagged.bed' \
      -out '~{cohort_prefix}_outlier_size_range_filtered_variants' \
      -t '~{svtype}' \
      -l ~{svtype_size_range_lower_cutoff}
      -hi ~{svtype_size_range_higher_cutoff}

    python3 '~{final_update_outlier_discovery_flag_script}' \
      -i '~{cohort_prefix}_outlier_flagged.vcf.gz' \
      -o '~{cohort_prefix}_outliers_flagged_step0.vcf.gz' \
      -f '~{cohort_prefix}_outlier_size_range_filtered_variants'
    tabix '~{cohort_prefix}_outlier_flagged_step0.vcf.gz'
  >>>

  output {
    File outlier_flagged_vcf = '~{cohort_prefix}_outlier_flagged_step0.vcf.gz'
    File outlier_flagged_vcf_index= '~{cohort_prefix}_outlier_flagged_step0.vcf.gz.tbi'
  }
}

task FlagOutlierVariantsStep1 {
  input {
    File step0_outlier_flagged_vcf
    File step0_outlier_flagged_vcf_index
    String cohort_prefix
    File outlier_samples
    String svtype
    Int min_svlen
    Int max_svlen
    Float fraction_of_outlier_samples

    File flag_outlier_variants_based_on_size_range_counts_script
    File proportion_of_outlier_samples_associated_with_variant_script

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(step0_outlier_flagged_vcf, 'GB') * 5.0 + 8.0)

  runtime {
    memory: '4 GB'
    cpu: 1
    bootDiskSizeGb: 16
    disks: 'local-disk ${disk_size_gb}  HDD'
    preemptible: 1
    maxRetries: 0
    docker: runtime_docker
  }

  command <<<
    set -euo pipefail

    svtk vcf2bed --include-filters \
      -i ALL \
      '~{step0_outlier_flagged_vcf}' \
      step0_outlier_flagged.bed
 
    awk -F'\t' '$5 == ~{svtype} && ($3 - $2) >= ~{min_svlen} && ($3 - $2) <= ~{max_svlen}' \
      step0_outlier_flagged.bed \
      > step0_outlier_flagged_filtered.bed 

~{cohort_prefix}_filtered_outlier_discovered_flag_added_to_size_range_outliers_FINAL.bed

    python3 ~{proportion_of_outlier_samples_associated_with_variant_script} \
      -o '~{outlier_samples}' \
      -i step0_outlier_flagged_filtered.bed \
      -out outlier_variants_in_size_range.tsv \
      -f ~{fraction_of_outlier_samples}

    python3 ~{flag_outlier_variants_based_on_size_range_counts_script} \
      -v outlier_variants_in_size_range.tsv \
      '~{step0_outlier_flagged_vcf}' \
      -o '~{cohort_prefix}_outlier_flagged_step1.vcf.gz'
    tabix '~{cohort_prefix}_outlier_flagged_step1.vcf.gz' 
  >>>

  output {
    File outlier_flagged_vcf = '~{cohort_prefix}_outlier_flagged_step1.vcf.gz'
    File outlier_flagged_vcf_index = '~{cohort_prefix}_outlier_flagged_step1.vcf.gz.tbi'
  }
}

task RemoveOutlierSamples {
  input {
    File final_flagged_vcf
    File final_flagged_vcf_index
    File outlier_samples
    String cohort_prefix

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(final_flagged_vcf, 'GB') + 8.0)

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

    bcftools view \
      -c 1
      --samples-file ^~{outlier_samples} \
      -o '~{cohort_prefix}_outliers_removed.vcf.gz' \
      -O z \
      '~{final_flagged_vcf}'
    tabix '~{cohort_prefix}_outliers_removed.vcf.gz'
  >>>

  output {
    File outliers_removed_vcf = '~{cohort_prefix}_outliers_removed.vcf.gz'
    File outliers_removed_vcf_index = '~{cohort_prefix}_outliers_removed.vcf.gz.tbi'
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
