version 1.0

workflow OutlierExclusion {
  input {

    # MakeJoinedRawCallsDB ------------------------------------------------------
    File joined_raw_calls_vcf
    File joined_raw_calls_vcf_index
    File make_tables_sql
    File duckdb_zip
    String bcftools_docker

    # CountSVsPerGenome -------------------------------------------------------
    String countsvs_svtype
    Int countsvs_max_svlen
    Int countsvs_min_svlen
    File templater_awk
    File count_svs_sql_tmpl
    String linux_docker

    # DetermineOutlierSamples -------------------------------------------------
    String cohort_prefix
    File wgd_scores
    File determine_outlier_samples_script
    Float min_wgd_score
    Float max_wgd_score
    Int iqr_multiplier
    String pandas_docker

    String filter_name
    String bcftools_filter

    File concordance_vcf
    File concordance_vcf_index

    String svtype
    Int min_svsize
    Int max_svsize
    # Int joinrawcalls_sv_counts_size_range_lower_cutoff
    # Int joinrawcalls_sv_counts_size_range_higher_cutoff
    Int iqr_multiplier

    File wgd_score_file
    Float min_wgd_score
    Float max_wgd_score
    # Float wgd_lower_cutoff
    # Float wgd_higher_cutoff
    Float fraction_of_outlier_samples

    Array[File] clustered_depth_vcfs
    Array[File] clustered_manta_vcfs
    Array[File] clustered_wham_vcfs
    Array[File] clustered_melt_vcfs
    #Array[File] clusterbatch_vcf_indexes

    # Scripts
    File determine_outlier_samples_script
    File identify_outlier_variants_script
    File identify_join_raw_calls_variants_script
    File reformat_vcf_header_script
    File update_outlier_discovery_flag_script
    File outlier_size_range_variants_script
    File final_update_outlier_discovery_flag_script
    File flag_outlier_variants_based_on_size_range_counts_script
    File proportion_of_outlier_samples_associated_with_variant_script

    String docker
    Int disk_size_gb
    Int machine_mem_mb
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
      runtime_docker = pandas_docker
  }

     call Concat_clusterbatch_vcfs {
         input:
            cohort_prefix = cohort_prefix,
            clustered_depth_vcfs = clustered_depth_vcfs,
            clustered_wham_vcfs = clustered_wham_vcfs,
            clustered_manta_vcfs = clustered_manta_vcfs,
            clustered_melt_vcfs = clustered_melt_vcfs,
            #clusterbatch_vcf_indexes = clusterbatch_vcf_indexes,
            disk_size_gb = disk_size_gb,
            docker = docker,
            machine_mem_mb = machine_mem_mb
     }

     call Get_outlier_variants {
         input:
            outlier_samples_file = Determine_outlier_samples.outlier_samples_list,
            merged_clusteredvcfs_bedfile = Concat_clusterbatch_vcfs.merged_clusteredvcfs_bed,
            join_raw_calls_bed = JoinRawCalls_sv_counts.join_raw_calls_bed,
            concordance_vcf = concordance_vcf,
            identify_outlier_variants_script = identify_outlier_variants_script,
            identify_join_raw_calls_variants_script = identify_join_raw_calls_variants_script,
            cohort_prefix = cohort_prefix,
            disk_size_gb = disk_size_gb,
            docker = docker,
            machine_mem_mb = machine_mem_mb

     }

     call Reformat_concordance_vcf_header {
         input:
            cohort_prefix = cohort_prefix,
            concordance_vcf = concordance_vcf,
            concordance_vcf_index = concordance_vcf_index,
            reformat_vcf_header_script = reformat_vcf_header_script,
            disk_size_gb = disk_size_gb,
            docker = docker,
            machine_mem_mb = machine_mem_mb

     }

     call Outlier_discovery_flag {
         input:
            update_outlier_discovery_flag_script = update_outlier_discovery_flag_script,
            outlier_size_range_variants_script = outlier_size_range_variants_script,
            final_update_outlier_discovery_flag_script = final_update_outlier_discovery_flag_script,
            reheadered_concordance_vcf = Reformat_concordance_vcf_header.reformatted_concordance_header_vcf,
            reheadered_concordance_vcf_index = Reformat_concordance_vcf_header.reformatted_concordance_header_vcf_index,
            concordance_final_outlier_variants = Get_outlier_variants.cohort_concordance_outlier_variants,
            cohort_prefix = cohort_prefix,
            outlier_samples_file = Determine_outlier_samples.outlier_samples_list,
            svtype_size_range_lower_cutoff = joinrawcalls_sv_counts_size_range_lower_cutoff,
            svtype_size_range_higher_cutoff = joinrawcalls_sv_counts_size_range_higher_cutoff,
            svtype = svtype,
            disk_size_gb = disk_size_gb,
            docker = docker,
            machine_mem_mb = machine_mem_mb
     }

     call Outlier_size_range_outlier_discovery_flagging {
         input:
            outlier_discovery_flagging_in_size_range_vcf = Outlier_discovery_flag.outlier_size_range_flagged_vcf,
            outlier_discovery_flagging_in_size_range_vcf_index = Outlier_discovery_flag.outlier_size_range_flagged_vcf_index,
            flag_outlier_variants_based_on_size_range_counts_script = flag_outlier_variants_based_on_size_range_counts_script,
            proportion_of_outlier_samples_associated_with_variant_script = proportion_of_outlier_samples_associated_with_variant_script,
            outlier_samples_file = Determine_outlier_samples.outlier_samples_list,
            cohort_prefix = cohort_prefix,
            disk_size_gb = disk_size_gb,
            svtype = svtype,
            svtype_size_range_lower_cutoff = joinrawcalls_sv_counts_size_range_lower_cutoff,
            svtype_size_range_higher_cutoff = joinrawcalls_sv_counts_size_range_higher_cutoff,
            fraction_of_outlier_samples = fraction_of_outlier_samples,
            docker = docker,
            machine_mem_mb = machine_mem_mb
     }

     call Remove_outlier_samples {
         input:
            final_flagged_vcf = Outlier_size_range_outlier_discovery_flagging.all_cohort_outlier_variants_flagged_vcf,
            final_flagged_vcf_index = Outlier_size_range_outlier_discovery_flagging.all_cohort_outlier_variants_flagged_vcf_index,
            outlier_samples_list = Determine_outlier_samples.outlier_samples_list,
            cohort_prefix = cohort_prefix,
            disk_size_gb = disk_size_gb,
            docker = docker,
            machine_mem_mb = machine_mem_mb
     }

     call ApplyManualFilter {
         input:
            cohort_prefix = cohort_prefix,
            vcf = Remove_outlier_samples.outlier_samples_removed_vcf,
            vcf_index = Remove_outlier_samples.outlier_samples_removed_vcf_index,
            filter_name = filter_name,
            bcftools_filter = bcftools_filter,
            disk_size_gb = disk_size_gb,
            docker = docker,
            machine_mem_mb = machine_mem_mb
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

    File make_tables_sql

    File duckdb_zip

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(join_raw_calls_vcf, 'GB') * 5.0)

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
      filtered.bcf > join_raw_calls_svlens.tsv

    bcftools query \
      --format '%ID\t%INFO/MEMBERS\n' \
      filtered.bcf \
      | awk -F'\t' '{split($2, a, /,/); for (i in a) print $1"\t"a[i]}' \
      > join_raw_calls_clusters.tsv

    unzip '~{duckdb_zip}'
    chmod u+x ./duckdb
    cat '~{make_tables_sql}' | ./duckdb joined_raw_calls.duckdb
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

    File templater_awk
    File count_svs_sql_tmpl

    File duckdb_zip

    String runtime_docker
  }

  Int disk_size_gb = ceil(size(joined_raw_calls_db, 'GB') + 8.0)

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

    printf 'svtype\t~{svtype}\n' > sql_params.tsv
    printf 'min_svlen\t~{min_svlin}\n' >> sql_params.tsv
    printf 'max_svlen\t~{max_svlen}\n' >> sql_params.tsv

    awk -f '~{templater_awk}' 'sql_params.tsv' '~{count_svs_sql_tmpl}' \
      | ./duckdb '~{joined_raw_calls_db}'
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
    bootDiskSizeGb: 15
    disks: 'local-disk ' + ceil(input_size * 1.3) + ' HDD'
    preemptible: 1
    maxRetries: 1
    docker: runtime_docker
  }

  output {
    File outlier_samples_list = '~{cohort_prefix}_outlier_sample_list.txt'
  }
}

task Concat_clusterbatch_vcfs {
    input {
        Array[File] clustered_depth_vcfs
        Array[File] clustered_wham_vcfs
        Array[File] clustered_manta_vcfs
        Array[File] clustered_melt_vcfs
        #Array[File] clusterbatch_vcf_indexes
        String cohort_prefix
        Int disk_size_gb
        String docker
        Int machine_mem_mb
    }
    
    Array[File] clusterbatch_vcfs = flatten([clustered_depth_vcfs, clustered_wham_vcfs, clustered_manta_vcfs, clustered_melt_vcfs])

    command <<<
        set -euo pipefail


        for vcf in ~{sep=' ' clusterbatch_vcfs}; do
   
            if [[ -f ${vcf} ]]; then
                echo "Reindexing VCF file: ${vcf}"
        
                if [[ -f "${vcf}.tbi" ]]; then
                    rm "${vcf}.tbi"
                fi
        
                tabix -p vcf ${vcf}
            else
                echo "VCF file not found: ${vcf}"
            fi
        done


        for vcf in ~{sep=' ' clusterbatch_vcfs}; do
            if [[ -f ${vcf} ]]; then
                echo "Converting VCF to BED: ${vcf}"
                svtk vcf2bed --include-filters -i ALL ${vcf} ${vcf}.bed
            else
                echo "VCF file not found: ${vcf}"
            fi
        done


        bed_files=$(for vcf in ~{sep=' ' clusterbatch_vcfs}; do echo "${vcf}.bed"; done)
        echo "Concatenating BED files: $bed_files"
        cat $bed_files > ~{cohort_prefix}_merged_clustered_vcfs.bed
        

    >>>

    runtime {
        memory: "~{machine_mem_mb} MiB"
        cpu: "1"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 1
        docker: docker
    }

    output {
        File merged_clusteredvcfs_bed = "~{cohort_prefix}_merged_clustered_vcfs.bed"
    }
}

task Get_outlier_variants {
    input {
        File outlier_samples_file
        File merged_clusteredvcfs_bedfile
        File join_raw_calls_bed
        File concordance_vcf
        File identify_outlier_variants_script
        File identify_join_raw_calls_variants_script
        Int disk_size_gb
        String cohort_prefix
        String docker
        Int machine_mem_mb
    }

    command <<<
        set -euo pipefail

        python3 ~{identify_outlier_variants_script} -o ~{outlier_samples_file} -v ~{merged_clusteredvcfs_bedfile} -out ~{cohort_prefix}_extracted_variants.tsv

        python3 ~{identify_join_raw_calls_variants_script} -e ~{cohort_prefix}_extracted_variants.tsv -b ~{join_raw_calls_bed} -out ~{cohort_prefix}_join_raw_calls_variants.txt

        svtk vcf2bed --include-filters -i ALL ~{concordance_vcf} ~{cohort_prefix}_concordance.bed

        cut -f1,2,3,4,35 ~{cohort_prefix}_concordance.bed > ~{cohort_prefix}_concordance_variants_and_truth_variants.bed

        cut -f35 ~{cohort_prefix}_concordance.bed > ~{cohort_prefix}_SV_concordance_variant_ids.txt

        cut -f4,5 ~{cohort_prefix}_concordance_variants_and_truth_variants.bed > concordance_filtered_variants_and_truth_variants.tsv

        comm -12 <(sort ~{cohort_prefix}_join_raw_calls_variants.txt) <(sort ~{cohort_prefix}_SV_concordance_variant_ids.txt) > ~{cohort_prefix}_filtered_outlier_variants.txt

        awk 'NR==FNR {ids[$1]; next} $2 in ids {print $1}' ~{cohort_prefix}_filtered_outlier_variants.txt concordance_filtered_variants_and_truth_variants.tsv > ~{cohort_prefix}_concordance_final_outlier_variant_ids.txt
    >>>

    runtime {
        memory: "~{machine_mem_mb} MiB"
        cpu: "1"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 1
        docker: docker
    }

    output {
        File cohort_concordance_outlier_variants = "~{cohort_prefix}_concordance_final_outlier_variant_ids.txt"
    }
}

task Reformat_concordance_vcf_header {
    input {
        String cohort_prefix
        File concordance_vcf
        File concordance_vcf_index
        File reformat_vcf_header_script
        Int disk_size_gb
        String docker
        Int machine_mem_mb
    }

    command <<<
        set -euo pipefail

        bcftools view -h ~{concordance_vcf} | grep '##FILTER' > ~{cohort_prefix}_concordance_header.txt

        echo '##FILTER=<ID=outlier_discovered,Description="Variant with only outlier samples associated">' >> ~{cohort_prefix}_concordance_header.txt

        python3 ~{reformat_vcf_header_script} -i ~{concordance_vcf} -hf ~{cohort_prefix}_concordance_header.txt -o ~{cohort_prefix}_reformatted_concordance_vcf_header.vcf.gz

        tabix -f ~{cohort_prefix}_reformatted_concordance_vcf_header.vcf.gz
    >>>

    runtime {
        memory: "~{machine_mem_mb} MiB"
        cpu: "1"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 1
        docker: docker
    }

    output {
        File reformatted_concordance_header_vcf = "~{cohort_prefix}_reformatted_concordance_vcf_header.vcf.gz"
        File reformatted_concordance_header_vcf_index = "~{cohort_prefix}_reformatted_concordance_vcf_header.vcf.gz.tbi"
    }
}

task Outlier_discovery_flag {
    input {
        File update_outlier_discovery_flag_script
        File outlier_size_range_variants_script
        File final_update_outlier_discovery_flag_script
        File reheadered_concordance_vcf
        File reheadered_concordance_vcf_index
        File concordance_final_outlier_variants
        String cohort_prefix
        File outlier_samples_file
        Int svtype_size_range_lower_cutoff
        Int svtype_size_range_higher_cutoff
        String svtype
        Int disk_size_gb
        String docker
        Int machine_mem_mb
    }

    command <<<
        set -euo pipefail

        python3 ~{update_outlier_discovery_flag_script} -i ~{reheadered_concordance_vcf} -v ~{concordance_final_outlier_variants} -o ~{cohort_prefix}_filtered_outlier_discovered_flag_added.vcf.gz
        tabix ~{cohort_prefix}_filtered_outlier_discovered_flag_added.vcf.gz

        svtk vcf2bed --include-filters -i ALL ~{cohort_prefix}_filtered_outlier_discovered_flag_added.vcf.gz ~{cohort_prefix}_filtered_outlier_discovered_flag_added.bed
        python3 ~{outlier_size_range_variants_script} -o ~{outlier_samples_file} -i ~{cohort_prefix}_filtered_outlier_discovered_flag_added.bed -out ~{cohort_prefix}_outlier_size_range_filtered_variants -t ~{svtype} -l ~{svtype_size_range_lower_cutoff} -hi ~{svtype_size_range_higher_cutoff}

        python3 ~{final_update_outlier_discovery_flag_script} -i ~{cohort_prefix}_filtered_outlier_discovered_flag_added.vcf.gz -o ~{cohort_prefix}_filtered_outlier_discovered_flag_added_FINAL.vcf.gz -f ~{cohort_prefix}_outlier_size_range_filtered_variants
        tabix ~{cohort_prefix}_filtered_outlier_discovered_flag_added_FINAL.vcf.gz
    >>>

    runtime {
        memory: "~{machine_mem_mb} MiB"
        cpu: "1"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 1
        docker: docker
    }

    output {
        File outlier_size_range_flagged_vcf = "~{cohort_prefix}_filtered_outlier_discovered_flag_added_FINAL.vcf.gz"
        File outlier_size_range_flagged_vcf_index = "~{cohort_prefix}_filtered_outlier_discovered_flag_added_FINAL.vcf.gz.tbi"
    }
}

task Outlier_size_range_outlier_discovery_flagging {
    input {
        File outlier_discovery_flagging_in_size_range_vcf
        File outlier_discovery_flagging_in_size_range_vcf_index
        File flag_outlier_variants_based_on_size_range_counts_script
        File proportion_of_outlier_samples_associated_with_variant_script
        File outlier_samples_file
        String cohort_prefix
        Int disk_size_gb
        String svtype
        Int svtype_size_range_lower_cutoff
        Int svtype_size_range_higher_cutoff
        Float fraction_of_outlier_samples
        String docker
        Int machine_mem_mb
    }

    command <<<
        set -euo pipefail

        svtk vcf2bed --include-filters -i ALL ~{outlier_discovery_flagging_in_size_range_vcf} ~{cohort_prefix}_filtered_outlier_discovered_flag_added_FINAL.bed

        awk -F'\t' '$5 == ~{svtype} && ($3 - $2) >= ~{svtype_size_range_lower_cutoff} && ($3 - $2) <= ~{svtype_size_range_higher_cutoff}' ~{cohort_prefix}_filtered_outlier_discovered_flag_added_FINAL.bed > ~{cohort_prefix}_filtered_outlier_discovered_flag_added_to_size_range_outliers_FINAL.bed

        python3 ~{proportion_of_outlier_samples_associated_with_variant_script} -o ~{outlier_samples_file} -i ~{cohort_prefix}_filtered_outlier_discovered_flag_added_to_size_range_outliers_FINAL.bed -out ~{cohort_prefix}_outlier_variants_in_size_range.tsv -f ~{fraction_of_outlier_samples}

        python3 ~{flag_outlier_variants_based_on_size_range_counts_script} -v ~{cohort_prefix}_outlier_variants_in_size_range.tsv -i ~{outlier_discovery_flagging_in_size_range_vcf} -o ~{cohort_prefix}_all_outlier_variants_flagged_with_outlier_discovered.vcf.gz

        tabix ~{cohort_prefix}_all_outlier_variants_flagged_with_outlier_discovered.vcf.gz
    >>>

    runtime {
        memory: "~{machine_mem_mb} MiB"
        cpu: "1"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 1
        docker: docker
    }

    output {
        File all_cohort_outlier_variants_flagged_vcf = "~{cohort_prefix}_all_outlier_variants_flagged_with_outlier_discovered.vcf.gz"
        File all_cohort_outlier_variants_flagged_vcf_index = "~{cohort_prefix}_all_outlier_variants_flagged_with_outlier_discovered.vcf.gz.tbi"
    }
}

task Remove_outlier_samples {
    input {
        File final_flagged_vcf
        File final_flagged_vcf_index
        File outlier_samples_list
        String cohort_prefix
        Int disk_size_gb
        String docker
        Int machine_mem_mb
    }

    command <<<
        set -euo pipefail

        bcftools view -c 1 --samples-file ^~{outlier_samples_list} -o ~{cohort_prefix}_FINAL_outlier_samples_removed.vcf.gz -O z ~{final_flagged_vcf}
        tabix ~{cohort_prefix}_FINAL_outlier_samples_removed.vcf.gz
    >>>

    runtime {
        memory: "~{machine_mem_mb} MiB"
        cpu: "1"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 1
        docker: docker
    }

    output {
        File outlier_samples_removed_vcf = "~{cohort_prefix}_FINAL_outlier_samples_removed.vcf.gz"
        File outlier_samples_removed_vcf_index = "~{cohort_prefix}_FINAL_outlier_samples_removed.vcf.gz.tbi"
    }
}

task ApplyManualFilter {
    input {
        String cohort_prefix
        File vcf
        File? vcf_index
        String filter_name
        String bcftools_filter
        Int disk_size_gb
        String docker
        Int machine_mem_mb
    }

    String hard_filtered_vcf_name = "~{cohort_prefix}.~{filter_name}.vcf.gz"

    command <<<
        set -euo pipefail

        bcftools view -e '~{bcftools_filter}' ~{vcf} -Oz -o "~{hard_filtered_vcf_name}"

        tabix "~{hard_filtered_vcf_name}"
    >>>

    runtime {
        memory: "~{machine_mem_mb} MiB"
        cpu: "1"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 1
        docker: docker
    }

    output {
        File hard_filtered_vcf = "~{hard_filtered_vcf_name}"
        File hard_filtered_vcf_index = "~{hard_filtered_vcf_name}.tbi"
    }
}


