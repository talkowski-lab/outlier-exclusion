version 1.0

import "OutlierExclusion.wdl" as oe

workflow OutlierExclusionScattered {
  input {
    # Split by contig VCFs
    Array[File] join_raw_calls_vcfs

    Array[Array[String]] svtypes_to_filter = [["DEL", 5000, 25000], ["DUP", 5000, 25000]]

    File? wgd_scores
    Float? min_wgd_score
    Float? max_wgd_score
    Float iqr_multiplier = 8.0

    File? outlier_samples

    Array[File]? clustered_depth_vcfs
    Array[File]? clustered_manta_vcfs
    Array[File]? clustered_wham_vcfs
    Array[File]? clustered_melt_vcfs
    Array[File]? clustered_scramble_vcfs
    Float min_outlier_sample_prop

    Array[File] filter_genotypes_vcfs

    String base_docker
    String pipeline_docker

    # Same length as filter_genotypes_vcfs to give each output VCF a different
    # prefix
    Array[String]? output_prefix_list
    # File with one line per prefix to give each output VCF a different prefix
    File? output_prefix_file
    # Single prefix to use for all output VCFs
    String? output_prefix
  }

  scatter (vcf in join_raw_calls_vcfs) {
    call oe.GetJoinRawCallsClusters {
      input:
        vcf_or_bcf = vcf,
        base_docker = base_docker
    }

    if (!defined(outlier_samples)) {
      call oe.ConvertVcfOrBcfToTsv {
        input:
          vcf_or_bcf = vcf,
          base_docker = base_docker
      }
    }
  }

  call oe.MakeJoinRawCallsClustersDb {
    input:
      clusters = GetJoinRawCallsClusters.clusters,
      base_docker = base_docker
  }

  if (!defined(outlier_samples)) {
    call oe.MakeSvCountsDb {
      input:
        tsvs = select_all(ConvertVcfOrBcfToTsv.tsv),
        filters = svtypes_to_filter,
        pipeline_docker = pipeline_docker
    }

    call oe.DetermineOutlierSamples {
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
    call oe.FormatOutlierSamples {
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
    call oe.DetermineOutlierVariants {
      input:
        clustered_vcfs = batch_vcfs,
        outlier_samples_db = outlier_samples_db,
        min_outlier_sample_prop = min_outlier_sample_prop,
        jrc_clusters_db = MakeJoinRawCallsClustersDb.jrc_clusters_db,
        pipeline_docker = pipeline_docker
    }
  }

  Array[String]? prefix_list = if (!defined(output_prefix_list) && defined(output_prefix_file)) then read_lines(select_first([output_prefix_file])) else output_prefix_list
  String static_prefix = select_first([output_prefix, ""])

  scatter (i in range(length(filter_genotypes_vcfs))) {
    call oe.FlagOutlierVariants {
      input:
        output_prefix = if defined(prefix_list) then select_first([prefix_list])[i] else static_prefix,
        filter_genotypes_vcf = filter_genotypes_vcfs[i],
        outlier_variants = DetermineOutlierVariants.outlier_variants,
        pipeline_docker = pipeline_docker
    }
  }

  output {
    Array[File] outlier_annotated_vcf = FlagOutlierVariants.outlier_annotated_vcf
    Array[File] outlier_annotated_vcf_index = FlagOutlierVariants.outlier_annotated_vcf_index
    File? sv_count_outlier_samples = DetermineOutlierSamples.outlier_samples
    File? wgd_outlier_samples = DetermineOutlierSamples.wgd_outlier_samples
  }
}
