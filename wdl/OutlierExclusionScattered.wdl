version 1.0

import "OutlierExclusion.wdl" as oe

workflow OutlierExclusionScattered {
  input {
    # Split by contig VCFs
    Array[File] join_raw_calls_vcfs

    Array[Array[String]] svtypes_to_filter = [["DEL", 5000, 250000], ["DUP", 5000, 25000]]

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

    String docker
    String output_prefix
  }

  scatter (vcf in join_raw_calls_vcfs) {
    call oe.GetJoinRawCallsClusters {
      input:
        vcf_or_bcf = vcf,
        runtime_docker = docker
    }

    if (!defined(outlier_samples)) {
      call oe.ConvertVcfOrBcfToTsv {
        input:
          vcf_or_bcf = vcf,
          runtime_docker = docker
      }
    }
  }

  call oe.MakeJoinRawCallsClustersDb {
    input:
      clusters = GetJoinRawCallsClusters.clusters,
      runtime_docker = docker
  }

  if (defined(ConvertVcfOrBcfToTsv.tsv)) {
    call oe.MakeSvDb {
      input:
        tsvs = select_all(ConvertVcfOrBcfToTsv.tsv),
        runtime_docker = docker
    }

    call oe.MakeSvCountsDb {
      input:
        sv_db = MakeSvDb.sv_db,
        filters = svtypes_to_filter,
        runtime_docker = docker
    }

    call oe.DetermineOutlierSamples {
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
    call oe.FormatOutlierSamples {
      input:
        outlier_samples = select_first([outlier_samples]),
        runtime_docker = docker
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
        runtime_docker = docker
    }
  }

  scatter (vcf in filter_genotypes_vcfs) {
    call oe.FlagOutlierVariants {
      input:
        output_prefix = output_prefix,
        filter_genotypes_vcf = vcf,
        outlier_variants = DetermineOutlierVariants.outlier_variants,
        runtime_docker = docker
    }
  }

  output {
    Array[File] outlier_annotated_vcf = FlagOutlierVariants.outlier_annotated_vcf
    Array[File] outlier_annotated_vcf_index = FlagOutlierVariants.outlier_annotated_vcf_index
    File? sv_count_outlier_samples = DetermineOutlierSamples.outlier_samples
    File? wgd_outlier_samples = DetermineOutlierSamples.wgd_outlier_samples
  }
}
