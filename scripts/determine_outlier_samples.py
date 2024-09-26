import pandas as pd
import numpy as np
import argparse

def detect_outliers(sv_counts_file, iqr_multiplier, sv_counts_specific_size_range_file, wgd_score_file, wgd_low_cutoff, wgd_high_cutoff, output_file):
    sample_sv_count_jrc = pd.read_csv(sv_counts_file, sep='\t', names=["Samples", "SVTYPE", "Count"])

    pivot_df_jrc = sample_sv_count_jrc.pivot_table(index='Samples', columns='SVTYPE', values='Count', fill_value=0)
    pivot_df_jrc = pivot_df_jrc.reset_index()

    pivot_df_jrc['Sum'] = pivot_df_jrc.iloc[:, 1:5].sum(axis=1)

    median_total_jrc = np.median(pivot_df_jrc['Sum'])
    median_del_jrc = np.median(pivot_df_jrc['DEL'])
    median_dup_jrc = np.median(pivot_df_jrc['DUP'])

    Q1_total_jrc = pivot_df_jrc['Sum'].quantile(0.25)
    Q3_total_jrc = pivot_df_jrc['Sum'].quantile(0.75)
    IQR_total_jrc = Q3_total_jrc - Q1_total_jrc
    IQR_6_total_jrc = iqr_multiplier * IQR_total_jrc
    outlier_cutoff_jrc = median_total_jrc + IQR_6_total_jrc

    outliers_jrc = pivot_df_jrc[pivot_df_jrc['Sum'] > outlier_cutoff_jrc]

    Q1_del_jrc = pivot_df_jrc['DEL'].quantile(0.25)
    Q3_del_jrc = pivot_df_jrc['DEL'].quantile(0.75)
    IQR_del_jrc = Q3_del_jrc - Q1_del_jrc
    IQR_6_del_jrc = iqr_multiplier * IQR_del_jrc
    outlier_cutoff_del_jrc = median_del_jrc + IQR_6_del_jrc
    outliers_del_jrc = pivot_df_jrc[pivot_df_jrc['DEL'] > outlier_cutoff_del_jrc]

    Q1_dup_jrc = pivot_df_jrc['DUP'].quantile(0.25)
    Q3_dup_jrc = pivot_df_jrc['DUP'].quantile(0.75)
    IQR_dup_jrc = Q3_dup_jrc - Q1_dup_jrc
    IQR_6_dup_jrc = iqr_multiplier * IQR_dup_jrc
    outlier_cutoff_dup_jrc = median_dup_jrc + IQR_6_dup_jrc
    outliers_dup_jrc = pivot_df_jrc[pivot_df_jrc['DUP'] > outlier_cutoff_dup_jrc]

    outlier_samples_total_jrc = list(outliers_jrc['Samples'])
    outlier_samples_del_jrc = list(outliers_del_jrc['Samples'])
    outlier_samples_dup_jrc = list(outliers_dup_jrc['Samples'])

    specific_size_range = pd.read_csv(sv_counts_specific_size_range_file, sep='\t', names=["Count", "Samples"])
    median_specific_size_range = np.median(specific_size_range['Count'])

    Q1_svtype_size_range = specific_size_range['Count'].quantile(0.25)
    Q3_svtype_size_range = specific_size_range['Count'].quantile(0.75)
    IQR_svtype_specific_size_range = Q3_svtype_size_range - Q1_svtype_size_range
    IQR_outlier_criteria = iqr_multiplier * IQR_svtype_specific_size_range

    outlier_cutoff_specific_size_range = median_specific_size_range + IQR_outlier_criteria
    outliers_specific_size_range = specific_size_range[specific_size_range['Count'] > outlier_cutoff_specific_size_range]
    outlier_specific_size_range_samples = list(outliers_specific_size_range['Samples'])
    
    wgd_scores = pd.read_csv(wgd_score_file, sep='\t')
    wgd_scores['score'] = pd.to_numeric(wgd_scores['score'], errors='coerce')
    
    low_cutoff = wgd_scores[wgd_scores['score'] < wgd_low_cutoff]
    high_cutoff = wgd_scores[wgd_scores['score'] > wgd_high_cutoff]
    combine_wgd_outliers = pd.concat([low_cutoff, high_cutoff])
    combine_wgd_outliers_samples = list(combine_wgd_outliers['#ID'])

    merge_large_dups_and_wgd = list(set(outlier_samples_total_jrc + outlier_samples_del_jrc + outlier_samples_dup_jrc + outlier_specific_size_range_samples + combine_wgd_outliers_samples))

    merged_outlier_list_jrc = pd.DataFrame(merge_large_dups_and_wgd, columns=['Sample'])
    merged_outlier_list_jrc.to_csv(output_file, index=False, header=False)

    print(f"Outliers saved to 'merged_outlier_samples.txt' with IQR multiplier {iqr_multiplier}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect outliers in SV count files.")
    parser.add_argument('-s', '--sv_counts_file', required=True, help='Path to the input SV counts file')
    parser.add_argument('-r', '--sv_counts_specific_size_range_file', required=True, help='Path to the dup_5_to_10kb counts file')
    parser.add_argument('-i', '--iqr_multiplier', type=float, default=8.0, help='IQR multiplier for outlier detection (default: 8.0)')
    parser.add_argument('-o', '--output_file', type=str, help='outlier samples textfile')
    parser.add_argument('-w', '--wgd_score_file', type=str, help='wgd score file')
    parser.add_argument('-l', '--wgd_lower_cutoff', type=float, help='wgd score lower cutoff')
    parser.add_argument('-hi', '--wgd_higher_cutoff', type=float, help='wgd score higher cutoff')

    args = parser.parse_args()

    detect_outliers(args.sv_counts_file, args.iqr_multiplier, args.sv_counts_specific_size_range_file, args.wgd_score_file, args.wgd_lower_cutoff, args.wgd_higher_cutoff, args.output_file)
