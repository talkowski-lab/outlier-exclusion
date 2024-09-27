import csv
import argparse

def filter_variants(outlier_samples_file, input_bed_file, output_variant_ids_file, svtype, svtype_size_range_lower_cutoff, svtype_size_range_higher_cutoff):
    
    outlier_samples = set()
    with open(outlier_samples_file, 'r') as f:
        for line in f:
            outlier_samples.add(line.strip())

    
    with open(input_bed_file, 'r') as infile, open(output_variant_ids_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')

        next(reader)  

        for row in reader:
            chrom, start, end, variant_id, svtype, samples = row[:6]

            try:
                start, end = int(start), int(end)
            except ValueError:
                continue

            size = end - start

            
            if svtype == svtype and svtype_size_range_lower_cutoff <= size <= svtype_size_range_higher_cutoff:
                sample_list = samples.split(',')
                outlier_only_samples = [sample for sample in sample_list if sample in outlier_samples]

                
                if len(outlier_only_samples) == len(sample_list):
                    outfile.write(f"{variant_id}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter variants with only outlier samples.")
    
    parser.add_argument('-o', '--outlier_samples_file', required=True, help='Path to the outlier samples file')
    parser.add_argument('-i', '--input_bed_file', required=True, help='Path to the input BED file')
    parser.add_argument('-out', '--output_variant_ids_file', required=True, help='Path to the output filtered variant IDs file')
    parser.add_argument('-t', '--svtype', required=True, help='svtype of variants causing outlier bumps')
    parser.add_argument('-l', '--svtype_size_range_lower_cutoff', required=True, type = int, help='lower cutoff of size range of svtype')
    parser.add_argument('-hi', '--svtype_size_range_higher_cutoff', required=True, type = int, help='higher cutoff of size range of svtype')




    args = parser.parse_args()

    filter_variants(args.outlier_samples_file, args.input_bed_file, args.output_variant_ids_file, args.svtype, args.svtype_size_range_lower_cutoff, args.svtype_size_range_higher_cutoff)
