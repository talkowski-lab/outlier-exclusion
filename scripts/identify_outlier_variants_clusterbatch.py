import csv
import argparse
import sys


csv.field_size_limit(sys.maxsize)

def extract_variants(outlier_samples_file, vcf_bed_file, output_file):
    
    outlier_samples = set()
    with open(outlier_samples_file, 'r') as f:
        for line in f:
            outlier_samples.add(line.strip())

    with open(vcf_bed_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')

        outfile.write("Variant_ID\tOutlier_Samples\n")

        for row in reader:
            variant_id = row[3]  
            samples = row[5].split(',')  

            outlier_only_samples = [sample for sample in samples if sample in outlier_samples]

            if len(outlier_only_samples) == len(samples):
                outfile.write(f"{variant_id}\t{','.join(outlier_only_samples)}\n")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Extract variants with all outlier samples.")
    parser.add_argument('-o', '--outlier_samples_file', required=True, help='Path to the outlier samples file')
    parser.add_argument('-v', '--vcf_bed_file', required=True, help='Path to the VCF BED file')
    parser.add_argument('-out', '--output_file', required=True, help='Path to the output TSV file')

    args = parser.parse_args()
    
    extract_variants(args.outlier_samples_file, args.vcf_bed_file, args.output_file)
