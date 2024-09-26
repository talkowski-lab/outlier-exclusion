import argparse
import csv

def extract_outlier_variants(outlier_samples_file, input_bed_file, output_tsv_file, fraction_of_outlier_samples):
    
    outlier_samples = set()
    with open(outlier_samples_file, 'r') as f:
        for line in f:
            outlier_samples.add(line.strip())

    
    with open(input_bed_file, 'r') as infile, open(output_tsv_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        
        
        outfile.write("Variant_ID\tOutlier_Samples\n")

        
        for row in reader:
            variant_id = row[3]  
            samples = row[5].split(',')  
            
            
            outlier_only_samples = [sample for sample in samples if sample in outlier_samples]
            
            
            outlier_ratio = len(outlier_only_samples) / len(samples)
            
           
            if outlier_ratio >= fraction_of_outlier_samples:
                outfile.write(f"{variant_id}\t{','.join(outlier_only_samples)}\n")

    print(f"Finished processing {input_bed_file}. Extracted variants saved in {output_tsv_file}.")

if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description="Extract outlier variants from BED file.")
    
    parser.add_argument('-o', '--outlier_samples_file', required=True, help='Path to the file containing outlier sample IDs')
    parser.add_argument('-i', '--input_bed_file', required=True, help='Path to the input BED file')
    parser.add_argument('-out', '--output_tsv_file', required=True, help='Path to the output TSV file')
    parser.add_argument('-f', '--fraction_of_outlier_samples', required=True, help='Path to the output TSV file')

    args = parser.parse_args()

    
    extract_outlier_variants(args.outlier_samples_file, args.input_bed_file, args.output_tsv_file, args.fraction_of_outlier_samples)
