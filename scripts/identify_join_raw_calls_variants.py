import csv
import argparse

def extract_variant_ids(extracted_variants_file, bed_file, output_file):
   
    variant_ids_set = set()
    with open(extracted_variants_file, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            variant_ids_set.add(row[0]) 

   
    with open(bed_file, 'r') as bedfile, open(output_file, 'w') as outfile:
        reader = csv.reader(bedfile, delimiter='\t')

        
        for row in reader:
            variant_id = row[3]  
            members = row[10].split(',')  
            
           
            if any(member in variant_ids_set for member in members):
                outfile.write(f"{variant_id}\n")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Extract variant IDs from BED based on extracted variants.")
    parser.add_argument('-e', '--extracted_variants_file', required=True, help='Path to the extracted variants TSV file')
    parser.add_argument('-b', '--bed_file', required=True, help='Path to the BED file containing variants and samples')
    parser.add_argument('-out', '--output_file', required=True, help='Path to the output file where matching variant IDs will be saved')

    
    args = parser.parse_args()

    
    extract_variant_ids(args.extracted_variants_file, args.bed_file, args.output_file)
