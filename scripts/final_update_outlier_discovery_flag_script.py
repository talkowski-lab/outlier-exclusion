import argparse
import pysam

def filter_variants(input_vcf, output_vcf, filtered_variant_ids_file):
    
    with open(filtered_variant_ids_file, "r") as f:
        variant_ids = set(line.strip() for line in f)

    
    with pysam.VariantFile(input_vcf, "r") as vcf_in, pysam.VariantFile(output_vcf, "w", header=vcf_in.header) as vcf_out:
        for record in vcf_in:
            variant_id = record.id

            
            if variant_id in variant_ids:
                if "." in record.filter:
                    record.filter.clear() 
                    record.filter.add('outlier_discovered')  

           
            vcf_out.write(record)

    print(f"Finished replacing PASS (.) with 'outlier_discovered' for matching variants in {output_vcf}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter variants and update VCF FILTER field for matching variant IDs.")
    parser.add_argument('-i', '--input_vcf', required=True, help="Path to the input VCF file (gzipped)")
    parser.add_argument('-o', '--output_vcf', required=True, help="Path to the output VCF file (gzipped)")
    parser.add_argument('-f', '--filtered_variant_ids_file', required=True, help="Path to the file containing filtered variant IDs")

    args = parser.parse_args()

    
    filter_variants(args.input_vcf, args.output_vcf, args.filtered_variant_ids_file)
