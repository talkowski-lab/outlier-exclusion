import argparse
import pysam


def filter_vcf_by_variant_ids(variant_ids_file, input_vcf, output_vcf):
    
    with open(variant_ids_file, "r") as f:
        variant_ids = set(line.strip() for line in f)

    
    with pysam.VariantFile(input_vcf, "r") as vcf_in, pysam.VariantFile(output_vcf, "w", header=vcf_in.header) as vcf_out:
        
        for record in vcf_in:
            variant_id = record.id
            
            
            if variant_id in variant_ids:
                
                if "." in record.filter:
                    record.filter.clear()  
                    record.filter.add('outlier_discovered')  
            
           
            vcf_out.write(record)

    print(f"Finished processing {input_vcf}. Filtered VCF saved as {output_vcf}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter VCF by variant IDs from a BED file.")
    
    parser.add_argument('-v', '--variant_ids_file', required=True, help='Path to the filtered BED file with variant IDs')
    parser.add_argument('-i', '--input_vcf', required=True, help='Path to the input VCF file')
    parser.add_argument('-o', '--output_vcf', required=True, help='Path to the output VCF file')

    args = parser.parse_args()

    filter_vcf_by_variant_ids(args.variant_ids_file, args.input_vcf, args.output_vcf)
