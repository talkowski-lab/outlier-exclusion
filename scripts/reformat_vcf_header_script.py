import csv
import pysam
import argparse

def update_vcf_with_filter_status(original_vcf, new_header_file, output_vcf):
    
    vcf_in = pysam.VariantFile(original_vcf)

    
    with open(new_header_file, 'r') as f:
        new_header_lines = f.readlines()

    
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

    
    for line in new_header_lines:
        if line.startswith("##FILTER"):
            vcf_out.header.add_line(line.strip())

    
    for record in vcf_in:
        vcf_out.write(record)

    
    vcf_in.close()
    vcf_out.close()

    
    pysam.tabix_index(output_vcf, preset="vcf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update VCF with new FILTER header lines and reindex.")
    parser.add_argument('-i', '--input_vcf', required=True, help="Path to the input VCF file")
    parser.add_argument('-hf', '--header_file', required=True, help="Path to the file containing the new header lines")
    parser.add_argument('-o', '--output_vcf', required=True, help="Path to the output VCF file")

    args = parser.parse_args()

    update_vcf_with_filter_status(args.input_vcf, args.header_file, args.output_vcf)
