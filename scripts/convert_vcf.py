# convert_vcf.py

import argparse

def convert_vcf(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Write new VCF header
        outfile.write("##fileformat=VCFv4.2\n")
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for line in infile:
            if line.startswith("#"):
                continue  # skip original headers

            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue  # skip malformed lines

            chrom, pos, vid, ref, alt = parts[:5]
            # fill in dummy fields for QUAL, FILTER, INFO
            outfile.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat a VCF by replacing the header and standardizing fields.")
    parser.add_argument("input_file", help="Path to input VCF file")
    parser.add_argument("output_file", help="Path to output VCF file")

    args = parser.parse_args()
    convert_vcf(args.input_file, args.output_file)
