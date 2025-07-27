import gzip
import argparse

def add_vcf_header(input_file, output_file):
    with gzip.open(input_file, 'rt') as infile:
        lines = infile.readlines()

    # Check if header already contains the fileformat line
    header_line = "##fileformat=VCFv4.2\n"
    has_fileformat = any(line.startswith("##fileformat=") for line in lines)

    if not has_fileformat:
        # Insert at top before other ## headers
        i = 0
        while i < len(lines) and lines[i].startswith("##"):
            i += 1
        lines.insert(0, header_line)

    with gzip.open(output_file, 'wt') as outfile:
        outfile.writelines(lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ensure ##fileformat=VCFv4.2 is present in a gzipped VCF file.")
    parser.add_argument("input_file", help="Path to input .vcf.gz file")
    parser.add_argument("output_file", help="Path to output .vcf.gz file with corrected header")

    args = parser.parse_args()
    add_vcf_header(args.input_file, args.output_file)
