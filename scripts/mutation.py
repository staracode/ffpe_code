import argparse
from cyvcf2 import VCF
from pyfaidx import Fasta
from collections import defaultdict
from Bio.Seq import reverse_complement
import tempfile
import shutil
import gzip

# Initialize 96 trinucleotide contexts
bases = ['A', 'C', 'G', 'T']
mutations = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
context_keys = [f"{m}@{x}_{y}" for m in mutations for x in bases for y in bases]

def fix_vcf_header_if_needed(vcf_path):
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, 'rt') as f:
        lines = f.readlines()

    if any(line.startswith("##fileformat=") for line in lines):
        return vcf_path  # already fine

    # Insert fileformat header
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz", mode="wt")
    temp.write("##fileformat=VCFv4.2\n")
    for line in lines:
        if not line.startswith("##fileformat="):
            temp.write(line)
    temp.close()
    return temp.name

def normalize_context(chrom, pos, ref, alt, fasta):
    """Return pyrimidine-normalized trinucleotide context and mutation key."""
    # Normalize chromosome name for fasta access
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom

    try:
        seq = fasta[chrom][pos - 2:pos + 1].seq.upper()  # 1-based indexing
    except KeyError:
        # chrom not found in fasta
        return None
    except:
        return None

    if len(seq) != 3:
        return None

    upstream, base, downstream = seq[0], seq[1], seq[2]

    if base != ref:
        return None

    if ref in ['C', 'T']:
        mut = f"{ref}>{alt}"
        return f"{mut}@{upstream}_{downstream}"
    else:
        ref_rc = reverse_complement(ref)
        alt_rc = reverse_complement(alt)
        mut = f"{ref_rc}>{alt_rc}"
        new_seq = reverse_complement(seq)
        return f"{mut}@{new_seq[0]}_{new_seq[2]}"

def count_mutations(vcf_file, fasta_file, output_file, format_type="tsv", sort_by="mutation", verbose=False):
    if verbose:
        print("Loading reference genome...")
    
    fasta = Fasta(fasta_file)
    fixed_vcf = fix_vcf_header_if_needed(vcf_file)
    vcf = VCF(fixed_vcf)
    mutation_counts = {key: 0 for key in context_keys}

    if verbose:
        print("Counting mutations...")
        variant_count = 0

    for variant in vcf:
        if len(variant.REF) != 1 or any(len(alt) != 1 for alt in variant.ALT):
            continue  # skip non-SNVs

        for alt in variant.ALT:
            key = normalize_context(variant.CHROM, variant.POS, variant.REF, alt, fasta)
            if key and key in mutation_counts:
                mutation_counts[key] += 1
        
        if verbose:
            variant_count += 1
            if variant_count % 1000 == 0:
                print(f"Processed {variant_count} variants...")

    if verbose:
        print(f"Total variants processed: {variant_count}")

    # Sort based on user preference
    if sort_by == "count":
        sorted_items = sorted(mutation_counts.items(), key=lambda x: x[1], reverse=True)
    else:  # sort_by == "mutation"
        sorted_items = sorted(mutation_counts.items())

    # Determine separator and header based on format
    if format_type == "csv":
        separator = ","
        header = "MutationType,Count\n"
    else:  # tsv
        separator = "\t"
        header = "MutationType\tCount\n"

    with open(output_file, "w") as out:
        out.write(header)
        for key, count in sorted_items:
            out.write(f"{key}{separator}{count}\n")

    if verbose:
        print(f"Results written to: {output_file}")
        print(f"Total mutation types: {len(mutation_counts)}")
        print(f"Total mutations counted: {sum(mutation_counts.values())}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count 96 trinucleotide contexts from a VCF file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python mutation.py --input sample.vcf --reference hg19.fa --output counts.tsv
  python mutation.py -i sample.vcf.gz -r hg19.fa -o counts.tsv --verbose
  python mutation.py --input sample.vcf --reference hg19.fa --output counts.tsv --format csv
        """
    )
    
    # Required arguments
    parser.add_argument("--input", "-i", 
                       required=True,
                       help="Input VCF file (SNVs only, can be gzipped)")
    parser.add_argument("--reference", "-r", 
                       required=True,
                       help="Indexed reference FASTA file")
    parser.add_argument("--output", "-o", 
                       required=True,
                       help="Output file for mutation counts")
    
    # Optional arguments
    parser.add_argument("--format", "-f",
                       choices=["tsv", "csv"],
                       default="tsv",
                       help="Output format (default: tsv)")
    parser.add_argument("--verbose", "-v",
                       action="store_true",
                       help="Enable verbose output")
    parser.add_argument("--sort-by", "-s",
                       choices=["mutation", "count"],
                       default="mutation",
                       help="Sort output by mutation type or count (default: mutation)")

    args = parser.parse_args()
    
    if args.verbose:
        print(f"Processing VCF file: {args.input}")
        print(f"Using reference: {args.reference}")
        print(f"Output file: {args.output}")
        print(f"Output format: {args.format}")
        print(f"Sort by: {args.sort_by}")
        print("-" * 50)
    
    count_mutations(args.input, args.reference, args.output)
