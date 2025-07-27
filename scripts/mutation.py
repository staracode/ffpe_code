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

def count_mutations(vcf_file, fasta_file, output_file):
    fasta = Fasta(fasta_file)
    fixed_vcf = fix_vcf_header_if_needed(vcf_file)
    vcf = VCF(fixed_vcf)
    mutation_counts = {key: 0 for key in context_keys}

    for variant in vcf:
        if len(variant.REF) != 1 or any(len(alt) != 1 for alt in variant.ALT):
            continue  # skip non-SNVs

        for alt in variant.ALT:
            key = normalize_context(variant.CHROM, variant.POS, variant.REF, alt, fasta)
            if key and key in mutation_counts:
                mutation_counts[key] += 1

    with open(output_file, "w") as out:
        out.write("MutationType\tCount\n")
        for key in sorted(mutation_counts.keys()):
            out.write(f"{key}\t{mutation_counts[key]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count 96 trinucleotide contexts from a VCF.")
    parser.add_argument("input_file", help="Input VCF file (SNVs only)")
    parser.add_argument("reference_fasta", help="Indexed reference FASTA file")
    parser.add_argument("output_file", help="Output TSV file for mutation counts")

    args = parser.parse_args()
    count_mutations(args.input_file, args.reference_fasta, args.output_file)
