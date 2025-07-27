import os
import subprocess
import glob

# Configuration
input_dir = "BRP"  # update this to your input folder
output_dir = "mutation_outputs"  # where to save .tsv files
fasta_file = "hg19.fa"
script_path = "scripts/mutation.py"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Match .vcf and .vcf.gz files
vcf_files = glob.glob(os.path.join(input_dir, "*.vcf")) + glob.glob(os.path.join(input_dir, "*.vcf.gz"))

for vcf_file in vcf_files:
    base = os.path.basename(vcf_file)
    base_no_ext = base.replace(".vcf.gz", "").replace(".vcf", "")
    output_file = os.path.join(output_dir, f"{base_no_ext}-mutation.tsv")

    print(f"Processing: {vcf_file} → {output_file}")

    cmd = ["python", script_path, vcf_file, fasta_file, output_file]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error processing {vcf_file}: {e}")
