import os
import subprocess
import glob
import pandas as pd

# === CONFIGURATION ===
input_dir = "BRP"  # folder containing .vcf or .vcf.gz
output_dir = "mutation_outputs"
fasta_file = "hg19.fa"
script_path = "scripts/mutation.py"
merged_output = "mutation_matrix.csv"

# === STEP 1: Run mutation.py on each VCF ===
os.makedirs(output_dir, exist_ok=True)
vcf_files = glob.glob(os.path.join(input_dir, "*.vcf")) + glob.glob(os.path.join(input_dir, "*.vcf.gz"))

for vcf_file in vcf_files:
    base = os.path.basename(vcf_file)
    sample_name = base.replace(".vcf.gz", "").replace(".vcf", "")
    output_file = os.path.join(output_dir, f"{sample_name}-mutation.tsv")

    print(f"Processing: {vcf_file} → {output_file}")

    cmd = ["python", script_path, vcf_file, fasta_file, output_file]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error processing {vcf_file}: {e}")

# === STEP 2: Build wide-format mutation matrix ===
all_files = glob.glob(os.path.join(output_dir, "*-mutation.tsv"))
matrix = pd.DataFrame()

for file in all_files:
    sample_name = os.path.basename(file).replace("-mutation.tsv", "")
    df = pd.read_csv(file, sep="\t")
    df = df.set_index("MutationType")["Count"]
    df.name = sample_name
    matrix = pd.concat([matrix, df], axis=1)

# Replace NaN with 0 (if a mutation is not present in a sample)
matrix = matrix.fillna(0).astype(int)

# === STEP 3: Save matrix ===
matrix.to_csv(merged_output)
print(f"\n✅ Wide-format mutation matrix saved to: {merged_output}")
