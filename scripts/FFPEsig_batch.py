#!/user/bin/python

## loading packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import getopt
from datetime import datetime, date
from pathlib import Path

# Import functions from the original FFPEsig.py
from FFPEsig import sig_extraction, correct_FFPE_profile, SBS96_plot, ffpe_sig_repaired, ffpe_sig_unrepaired

def process_single_sample(df, sample_id, label, output_dir, ite=100, precision=0.95):
    """
    Process a single sample with FFPE signature correction
    """
    print(f"Processing sample: {sample_id}")
    
    # Get the sample data
    sample_data = df[sample_id].to_numpy()
    
    # Choose the appropriate FFPE signature
    if label == "Unrepaired":
        ffpe_sig = ffpe_sig_unrepaired
    elif label == "Repaired":
        ffpe_sig = ffpe_sig_repaired
    else:
        raise ValueError("Label must be 'Unrepaired' or 'Repaired'")
    
    # Correct the profile
    corrected_profile, corrected_solutions = correct_FFPE_profile(
        V=sample_data,
        W1=ffpe_sig,
        sample_id=sample_id,
        ite=ite,
        precision=precision
    )
    
    # Create output filenames
    output1 = os.path.join(output_dir, f"{sample_id}_corrected_profile.csv")
    output2 = os.path.join(output_dir, f"{sample_id}_all_solutions.csv")
    plot_bef = os.path.join(output_dir, f"{sample_id}_before_correction.pdf")
    plot_aft = os.path.join(output_dir, f"{sample_id}_after_correction.pdf")
    
    # Save corrected profile
    df_tmp = pd.DataFrame()
    df_tmp[sample_id] = corrected_profile
    df_tmp.to_csv(output1, encoding='utf-8', index=False)
    
    # Save all solutions
    corrected_solutions.to_csv(output2, encoding='utf-8', index=False)
    
    # Create plots
    SBS96_plot(sample_data, label=f"{sample_id}\n", 
               name="Original", file=plot_bef,
               height=2, width=10, s=11)
    
    SBS96_plot(corrected_profile, label=f"{sample_id}\n", 
               name="Corrected", file=plot_aft,
               height=2, width=10, s=11)
    
    print(f"✓ Completed {sample_id}")
    return corrected_profile

def process_all_samples(input_file, label, output_dir, ite=100, precision=0.95):
    """
    Process all samples in the mutation matrix
    """
    print(f"Loading mutation matrix: {input_file}")
    df = pd.read_csv(input_file, index_col=0)
    
    # Get all sample columns (exclude the index column)
    sample_columns = df.columns.tolist()
    print(f"Found {len(sample_columns)} samples: {sample_columns}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each sample
    corrected_profiles = {}
    for sample_id in sample_columns:
        try:
            corrected_profile = process_single_sample(
                df, sample_id, label, output_dir, ite, precision
            )
            corrected_profiles[sample_id] = corrected_profile
        except Exception as e:
            print(f"❌ Error processing {sample_id}: {e}")
            continue
    
    # Create combined corrected matrix
    if corrected_profiles:
        corrected_df = pd.DataFrame(corrected_profiles)
        combined_output = os.path.join(output_dir, "all_samples_corrected_matrix.csv")
        corrected_df.to_csv(combined_output, index=True)
        print(f"✓ Combined corrected matrix saved to: {combined_output}")
    
    return corrected_profiles

def main():
    today = date.today()
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(f"Start FFPEsig Batch Processing at {current_time} on {today}")

    usage = """    
FFPEsig Batch Processing - Process all samples in a mutation matrix
------------------------
To run FFPEsig batch processing:

    python FFPEsig_batch.py [--input|-i] <Path-to-mutation-matrix> [--label|-l] <Unrepaired|Repaired> [--output_dir|-o] <Path-of-output-folder>

------------------------    
Example:

    python FFPEsig_batch.py --input results/mutation_matrix.csv --label Unrepaired --output_dir mutation_corrected_outputs_unrepaired

    OR
    
    python FFPEsig_batch.py -i results/mutation_matrix.csv -l Unrepaired -o mutation_corrected_outputs_unrepaired
------------------------    
Note:

    1. Input file, [--input|-i], must be a CSV format mutation matrix with sample IDs as columns;
    2. Label option, [--label|-l], must be either <Unrepaired|Repaired>;
    3. All samples in the matrix will be processed automatically.
    
    """
    
    # Get full command-line arguments
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    
    short_options = "h:i:l:o:t:p:"
    long_options = ["help", "input=", "label=", "output_dir=", "iteration=", "precision="]
    
    try:
        arguments, values = getopt.getopt(argument_list, short_options, long_options)
    except getopt.error as err:
        print(str(err))
        sys.exit(2)
        
    arg_length = 0
    ite = 100
    prec = 0.95
    file = None
    label = None
    out_dir = None
    
    # Evaluate given options
    for current_argument, current_value in arguments:
        if current_argument in ("-h", "--help"):
            print(usage)
            sys.exit(2)
        elif current_argument in ("-i", "--input"):
            file = current_value
            arg_length += 1
            print(f'\tInput file: {file}')
        elif current_argument in ("-l", "--label"):
            label = current_value
            arg_length += 1
            print(f"\tCorrecting {label} FFPE samples")
        elif current_argument in ("-o", "--output_dir"):
            out_dir = current_value
            arg_length += 1
            print(f"\tOutput folder: {out_dir}")
        elif current_argument in ("-t", "--iteration"):
            ite = int(current_value)
            print(f"\tIterations: {ite}")
        elif current_argument in ("-p", "--precision"):
            prec = float(current_value)
            print(f"\tPrecision: {prec}")

    if arg_length != 3:
        print("Three correct arguments are required to run FFPEsig batch processing.")
        print(usage)
        sys.exit(2)
    
    # Process all samples
    corrected_profiles = process_all_samples(file, label, out_dir, ite, prec)
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(f'Batch processing completed at {current_time} on {today}')
    print(f'Processed {len(corrected_profiles)} samples successfully')

if __name__ == "__main__":
    main() 