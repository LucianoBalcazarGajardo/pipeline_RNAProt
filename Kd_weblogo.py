import os
import re
import pandas as pd
import subprocess
import argparse
from Bio import SeqIO

def generate_weblogo(input_csv, master_fasta, out_prefix, title):
    print(f"\n--- WEBLOGO GENERATION ---")
    
    # 1. Load the filtered CSV (Using ";" as separator matching your Kd scripts)
    try:
        df = pd.read_csv(input_csv, sep=";")
    except Exception as e:
        print(f"[ERROR] Could not read {input_csv}: {e}")
        return
        
    # 2. Load Master FASTA
    rna_dict = {}
    try:
        for record in SeqIO.parse(master_fasta, "fasta"):
            rna_dict[record.id] = str(record.seq)
    except FileNotFoundError:
        print(f"[ERROR] Master FASTA file not found at {master_fasta}")
        return
        
    # 3. Extract ARN IDs and create a temporary FASTA for WebLogo
    fasta_out_path = f"{out_prefix}.fasta"
    logo_path = f"{out_prefix}.png"
    
    matches_found = 0
    with open(fasta_out_path, 'w') as f_out:
        for idx, row in df.iterrows():
            # Convert row to string to robustly search for "ARN" followed by numbers
            row_str = " ".join(row.astype(str))
            match = re.search(r'(ARN\d+)', row_str)
            if match:
                arn_id = match.group(1)
                if arn_id in rna_dict:
                    f_out.write(f">{arn_id}\n{rna_dict[arn_id]}\n")
                    matches_found += 1

    print(f"Matched {matches_found} RNA sequences from the CSV.")

    if matches_found == 0:
        print("[!] No RNA IDs matched with the Master FASTA. WebLogo aborted.")
        return

    # 4. Execute WebLogo3
    print("Executing WebLogo3...")
    try:
        command = [
            "weblogo", 
            "-f", fasta_out_path, 
            "-o", logo_path, 
            "-F", "png", 
            "--resolution", "300", 
            "--title", title,
            "--title-fontsize", "10",
            "--size", "large",
            "--units", "probability",    
            "--sequence-type", "rna"     
        ]
        subprocess.run(command, check=True)
        print(f"[OK] WebLogo saved at: {logo_path}\n")
    except Exception as e:
        print(f"[ERROR] Executing WebLogo: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates a WebLogo from a filtered Kd CSV.")
    parser.add_argument("-i", "--input", required=True, help="Filtered CSV file")
    parser.add_argument("-fasta", "--master_fasta", required=True, help="Path to Master FASTA file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output prefix for FASTA and PNG files (without extension)")
    parser.add_argument("-t", "--title", required=True, help="WebLogo title")
    
    args = parser.parse_args()
    generate_weblogo(args.input, args.master_fasta, args.out_prefix, args.title)
