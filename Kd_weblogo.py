import os
import pandas as pd
import subprocess
import argparse

def generate_weblogo_from_csv(input_csv, seq_col, out_prefix, title):
    print(f"\n--- WEBLOGO GENERATION ---")
    
    # 1. Load the filtered CSV
    try:
        df = pd.read_csv(input_csv, sep=";")
    except Exception as e:
        print(f"[ERROR] Could not read {input_csv}: {e}")
        return
        
    # 2. Identify the sequence column
    actual_seq_col = None
    if seq_col in df.columns:
        actual_seq_col = seq_col
    else:
        # Auto-detect common column names for sequences
        possible_cols = [c for c in df.columns if any(x in c.lower() for x in ['secuencia', 'sequence', 'seq', 'rna', 'arn'])]
        if possible_cols:
            actual_seq_col = possible_cols[0]
            print(f"[*] Sequence column auto-detected as: '{actual_seq_col}'")
        else:
            print(f"[ERROR] No sequence column found. Available columns: {list(df.columns)}")
            return
    
    # 3. Create a temporary FASTA for WebLogo directly from CSV rows
    fasta_out_path = f"{out_prefix}.fasta"
    logo_path = f"{out_prefix}.png"
    
    valid_seqs = 0
    with open(fasta_out_path, 'w') as f_out:
        for idx, row in df.iterrows():
            seq = str(row[actual_seq_col]).strip()
            if seq != 'nan' and len(seq) > 0:
                # Create a simple ID for the FASTA entry based on the row index
                seq_id = f"Model_Top_{idx+1}"
                f_out.write(f">{seq_id}\n{seq}\n")
                valid_seqs += 1

    print(f"Extracted {valid_seqs} sequences directly from the CSV.")

    if valid_seqs == 0:
        print("[!] No valid sequences found. WebLogo aborted.")
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
    parser = argparse.ArgumentParser(description="Generates a WebLogo extracting sequences directly from a CSV.")
    parser.add_argument("-i", "--input", required=True, help="Filtered CSV file")
    parser.add_argument("-seq", "--seq_col", default="Secuencia", help="Name of the sequence column (default: Secuencia)")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output prefix for FASTA and PNG files")
    parser.add_argument("-t", "--title", required=True, help="WebLogo title")
    
    args = parser.parse_args()
    generate_weblogo_from_csv(args.input, args.seq_col, args.out_prefix, args.title)
