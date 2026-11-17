import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def analyze_percentile(input_csv, img_out, csv_out, txt_out, percentile, title):
    # 1. Load CSV and find Kd column
    try:
        df = pd.read_csv(input_csv, sep=";")
    except Exception as e:
        print(f"[ERROR] Could not read {input_csv}: {e}")
        return

    kd_cols = [col for col in df.columns if 'Kd' in col or 'kd' in col.lower()]
    if not kd_cols:
        print("[ERROR] No Kd column found in the CSV.")
        return
    target_col = kd_cols[0]

    # 2. Filter errors and convert to numeric
    total_models = len(df)
    df_valid = df[~df[target_col].astype(str).str.contains("ERROR", na=False)].copy()
    df_valid[target_col] = df_valid[target_col].astype(float)
    valid_models = len(df_valid)

    if valid_models == 0:
        print("[!] No viable models to analyze.")
        return

    # 3. Sort and extract the percentile (Top X%)
    df_sorted = df_valid.sort_values(by=target_col, ascending=True)
    num_top = max(1, int((percentile / 100.0) * valid_models))
    df_top = df_sorted.head(num_top)
    
    kd_cutoff = df_top[target_col].max()

    # Save CSV with filtered models
    os.makedirs(os.path.dirname(os.path.abspath(csv_out)), exist_ok=True)
    df_top.to_csv(csv_out, sep=";", index=False)

    # 4. Prepare and save statistics in .txt
    info_text = (
        f"--- PERCENTILE ANALYSIS (Top {percentile}%) ---\n"
        f"Total source models: {total_models}\n"
        f"Viable models:       {valid_models}\n"
        f"Kd Cutoff (Limit):   {kd_cutoff:.2e} M\n"
        f"Extracted models:    {num_top}\n"
    )
    print(f"\n{info_text}")
    
    os.makedirs(os.path.dirname(os.path.abspath(txt_out)), exist_ok=True)
    with open(txt_out, "w") as f:
        f.write(info_text)

    # 5. Generate Logarithmic Histogram
    plt.figure(figsize=(12, 7))
    min_kd, max_kd = df_valid[target_col].min(), df_valid[target_col].max()
    if min_kd <= 0: min_kd = 1e-15
    bins = np.logspace(np.log10(min_kd), np.log10(max_kd), 40)
    
    plt.hist(df_valid[target_col], bins=bins, color='lightgray', edgecolor='black', alpha=0.7)
    plt.xscale('log')
    
    # Demarcate and shade the percentile region
    plt.axvline(kd_cutoff, color='red', linestyle='-', linewidth=2.5, 
                label=f'Top {percentile}% Cutoff ($\leq$ {kd_cutoff:.1e} M)')
    plt.axvspan(min_kd, kd_cutoff, color='red', alpha=0.2, label='Selected Region')

    # Customization
    plt.title(title, fontsize=15, fontweight='bold')
    plt.xlabel('Dissociation Constant ($K_d$) in Molar [Log Scale]', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.grid(axis='y', alpha=0.5, linestyle='--')
    plt.legend(loc='upper right')

    plt.savefig(img_out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[OK] Plot saved at: {img_out}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extracts and plots the Top X% of models by Kd.")
    parser.add_argument("-i", "--input", required=True, help="Original CSV file")
    parser.add_argument("-o", "--img_out", required=True, help="Output path for PNG image")
    parser.add_argument("-c", "--csv_out", required=True, help="Output path for the new filtered CSV")
    parser.add_argument("-txt", "--txt_out", required=True, help="Output path for the statistics TXT file")
    parser.add_argument("-p", "--percentile", type=float, default=5.0, help="Percentile to extract (e.g., 5 for Top 5%)")
    parser.add_argument("-t", "--title", default="$K_d$ Percentile Analysis", help="Plot title")
    
    args = parser.parse_args()
    analyze_percentile(args.input, args.img_out, args.csv_out, args.txt_out, args.percentile, args.title)
