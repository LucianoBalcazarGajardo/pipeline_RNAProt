import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def plot_kd(input_csv, img_out, title):
    # 1. Load CSV
    try:
        df = pd.read_csv(input_csv, sep=";")
    except Exception as e:
        print(f"[ERROR] Could not read {input_csv}: {e}")
        return

    # Find Kd column
    kd_cols = [col for col in df.columns if 'Kd' in col or 'kd' in col.lower()]
    if not kd_cols:
        print("[ERROR] No Kd column found in the CSV.")
        return
    target_col = kd_cols[0]

    # 2. FILTER ERRORS AND CONVERT TO NUMERIC
    total_models = len(df)
    df_valid = df[~df[target_col].astype(str).str.contains("ERROR", na=False)].copy()
    df_valid[target_col] = df_valid[target_col].astype(float)
    valid_models = len(df_valid)

    if valid_models == 0:
        print("[!] No viable models to plot.")
        return

    print(f"\n--- Kd DISTRIBUTION ---")
    print(f"Total models processed: {total_models}")
    print(f"Plotted models:         {valid_models}")
    print(f"Best absolute Kd (Min): {df_valid[target_col].min():.2e} M")
    print(f"Worst absolute Kd (Max):{df_valid[target_col].max():.2e} M\n")

    # 3. GENERATE RAW LOGARITHMIC HISTOGRAM (NO LINES)
    plt.figure(figsize=(12, 7))
    
    min_kd, max_kd = df_valid[target_col].min(), df_valid[target_col].max()
    if min_kd <= 0: min_kd = 1e-15 # Prevent log(0)
    bins = np.logspace(np.log10(min_kd), np.log10(max_kd), 40)
    
    plt.hist(df_valid[target_col], bins=bins, color='mediumseagreen', edgecolor='black', alpha=0.8)
    plt.xscale('log')
    
    # Basic customization
    plt.title(title, fontsize=15, fontweight='bold')
    plt.xlabel('Dissociation Constant ($K_d$) in Molar [Log Scale]', fontsize=12)
    plt.ylabel('Number of Models (Frequency)', fontsize=12)
    plt.grid(axis='y', alpha=0.5, linestyle='--')

    # Save image
    os.makedirs(os.path.dirname(os.path.abspath(img_out)), exist_ok=True)
    plt.savefig(img_out, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[OK] Raw histogram saved at: {img_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Raw Kd histogram generator.")
    parser.add_argument("-i", "--input", required=True, help="Path to original CSV file")
    parser.add_argument("-o", "--img_out", required=True, help="Output path for PNG image")
    parser.add_argument("-t", "--title", default="Raw $K_d$ Distribution", help="Plot title")
    
    args = parser.parse_args()
    plot_raw_kd(args.input, args.img_out, args.title)
