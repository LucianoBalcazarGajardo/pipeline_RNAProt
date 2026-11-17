import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def analyze_kd(input_csv, img_out, csv_out, top_percent, thresholds, title):
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
    df_valido = df[~df[target_col].astype(str).str.contains("ERROR", na=False)].copy()
    df_valido[target_col] = df_valido[target_col].astype(float)
    valid_models = len(df_valido)

    if valid_models == 0:
        print("[!] No viable models to analyze.")
        return

    # 3. SORT BY Kd AND EXTRACT TOP X%
    df_sorted = df_valido.sort_values(by=target_col, ascending=True)
    num_top = max(1, int((top_percent / 100.0) * valid_models))
    df_top = df_sorted.head(num_top)
    
    # Save new CSV
    os.makedirs(os.path.dirname(os.path.abspath(csv_out)), exist_ok=True)
    df_top.to_csv(csv_out, sep=";", index=False)
    
    kd_cutoff = df_top[target_col].max()

    print(f"\n--- BINDING AFFINITY ANALYSIS (Kd) ---")
    print(f"Total models processed: {total_models}")
    print(f"Valid models analyzed:  {valid_models}")
    print(f"Best absolute Kd (Min): {df_sorted[target_col].min():.2e} M")
    print(f"Worst absolute Kd (Max):{df_sorted[target_col].max():.2e} M")
    print(f"Top {top_percent}% cutoff:      {kd_cutoff:.2e} M ({num_top} models saved)")
    print(f"Top file saved at:      {csv_out}\n")

    # 4. GENERATE LOGARITHMIC HISTOGRAM
    plt.figure(figsize=(12, 7))
    
    min_kd, max_kd = df_valido[target_col].min(), df_valido[target_col].max()
    if min_kd <= 0: min_kd = 1e-15 # Prevent log(0)
    bins = np.logspace(np.log10(min_kd), np.log10(max_kd), 40)
    
    plt.hist(df_valido[target_col], bins=bins, color='mediumseagreen', edgecolor='black', alpha=0.75)
    plt.xscale('log')
    
    # 5. DRAW VERTICAL LINES
    plt.axvline(kd_cutoff, color='red', linestyle='-', linewidth=2.5, 
                label=f'Top {top_percent}% Cutoff ($\leq$ {kd_cutoff:.1e} M)')
    
    plt.axvspan(min_kd, kd_cutoff, color='red', alpha=0.1)

    if thresholds:
        colors = ['blue', 'orange', 'purple', 'brown']
        for i, ref in enumerate(thresholds):
            c = colors[i % len(colors)]
            plt.axvline(ref, color=c, linestyle='--', linewidth=2, 
                        label=f'Ref: {ref:.1e} M')

    # Plot customization
    plt.title(title, fontsize=15, fontweight='bold')
    plt.xlabel('Dissociation Constant ($K_d$) in Molar [Log Scale]', fontsize=12)
    plt.ylabel('Number of Models (Frequency)', fontsize=12)
    plt.grid(axis='y', alpha=0.5, linestyle='--')
    plt.legend(loc='upper right', fontsize=10)

    # Save image
    os.makedirs(os.path.dirname(os.path.abspath(img_out)), exist_ok=True)
    plt.savefig(img_out, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[OK] Histogram saved at: {img_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Kd Top analyzer and extractor for PRA-Pred.")
    parser.add_argument("-i", "--input", required=True, help="Original CSV file")
    parser.add_argument("-o", "--img_out", required=True, help="Output path for PNG histogram")
    parser.add_argument("-c", "--csv_out", required=True, help="Output path for new Top % CSV")
    parser.add_argument("-p", "--percent", type=float, default=5.0, help="Top percentage to extract (e.g., 5 for 5%)")
    parser.add_argument("-l", "--lines", nargs='+', type=float, help="Kd values for vertical reference lines (e.g., 1e-9 1e-6)")
    parser.add_argument("-t", "--title", default="Dissociation Constant ($K_d$) Distribution", help="Plot title")
    
    args = parser.parse_args()
    analyze_kd(args.input, args.img_out, args.csv_out, args.percent, args.lines, args.title)
