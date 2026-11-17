import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def analyze_range(input_csv, img_out, csv_out, txt_out, min_val, max_val, title):
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

    # If no limits are provided, use the real min/max of the dataset
    real_min = df_valid[target_col].min()
    real_max = df_valid[target_col].max()
    
    bound_min = min_val if min_val is not None else real_min
    bound_max = max_val if max_val is not None else real_max

    # 3. Filter by the specified range
    df_range = df_valid[(df_valid[target_col] >= bound_min) & (df_valid[target_col] <= bound_max)]
    models_in_range = len(df_range)

    # Save CSV with filtered models
    os.makedirs(os.path.dirname(os.path.abspath(csv_out)), exist_ok=True)
    df_range.to_csv(csv_out, sep=";", index=False)

    # 4. Prepare and save statistics in .txt
    info_text = (
        f"--- Kd RANGE ANALYSIS ---\n"
        f"Total source models: {total_models}\n"
        f"Selected range:      [{bound_min:.2e} M, {bound_max:.2e} M]\n"
        f"Models in range:     {models_in_range} ({(models_in_range/valid_models)*100:.1f}% of viable)\n"
    )
    print(f"\n{info_text}")
    
    os.makedirs(os.path.dirname(os.path.abspath(txt_out)), exist_ok=True)
    with open(txt_out, "w") as f:
        f.write(info_text)

    if models_in_range == 0:
        print("[!] No models in this range. No plot will be generated.")
        return

    # 5. Generate Logarithmic Histogram
    plt.figure(figsize=(12, 7))
    if real_min <= 0: real_min = 1e-15
    bins = np.logspace(np.log10(real_min), np.log10(real_max), 40)
    
    plt.hist(df_valid[target_col], bins=bins, color='lightgray', edgecolor='black', alpha=0.7)
    plt.xscale('log')
    
    # Demarcate and shade the range region
    plt.axvline(bound_min, color='blue', linestyle='--', linewidth=2)
    plt.axvline(bound_max, color='blue', linestyle='--', linewidth=2)
    plt.axvspan(bound_min, bound_max, color='blue', alpha=0.2, label=f'Range: {bound_min:.1e} to {bound_max:.1e}')

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
    parser = argparse.ArgumentParser(description="Extracts and plots models within a Kd range.")
    parser.add_argument("-i", "--input", required=True, help="Original CSV file")
    parser.add_argument("-o", "--img_out", required=True, help="Output path for PNG image")
    parser.add_argument("-c", "--csv_out", required=True, help="Output path for the new filtered CSV")
    parser.add_argument("-txt", "--txt_out", required=True, help="Output path for the statistics TXT file")
    parser.add_argument("-min", "--min_val", type=float, help="Lower limit of the Kd range (e.g., 1e-9)")
    parser.add_argument("-max", "--max_val", type=float, help="Upper limit of the Kd range (e.g., 1e-6)")
    parser.add_argument("-t", "--title", default="$K_d$ Range Analysis", help="Plot title")
    
    args = parser.parse_args()
    analyze_range(args.input, args.img_out, args.csv_out, args.txt_out, args.min_val, args.max_val, args.title)
