import os
import subprocess
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
from matplotlib.patches import ConnectionPatch, Rectangle


def generate_single_weblogo(df_subset, seq_col, tmp_fasta, tmp_png):
    valid_seqs = 0
    with open(tmp_fasta, "w") as f_out:
        for idx, row in df_subset.iterrows():
            seq = str(row[seq_col]).strip()
            if seq != "nan" and len(seq) > 0:
                f_out.write(f">Seq_{idx+1}\n{seq}\n")
                valid_seqs += 1

    if valid_seqs == 0:
        return False

    command = [
        "weblogo",
        "-f", tmp_fasta,
        "-o", tmp_png,
        "-F", "png",
        "--resolution", "600",
        "--size", "large",
        "--units", "probability",
        "--sequence-type", "rna",
        "--fineprint", ""
    ]
    subprocess.run(command, check=True)
    return True


def plot_composite(input_csv, img_out, percentiles, seq_col, title):
    # 1. Load CSV and find Kd column
    try:
        df = pd.read_csv(input_csv, sep=";")
        if len(df.columns) == 1:
            df = pd.read_csv(input_csv, sep=",")
    except Exception as e:
        print(f"[ERROR] Could not read {input_csv}: {e}")
        return

    kd_cols = [col for col in df.columns if "kd" in col.lower()]
    if not kd_cols:
        print("[ERROR] No Kd column found in the CSV.")
        return
    target_col = kd_cols[0]

    if seq_col not in df.columns:
        print(f"[ERROR] Sequence column '{seq_col}' not found. Available: {list(df.columns)}")
        return

    # 2. Filter errors and sort by Kd
    df_valid = df[~df[target_col].astype(str).str.contains("ERROR", na=False)].copy()
    df_valid[target_col] = df_valid[target_col].astype(float)
    valid_models = len(df_valid)

    if valid_models == 0:
        print("[!] No viable models to plot.")
        return

    df_sorted = df_valid.sort_values(by=target_col, ascending=True)
    percentiles_sorted = sorted(percentiles)  # Smallest at bottom, largest at top
    n_pct = len(percentiles_sorted)

    out_dir = os.path.dirname(os.path.abspath(img_out))
    os.makedirs(out_dir, exist_ok=True)

    # 3. Generate WebLogo PNG for each percentile
    logo_files = []
    cutoffs = []
    print("\n--- GENERATING COMPOSITE FIGURE ---")
    for pct in percentiles_sorted:
        num_top = max(1, int((pct / 100.0) * valid_models))
        df_top = df_sorted.head(num_top)
        kd_cutoff = df_top[target_col].max()
        cutoffs.append(kd_cutoff)

        tmp_fasta = os.path.join(out_dir, f"_tmp_top{pct:g}.fasta")
        tmp_png = os.path.join(out_dir, f"_tmp_top{pct:g}.png")
        print(f"[*] Top {pct:g}%: Kd <= {kd_cutoff:.2e} M ({num_top} models) -> Generating WebLogo...")
        ok = generate_single_weblogo(df_top, seq_col, tmp_fasta, tmp_png)
        logo_files.append(tmp_png if ok else None)

        if os.path.exists(tmp_fasta):
            os.remove(tmp_fasta)

    # 4. Layout: balanced proportions with larger WebLogos
    fig = plt.figure(figsize=(15, max(7, 2.7 * n_pct)))
    gs = gridspec.GridSpec(
        n_pct, 2,
        width_ratios=[1.35, 1.65],
        wspace=0.16,
        hspace=0.15
    )

    # Right panel: Classic histogram with exact cutoff-aligned bar coloring
    ax_hist = fig.add_subplot(gs[:, 1])
    min_kd, max_kd = df_valid[target_col].min(), df_valid[target_col].max()
    if min_kd <= 0:
        min_kd = 1e-15

    log_min, log_max = np.log10(min_kd), np.log10(max_kd)
    bins = np.logspace(log_min, log_max, 40)

    # Draw only black edges of histogram bars (above colored fills)
    counts, bin_edges, _ = ax_hist.hist(
        df_valid[target_col], bins=bins, facecolor="none", edgecolor="black", zorder=3
    )

    # Define exact color intervals based on each Kd cutoff
    colors = ["#d62728", "#ff7f0e", "#1f77b4", "#9467bd", "#2ca02c"]
    intervals = []
    prev_cut = 0.0
    for i, kd_cut in enumerate(cutoffs):
        intervals.append((prev_cut, kd_cut, colors[i % len(colors)], 0.55))
        prev_cut = kd_cut
    intervals.append((prev_cut, np.inf, "lightgray", 0.8))

    # Fill each bar, splitting precisely where each cutoff line falls
    for left_edge, right_edge, height in zip(bin_edges[:-1], bin_edges[1:], counts):
        if height == 0:
            continue
        for a, b, color, alpha in intervals:
            sub_left = max(left_edge, a)
            sub_right = min(right_edge, b)
            if sub_right > sub_left:
                ax_hist.add_patch(
                    Rectangle(
                        (sub_left, 0), sub_right - sub_left, height,
                        facecolor=color, alpha=alpha, edgecolor="none", zorder=2
                    )
                )

    ax_hist.set_xscale("log")
    max_y = max(counts) if len(counts) > 0 else 1
    ylim_max = max_y * 1.25
    ax_hist.set_ylim(0, ylim_max)
    ax_hist.set_xlabel("Dissociation Constant ($K_d$) in Molar [Log Scale]", fontsize=12)
    ax_hist.set_ylabel("Number of Models (Frequency)", fontsize=12)
    ax_hist.grid(axis="y", alpha=0.4, linestyle="--", zorder=1)

    # Left panels: WebLogos (smallest percentile at bottom, largest at top)
    logo_axes = []
    for i in range(n_pct):
        row_idx = n_pct - 1 - i
        ax_logo = fig.add_subplot(gs[row_idx, 0])
        if logo_files[i] and os.path.exists(logo_files[i]):
            img = mpimg.imread(logo_files[i])
            ax_logo.imshow(img)
        ax_logo.axis("off")
        logo_axes.append(ax_logo)

    fig.canvas.draw()
    hist_box = ax_hist.get_position()

    # 5. L-shaped arrows from each cutoff line to its corresponding WebLogo
    for i, (pct, kd_cutoff, ax_logo) in enumerate(zip(percentiles_sorted, cutoffs, logo_axes)):
        color = colors[i % len(colors)]
        logo_box = ax_logo.get_position()
        logo_mid_y = (logo_box.y0 + logo_box.y1) / 2.0

        frac_y = (logo_mid_y - hist_box.y0) / (hist_box.y1 - hist_box.y0)
        y_data = frac_y * ylim_max

        ax_hist.plot(
            [kd_cutoff, kd_cutoff], [0, y_data],
            color=color, linestyle="-", linewidth=1.8, zorder=4
        )

        arrow = ConnectionPatch(
            xyA=(kd_cutoff, y_data), coordsA=ax_hist.transData,
            xyB=(1.02, 0.5), coordsB=ax_logo.transAxes,
            arrowstyle="->", color=color, lw=2.0, mutation_scale=16, clip_on=False
        )
        fig.add_artist(arrow)

        ax_hist.text(
            kd_cutoff,
            y_data + ylim_max * 0.02,
            f" {pct:g}%",
            color=color,
            fontsize=11,
            fontweight="bold",
            ha="right",
            va="bottom",
            zorder=5
        )

    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.98)
    plt.savefig(img_out, dpi=600, bbox_inches="tight")
    plt.close()

    for f in logo_files:
        if f and os.path.exists(f):
            os.remove(f)

    print(f"\n[OK] Composite figure saved at: {img_out}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates a composite Kd histogram + stacked WebLogos figure.")
    parser.add_argument("-i", "--input", required=True, help="Path to original CSV file")
    parser.add_argument("-o", "--img_out", required=True, help="Output path for the composite PNG")
    parser.add_argument("-p", "--percentiles", type=float, nargs="+", default=[1, 5, 10],
                        help="Percentiles to plot (e.g., -p 1 5 10)")
    parser.add_argument("-seq", "--seq_col", default="RNA_Sequence", help="Sequence column name")
    parser.add_argument("-t", "--title", default="Predicted $K_d$ Distribution & RNA Motifs", help="Figure title")

    args = parser.parse_args()
    plot_composite(args.input, args.img_out, args.percentiles, args.seq_col, args.title)
