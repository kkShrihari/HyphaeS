#!/usr/bin/env python3
"""
HyphaeS — mod1 QC Plot Generator
Generates publication-ready plots from mod1 QC outputs.

Usage:
    python mod1_QC_plots.py --results results/mod1_QC --outdir results/mod1_QC/plots

Plots generated:
    1. Read retention waterfall
    2. Nonpareil coverage curve
    3. SeqKit stats summary (raw vs clean)
    4. GC content comparison
    5. Combined QC summary dashboard
"""

import os
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from pathlib import Path
import glob

# ── Style ─────────────────────────────────────────────────────────────────────
COLORS = {
    "raw":      "#4C72B0",
    "trimmed":  "#55A868",
    "nophix":   "#C44E52",
    "lowcomp":  "#8172B2",
    "clean":    "#CCB974",
    "accent":   "#64B5CD",
    "bg":       "#F8F9FA",
    "grid":     "#E0E0E0",
    "text":     "#2D3436",
}

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        11,
    "axes.titlesize":   13,
    "axes.titleweight": "bold",
    "axes.labelsize":   11,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.facecolor":   COLORS["bg"],
    "figure.facecolor": "white",
    "grid.color":       COLORS["grid"],
    "grid.linewidth":   0.8,
    "xtick.labelsize":  10,
    "ytick.labelsize":  10,
})

# ── Helpers ───────────────────────────────────────────────────────────────────
def find_samples(results_dir):
    samples = []
    for f in glob.glob(f"{results_dir}/09_retention/*_retention.tsv"):
        s = Path(f).stem.replace("_retention", "")
        samples.append(s)
    return sorted(samples)


def load_retention(results_dir, samples):
    dfs = []
    for s in samples:
        f = f"{results_dir}/09_retention/{s}_retention.tsv"
        if os.path.exists(f):
            dfs.append(pd.read_csv(f, sep="\t"))
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def load_seqkit(results_dir, samples):
    dfs = []
    for s in samples:
        f = f"{results_dir}/08_seqkit/{s}_seqkit.tsv"
        if os.path.exists(f):
            df = pd.read_csv(f, sep="\t")
            df["sample"] = s
            dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def load_npo(results_dir, sample):
    f = f"{results_dir}/07_nonpareil/{sample}.npo"
    if not os.path.exists(f):
        return None, {}
    meta = {}
    rows = []
    with open(f) as fh:
        for line in fh:
            if line.startswith("# @"):
                k, v = line.strip().lstrip("# @").split(": ", 1)
                meta[k] = v
            elif not line.startswith("#"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        rows.append([float(p) for p in parts])
                    except ValueError:
                        pass
    df = pd.DataFrame(rows, columns=["seqs", "redundancy", "stderr", "q2.5", "q97.5"][:len(rows[0])] if rows else [])
    return df, meta


def load_phix_stats(results_dir, sample):
    f = f"{results_dir}/03_nophix/{sample}_phix_stats.txt"
    if not os.path.exists(f):
        return 0.0
    with open(f) as fh:
        for line in fh:
            m = re.search(r"Total Removed.*?([\d.]+)%", line)
            if m:
                return float(m.group(1))
    return 0.0


def load_lowcomp_stats(results_dir, sample):
    f = f"{results_dir}/03b_lowcomp/{sample}_lowcomp_stats.txt"
    if not os.path.exists(f):
        return 0.0
    with open(f) as fh:
        for line in fh:
            m = re.search(r"Total Removed.*?([\d.]+)%", line)
            if m:
                return float(m.group(1))
    return 0.0


# ── Plot 1: Read Retention Waterfall ─────────────────────────────────────────
def plot_retention_waterfall(retention_df, outdir):
    if retention_df.empty:
        print("  [SKIP] No retention data found.")
        return

    fig, axes = plt.subplots(1, len(retention_df), figsize=(max(8, 5 * len(retention_df)), 6),
                              squeeze=False)

    for idx, (_, row) in enumerate(retention_df.iterrows()):
        ax = axes[0][idx]
        sample = row["sample"]

        raw      = row["raw_reads"]
        trimmed  = row["trimmed_reads"] if pd.notna(row["trimmed_reads"]) else raw
        clean    = row["clean_reads"]

        # Estimate intermediate steps
        phix_lost    = trimmed * (row["phix_pct"] / 100)
        lowcomp_lost = trimmed * (row["lowcomp_pct"] / 100)
        host_lost    = trimmed * (row["host_pct"] / 100)

        steps  = ["Raw", "Trimmed", "No PhiX", "No LowComp", "Clean"]
        values = [raw, trimmed, trimmed - phix_lost,
                  trimmed - phix_lost - lowcomp_lost, clean]
        colors = [COLORS["raw"], COLORS["trimmed"], COLORS["nophix"],
                  COLORS["lowcomp"], COLORS["clean"]]

        bars = ax.bar(steps, values, color=colors, edgecolor="white",
                      linewidth=1.5, zorder=3, width=0.6)

        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + raw * 0.01,
                    f"{int(val):,}", ha="center", va="bottom",
                    fontsize=9, color=COLORS["text"], fontweight="bold")

        pct = clean / raw * 100
        ax.set_title(f"{sample}\n{pct:.1f}% retained", fontsize=12, pad=10)
        ax.set_ylabel("Read Count" if idx == 0 else "")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x/1e6:.1f}M" if x >= 1e6 else f"{int(x):,}"))
        ax.set_ylim(0, raw * 1.15)
        ax.grid(axis="y", zorder=0)
        ax.tick_params(axis="x", rotation=30)

    fig.suptitle("Read Retention Through QC Pipeline", fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    out = f"{outdir}/01_read_retention_waterfall.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ── Plot 2: Nonpareil Coverage Curve ─────────────────────────────────────────
def plot_nonpareil(results_dir, samples, outdir):
    fig, ax = plt.subplots(figsize=(9, 6))

    palette = [COLORS["raw"], COLORS["trimmed"], COLORS["nophix"],
               COLORS["lowcomp"], COLORS["clean"], COLORS["accent"]]

    plotted = False
    for i, sample in enumerate(samples):
        npo_df, meta = load_npo(results_dir, sample)
        if npo_df is None or npo_df.empty:
            continue

        color = palette[i % len(palette)]
        R = float(meta.get("R", npo_df["seqs"].max()))

        ax.plot(npo_df["seqs"], npo_df["redundancy"],
                color=color, linewidth=2.5, label=sample, zorder=3)

        if "q2.5" in npo_df.columns and "q97.5" in npo_df.columns:
            ax.fill_between(npo_df["seqs"],
                            npo_df["q2.5"], npo_df["q97.5"],
                            alpha=0.15, color=color)

        # Mark current sequencing depth
        ax.axvline(x=R, color=color, linestyle="--", linewidth=1.2, alpha=0.7)
        max_cov = npo_df["redundancy"].max()
        ax.annotate(f"{max_cov:.1%}",
                    xy=(R, max_cov),
                    xytext=(R * 1.05, max_cov),
                    fontsize=9, color=color, fontweight="bold")
        plotted = True

    if not plotted:
        print("  [SKIP] No nonpareil data found.")
        plt.close()
        return

    # Reference lines
    ax.axhline(y=0.95, color="gray", linestyle=":", linewidth=1.5,
               label="95% coverage threshold")
    ax.axhline(y=1.0, color="#2D3436", linestyle="-", linewidth=1, alpha=0.3)

    ax.set_xlabel("Sequenced Reads", fontsize=12)
    ax.set_ylabel("Estimated Coverage (Redundancy)", fontsize=12)
    ax.set_title("Nonpareil Sequencing Coverage Curve", fontsize=14, pad=12)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0%}"))
    ax.legend(frameon=True, framealpha=0.9, fontsize=10)
    ax.grid(True, zorder=0)

    plt.tight_layout()
    out = f"{outdir}/02_nonpareil_coverage.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ── Plot 3: SeqKit Stats — Raw vs Clean ──────────────────────────────────────
def plot_seqkit_stats(seqkit_df, outdir):
    if seqkit_df.empty:
        print("  [SKIP] No seqkit data found.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    metrics = [
        ("num_seqs",    "Read Count",       lambda x: f"{x/1e6:.1f}M"),
        ("avg_len",     "Average Length",   lambda x: f"{x:.0f}bp"),
        ("gc_content",  "GC Content (%)",   lambda x: f"{x:.1f}%"),
    ]

    for ax, (col, label, fmt) in zip(axes, metrics):
        if col not in seqkit_df.columns:
            continue

        raw_data   = seqkit_df[seqkit_df["file"].str.contains("R1.fastq|R1.fq", na=False) &
                               ~seqkit_df["file"].str.contains("trimmed|nohost|04_nohost", na=False)]
        clean_data = seqkit_df[seqkit_df["file"].str.contains("04_nohost", na=False) &
                               seqkit_df["file"].str.contains("R1", na=False)]

        samples = seqkit_df["sample"].unique()
        x = np.arange(len(samples))
        w = 0.35

        raw_vals   = [raw_data[raw_data["sample"] == s][col].values[0]
                      if len(raw_data[raw_data["sample"] == s]) > 0 else 0
                      for s in samples]
        clean_vals = [clean_data[clean_data["sample"] == s][col].values[0]
                      if len(clean_data[clean_data["sample"] == s]) > 0 else 0
                      for s in samples]

        ax.bar(x - w/2, raw_vals,   w, label="Raw",   color=COLORS["raw"],   alpha=0.85, edgecolor="white")
        ax.bar(x + w/2, clean_vals, w, label="Clean", color=COLORS["clean"], alpha=0.85, edgecolor="white")

        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=30, ha="right")
        ax.set_title(label)
        ax.grid(axis="y", zorder=0)
        ax.legend(fontsize=9)

        for i, (rv, cv) in enumerate(zip(raw_vals, clean_vals)):
            ax.text(i - w/2, rv * 1.02, fmt(rv), ha="center", va="bottom", fontsize=8)
            ax.text(i + w/2, cv * 1.02, fmt(cv), ha="center", va="bottom", fontsize=8)

    fig.suptitle("SeqKit Stats: Raw vs Clean Reads", fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = f"{outdir}/03_seqkit_raw_vs_clean.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ── Plot 4: Contamination Summary ────────────────────────────────────────────
def plot_contamination(retention_df, outdir):
    if retention_df.empty:
        print("  [SKIP] No retention data for contamination plot.")
        return

    fig, ax = plt.subplots(figsize=(max(7, 2.5 * len(retention_df)), 5))

    contam_types = ["phix_pct", "lowcomp_pct", "host_pct"]
    contam_labels = ["PhiX", "Low Complexity", "Host (Human)"]
    contam_colors = [COLORS["nophix"], COLORS["lowcomp"], COLORS["accent"]]

    x = np.arange(len(retention_df))
    w = 0.25

    for i, (col, label, color) in enumerate(zip(contam_types, contam_labels, contam_colors)):
        vals = retention_df[col].fillna(0).values
        bars = ax.bar(x + i * w - w, vals, w, label=label,
                      color=color, alpha=0.85, edgecolor="white", zorder=3)
        for bar, val in zip(bars, vals):
            if val > 0.01:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.02,
                        f"{val:.2f}%", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(retention_df["sample"], rotation=30, ha="right")
    ax.set_ylabel("Contamination (%)")
    ax.set_title("Contamination Sources per Sample", fontsize=14, pad=12)
    ax.legend(frameon=True, fontsize=10)
    ax.grid(axis="y", zorder=0)
    ax.set_ylim(0, max(retention_df[contam_types].max().max() * 1.4, 1))

    plt.tight_layout()
    out = f"{outdir}/04_contamination_summary.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ── Plot 5: QC Summary Dashboard ─────────────────────────────────────────────
def plot_dashboard(retention_df, outdir):
    if retention_df.empty:
        print("  [SKIP] No data for dashboard.")
        return

    fig = plt.figure(figsize=(14, 8))
    fig.patch.set_facecolor("white")
    gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.4)

    # Panel 1: Retention %
    ax1 = fig.add_subplot(gs[0, 0])
    samples = retention_df["sample"].tolist()
    pcts = retention_df["pct_retained"].tolist()
    colors = [COLORS["clean"] if p >= 90 else COLORS["nophix"] for p in pcts]
    bars = ax1.barh(samples, pcts, color=colors, edgecolor="white", height=0.5)
    ax1.axvline(90, color="gray", linestyle="--", linewidth=1.2)
    ax1.set_xlim(0, 105)
    ax1.set_xlabel("% Reads Retained")
    ax1.set_title("Read Retention")
    for bar, pct in zip(bars, pcts):
        ax1.text(pct + 0.5, bar.get_y() + bar.get_height() / 2,
                 f"{pct:.1f}%", va="center", fontsize=9)

    # Panel 2: Raw vs Clean read counts
    ax2 = fig.add_subplot(gs[0, 1])
    x = np.arange(len(samples))
    ax2.bar(x - 0.2, retention_df["raw_reads"] / 1e6,   0.35, label="Raw",   color=COLORS["raw"],   alpha=0.85)
    ax2.bar(x + 0.2, retention_df["clean_reads"] / 1e6, 0.35, label="Clean", color=COLORS["clean"], alpha=0.85)
    ax2.set_xticks(x)
    ax2.set_xticklabels(samples, rotation=30, ha="right", fontsize=9)
    ax2.set_ylabel("Reads (M)")
    ax2.set_title("Raw vs Clean Reads")
    ax2.legend(fontsize=9)
    ax2.grid(axis="y")

    # Panel 3: Host contamination
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.bar(samples, retention_df["host_pct"], color=COLORS["accent"], alpha=0.85, edgecolor="white")
    ax3.set_ylabel("Host Reads (%)")
    ax3.set_title("Host Contamination")
    ax3.tick_params(axis="x", rotation=30)
    ax3.grid(axis="y")

    # Panel 4: PhiX + LowComp
    ax4 = fig.add_subplot(gs[1, 0])
    x = np.arange(len(samples))
    ax4.bar(x - 0.2, retention_df["phix_pct"],    0.35, label="PhiX",        color=COLORS["nophix"], alpha=0.85)
    ax4.bar(x + 0.2, retention_df["lowcomp_pct"], 0.35, label="Low Complexity", color=COLORS["lowcomp"], alpha=0.85)
    ax4.set_xticks(x)
    ax4.set_xticklabels(samples, rotation=30, ha="right", fontsize=9)
    ax4.set_ylabel("Reads Removed (%)")
    ax4.set_title("PhiX & Low Complexity")
    ax4.legend(fontsize=9)
    ax4.grid(axis="y")

    # Panel 5: Nonpareil redundancy
    ax5 = fig.add_subplot(gs[1, 1])
    nr = retention_df["nonpareil_redundancy"].replace("NA", np.nan).astype(float)
    colors5 = [COLORS["clean"] if v >= 0.95 else COLORS["nophix"] if v < 0.5 else COLORS["trimmed"]
               for v in nr.fillna(0)]
    ax5.bar(samples, nr, color=colors5, alpha=0.85, edgecolor="white")
    ax5.axhline(0.95, color="gray", linestyle="--", linewidth=1.2, label="95% threshold")
    ax5.set_ylabel("Nonpareil Redundancy")
    ax5.set_title("Sequencing Coverage")
    ax5.set_ylim(0, 1.05)
    ax5.tick_params(axis="x", rotation=30)
    ax5.legend(fontsize=9)
    ax5.grid(axis="y")

    # Panel 6: Summary table
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis("off")
    table_data = []
    for _, row in retention_df.iterrows():
        table_data.append([
            row["sample"][:12],
            f"{int(row['raw_reads']):,}",
            f"{row['pct_retained']:.1f}%",
            f"{row['host_pct']:.2f}%",
        ])
    table = ax6.table(
        cellText=table_data,
        colLabels=["Sample", "Raw Reads", "Retained", "Host%"],
        cellLoc="center", loc="center",
        bbox=[0, 0, 1, 1]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor(COLORS["raw"])
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#F0F4F8")
    ax6.set_title("Summary Table", fontsize=11, fontweight="bold", pad=8)

    fig.suptitle("HyphaeS — mod1 QC Summary Dashboard",
                 fontsize=16, fontweight="bold", y=1.01)

    out = f"{outdir}/05_QC_dashboard.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="HyphaeS mod1 QC Plot Generator")
    parser.add_argument("--results", default="results/mod1_QC",
                        help="Path to mod1_QC results directory")
    parser.add_argument("--outdir", default="results/mod1_QC/plots",
                        help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"\nHyphaeS — Generating mod1 QC plots")
    print(f"  Results : {args.results}")
    print(f"  Output  : {args.outdir}\n")

    samples = find_samples(args.results)
    if not samples:
        print("ERROR: No samples found. Check results directory.")
        return

    print(f"  Samples found: {', '.join(samples)}\n")

    retention_df = load_retention(args.results, samples)
    seqkit_df    = load_seqkit(args.results, samples)

    print("Generating plots...")
    plot_retention_waterfall(retention_df, args.outdir)
    plot_nonpareil(args.results, samples, args.outdir)
    plot_seqkit_stats(seqkit_df, args.outdir)
    plot_contamination(retention_df, args.outdir)
    plot_dashboard(retention_df, args.outdir)

    print(f"\nAll plots saved to: {args.outdir}/")
    print("Files:")
    for f in sorted(os.listdir(args.outdir)):
        print(f"  {f}")


if __name__ == "__main__":
    main()
