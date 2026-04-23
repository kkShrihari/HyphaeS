#!/usr/bin/env python3
"""
HyphaeS — mod2 Assembly Plot Generator
Generates publication-ready plots from mod2 assembly outputs.

Usage:
    python mod2_assembly_plots.py --results results/mod2_assembly --outdir results/mod2_assembly/plots

Plots generated:
    1. Contig length distribution histogram
    2. Assembly stats bar chart (N50, total length, contig count)
    3. Per-sample read mapping rate
    4. Coverage distribution per contig
    5. Assembly summary dashboard
"""

import os
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import glob

COLORS = {
    "primary":  "#2E86AB",
    "secondary":"#A23B72",
    "accent":   "#F18F01",
    "green":    "#44BBA4",
    "red":      "#E94F37",
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
})

def load_contigs(fasta_path):
    lengths = []
    current_len = 0
    try:
        with open(fasta_path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if current_len > 0:
                        lengths.append(current_len)
                    current_len = 0
                else:
                    current_len += len(line)
            if current_len > 0:
                lengths.append(current_len)
    except FileNotFoundError:
        pass
    return sorted(lengths, reverse=True)

def calc_n50(lengths):
    if not lengths:
        return 0
    total = sum(lengths)
    cumsum = 0
    for l in sorted(lengths, reverse=True):
        cumsum += l
        if cumsum >= total / 2:
            return l
    return 0

def calc_n75(lengths):
    if not lengths:
        return 0
    total = sum(lengths)
    cumsum = 0
    for l in sorted(lengths, reverse=True):
        cumsum += l
        if cumsum >= total * 0.75:
            return l
    return 0

def load_coverm(results_dir):
    f = f"{results_dir}/15_coverm/coverage_table.tsv"
    if not os.path.exists(f):
        return None
    return pd.read_csv(f, sep="\t")

def load_bowtie2_logs(results_dir, samples):
    rates = {}
    for s in samples:
        log = f"{results_dir}/logs/14_mapping/{s}.log"
        if os.path.exists(log):
            with open(log) as fh:
                for line in fh:
                    m = re.search(r"([\d.]+)% overall alignment rate", line)
                    if m:
                        rates[s] = float(m.group(1))
    return rates

def find_samples(results_dir):
    bams = glob.glob(f"{results_dir}/14_mapping/*.sorted.bam")
    return [Path(b).stem.replace(".sorted", "") for b in bams]

def plot_contig_distribution(results_dir, outdir):
    fasta = f"{results_dir}/12_filtered/contigs_filtered.fasta"
    lengths = load_contigs(fasta)
    if not lengths:
        print("  [SKIP] No contigs found.")
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax1 = axes[0]
    bins = np.logspace(np.log10(min(lengths)), np.log10(max(lengths)), 40)
    ax1.hist(lengths, bins=bins, color=COLORS["primary"], edgecolor="white", linewidth=0.5, alpha=0.85, zorder=3)
    ax1.set_xscale("log")
    ax1.set_xlabel("Contig Length (bp)")
    ax1.set_ylabel("Count")
    ax1.set_title("Contig Length Distribution")
    ax1.grid(True, zorder=0)
    ax1.axvline(1000,  color=COLORS["accent"], linestyle="--", linewidth=1.5, label="1 kb")
    ax1.axvline(10000, color=COLORS["red"],    linestyle="--", linewidth=1.5, label="10 kb")
    ax1.legend(fontsize=9)
    ax2 = axes[1]
    sorted_lens = sorted(lengths, reverse=True)
    cumlen = np.cumsum(sorted_lens)
    ax2.plot(range(1, len(sorted_lens)+1), cumlen / 1e6, color=COLORS["primary"], linewidth=2.5, zorder=3)
    ax2.fill_between(range(1, len(sorted_lens)+1), cumlen / 1e6, alpha=0.15, color=COLORS["primary"])
    n50 = calc_n50(lengths)
    n50_idx = next((i for i, l in enumerate(sorted_lens) if l <= n50), len(sorted_lens))
    ax2.axvline(n50_idx, color=COLORS["accent"], linestyle="--", linewidth=1.5, label=f"N50 = {n50:,} bp")
    ax2.set_xlabel("Contig Rank")
    ax2.set_ylabel("Cumulative Length (Mb)")
    ax2.set_title("Cumulative Assembly Length")
    ax2.legend(fontsize=9)
    ax2.grid(True, zorder=0)
    fig.suptitle(f"Contig Distribution — {len(lengths):,} contigs, {sum(lengths)/1e6:.1f} Mb total", fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = f"{outdir}/01_contig_distribution.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

def plot_assembly_stats(results_dir, outdir):
    fasta = f"{results_dir}/12_filtered/contigs_filtered.fasta"
    lengths = load_contigs(fasta)
    if not lengths:
        print("  [SKIP] No assembly stats.")
        return
    n50 = calc_n50(lengths)
    n75 = calc_n75(lengths)
    total = sum(lengths)
    count = len(lengths)
    maxlen = max(lengths)
    meanlen = np.mean(lengths)
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    ax1 = axes[0]
    metrics = ["Total\nLength (Mb)", "N50\n(kb)", "Largest\nContig (kb)", "Mean\nLength (bp)"]
    values  = [total/1e6, n50/1e3, maxlen/1e3, meanlen]
    colors  = [COLORS["primary"], COLORS["accent"], COLORS["green"], COLORS["secondary"]]
    bars = ax1.bar(metrics, values, color=colors, edgecolor="white", width=0.5, zorder=3)
    for bar, val in zip(bars, values):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.1f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax1.set_title("Assembly Key Metrics")
    ax1.grid(axis="y", zorder=0)
    ax2 = axes[1]
    cat_counts = [
        sum(1 for l in lengths if l > 50000),
        sum(1 for l in lengths if 10000 < l <= 50000),
        sum(1 for l in lengths if 5000  < l <= 10000),
        sum(1 for l in lengths if 1000  <= l <= 5000),
    ]
    cat_labels = [">50kb", "10-50kb", "5-10kb", "1-5kb"]
    filtered = [(c, l) for c, l in zip(cat_counts, cat_labels) if c > 0]
    if filtered:
        counts_f, labels_f = zip(*filtered)
        ax2.pie(counts_f, labels=labels_f,
                colors=[COLORS["red"], COLORS["accent"], COLORS["green"], COLORS["primary"]][:len(counts_f)],
                autopct="%1.1f%%", startangle=90, wedgeprops={"edgecolor": "white", "linewidth": 1.5})
    ax2.set_title("Contig Size Categories")
    ax3 = axes[2]
    ax3.axis("off")
    table_data = [
        ["Total contigs",   f"{count:,}"],
        ["Total length",    f"{total/1e6:.2f} Mb"],
        ["N50",             f"{n50:,} bp"],
        ["N75",             f"{n75:,} bp"],
        ["Largest contig",  f"{maxlen:,} bp"],
        ["Mean length",     f"{meanlen:.0f} bp"],
        [">10 kb contigs",  f"{sum(1 for l in lengths if l>10000):,}"],
    ]
    table = ax3.table(cellText=table_data, colLabels=["Metric", "Value"], cellLoc="left", loc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor(COLORS["primary"])
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#F0F4F8")
    ax3.set_title("Assembly Statistics", fontsize=11, fontweight="bold", pad=8)
    fig.suptitle("Assembly Quality Assessment", fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = f"{outdir}/02_assembly_stats.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

def plot_mapping_rates(results_dir, samples, outdir):
    rates = load_bowtie2_logs(results_dir, samples)
    if not rates:
        print("  [SKIP] No mapping rate logs found.")
        return
    fig, ax = plt.subplots(figsize=(max(7, 2.5*len(rates)), 5))
    colors = [COLORS["green"] if v >= 80 else COLORS["accent"] if v >= 50 else COLORS["red"] for v in rates.values()]
    bars = ax.bar(list(rates.keys()), list(rates.values()), color=colors, edgecolor="white", width=0.5, zorder=3)
    ax.axhline(80, color="gray", linestyle="--", linewidth=1.2, label="80% threshold")
    ax.set_ylim(0, 105)
    ax.set_ylabel("Mapping Rate (%)")
    ax.set_title("Per-Sample Read Mapping to Co-Assembly", fontsize=14, pad=12)
    ax.legend(fontsize=10)
    ax.grid(axis="y", zorder=0)
    ax.tick_params(axis="x", rotation=30)
    for bar, val in zip(bars, rates.values()):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, f"{val:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
    plt.tight_layout()
    out = f"{outdir}/03_mapping_rates.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

def plot_coverage_distribution(results_dir, outdir):
    cov_df = load_coverm(results_dir)
    if cov_df is None or cov_df.empty:
        print("  [SKIP] No coverage data found.")
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    cov_cols = [c for c in cov_df.columns if c != "Contig"]
    palette = [COLORS["primary"], COLORS["accent"], COLORS["green"], COLORS["secondary"], COLORS["red"]]
    ax1 = axes[0]
    for i, col in enumerate(cov_cols):
        vals = cov_df[col].dropna()
        vals = vals[vals > 0]
        if len(vals) > 0:
            ax1.hist(np.log10(vals + 0.01), bins=50, alpha=0.6, color=palette[i % len(palette)], label=col.replace(" Mean", ""), edgecolor="none")
    ax1.set_xlabel("log10(Coverage + 0.01)")
    ax1.set_ylabel("Contig Count")
    ax1.set_title("Coverage Distribution per Sample")
    ax1.legend(fontsize=9)
    ax1.grid(True, zorder=0)
    ax2 = axes[1]
    if cov_cols:
        mean_cov = cov_df[cov_cols].mean()
        ax2.bar(range(len(mean_cov)), mean_cov.values, color=palette[:len(mean_cov)], edgecolor="white", zorder=3)
        ax2.set_xticks(range(len(mean_cov)))
        ax2.set_xticklabels([c.replace(" Mean", "") for c in mean_cov.index], rotation=30, ha="right")
        ax2.set_ylabel("Mean Coverage (x)")
        ax2.set_title("Mean Contig Coverage per Sample")
        ax2.grid(axis="y", zorder=0)
        for i, val in enumerate(mean_cov.values):
            ax2.text(i, val * 1.02, f"{val:.1f}x", ha="center", va="bottom", fontsize=9)
    fig.suptitle("Coverage Analysis", fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = f"{outdir}/04_coverage_distribution.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

def plot_dashboard(results_dir, samples, outdir):
    fasta   = f"{results_dir}/12_filtered/contigs_filtered.fasta"
    lengths = load_contigs(fasta)
    rates   = load_bowtie2_logs(results_dir, samples)
    cov_df  = load_coverm(results_dir)
    if not lengths:
        print("  [SKIP] No data for dashboard.")
        return
    fig = plt.figure(figsize=(16, 9))
    gs  = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.4)
    n50 = calc_n50(lengths)
    n75 = calc_n75(lengths)
    ax1 = fig.add_subplot(gs[0, 0])
    bins = np.logspace(np.log10(min(lengths)), np.log10(max(lengths)), 30)
    ax1.hist(lengths, bins=bins, color=COLORS["primary"], edgecolor="white", linewidth=0.5, alpha=0.85, zorder=3)
    ax1.set_xscale("log")
    ax1.set_xlabel("Length (bp)")
    ax1.set_ylabel("Count")
    ax1.set_title("Contig Length Distribution")
    ax1.grid(True, zorder=0)
    ax2 = fig.add_subplot(gs[0, 1])
    stats = ["Contigs", "Total (Mb)", "N50 (kb)", "Max (kb)"]
    vals  = [len(lengths), sum(lengths)/1e6, n50/1e3, max(lengths)/1e3]
    bars = ax2.bar(stats, vals, color=[COLORS["primary"], COLORS["accent"], COLORS["green"], COLORS["secondary"]], edgecolor="white", width=0.5, zorder=3)
    for bar, val in zip(bars, vals):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.1f}", ha="center", va="bottom", fontsize=9)
    ax2.set_title("Key Metrics")
    ax2.grid(axis="y", zorder=0)
    ax2.tick_params(axis="x", rotation=20)
    ax3 = fig.add_subplot(gs[0, 2])
    cat_counts = [sum(1 for l in lengths if l > 50000), sum(1 for l in lengths if 10000<l<=50000), sum(1 for l in lengths if 5000<l<=10000), sum(1 for l in lengths if 1000<=l<=5000)]
    cat_labels = [">50kb", "10-50kb", "5-10kb", "1-5kb"]
    filtered = [(c, lb) for c, lb in zip(cat_counts, cat_labels) if c > 0]
    if filtered:
        counts_f, labels_f = zip(*filtered)
        ax3.pie(counts_f, labels=labels_f, colors=[COLORS["red"], COLORS["accent"], COLORS["green"], COLORS["primary"]][:len(counts_f)], autopct="%1.1f%%", startangle=90, wedgeprops={"edgecolor":"white","linewidth":1.5})
    ax3.set_title("Size Categories")
    ax4 = fig.add_subplot(gs[1, 0])
    if rates:
        ax4.bar(list(rates.keys()), list(rates.values()), color=COLORS["green"], edgecolor="white", width=0.5, zorder=3)
        ax4.axhline(80, color="gray", linestyle="--", linewidth=1.2)
        ax4.set_ylim(0, 105)
        ax4.set_ylabel("Mapping Rate (%)")
        ax4.set_title("Read Mapping Rate")
        ax4.tick_params(axis="x", rotation=30)
        ax4.grid(axis="y", zorder=0)
    else:
        ax4.text(0.5, 0.5, "No mapping\ndata yet", ha="center", va="center", transform=ax4.transAxes, fontsize=12, color="gray")
        ax4.set_title("Read Mapping Rate")
    ax5 = fig.add_subplot(gs[1, 1])
    if cov_df is not None and not cov_df.empty:
        cov_cols = [c for c in cov_df.columns if c != "Contig"]
        if cov_cols:
            vals2 = cov_df[cov_cols[0]].dropna()
            vals2 = vals2[vals2 > 0]
            ax5.hist(np.log10(vals2 + 0.01), bins=40, color=COLORS["primary"], edgecolor="none", alpha=0.85)
            ax5.set_xlabel("log10(Coverage)")
            ax5.set_ylabel("Contig Count")
    ax5.set_title("Coverage Distribution")
    ax5.grid(True, zorder=0)
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis("off")
    table_data = [["Total contigs", f"{len(lengths):,}"], ["Total length", f"{sum(lengths)/1e6:.2f} Mb"], ["N50", f"{n50:,} bp"], ["N75", f"{n75:,} bp"], ["Largest", f"{max(lengths):,} bp"], [">10 kb", f"{sum(1 for l in lengths if l>10000):,}"]]
    table = ax6.table(cellText=table_data, colLabels=["Metric", "Value"], cellLoc="left", loc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor(COLORS["primary"])
            cell.set_text_props(color="white", fontweight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#F0F4F8")
    ax6.set_title("Assembly Summary", fontsize=11, fontweight="bold", pad=8)
    fig.suptitle("HyphaeS — mod2 Assembly Summary Dashboard", fontsize=16, fontweight="bold", y=1.01)
    out = f"{outdir}/05_assembly_dashboard.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

def main():
    parser = argparse.ArgumentParser(description="HyphaeS mod2 Assembly Plot Generator")
    parser.add_argument("--results", default="results/mod2_assembly")
    parser.add_argument("--outdir",  default="results/mod2_assembly/plots")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    print(f"\nHyphaeS — Generating mod2 Assembly plots")
    print(f"  Results : {args.results}")
    print(f"  Output  : {args.outdir}\n")
    samples = find_samples(args.results)
    print(f"  Samples found: {', '.join(samples) if samples else 'none yet'}\n")
    print("Generating plots...")
    plot_contig_distribution(args.results, args.outdir)
    plot_assembly_stats(args.results, args.outdir)
    plot_mapping_rates(args.results, samples, args.outdir)
    plot_coverage_distribution(args.results, args.outdir)
    plot_dashboard(args.results, samples, args.outdir)
    print(f"\nAll plots saved to: {args.outdir}/")

if __name__ == "__main__":
    main()
