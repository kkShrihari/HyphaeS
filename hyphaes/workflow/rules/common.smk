# =============================================================================
# common.smk — Shared variables, paths, and helper functions
# Used by all 5 modules
# =============================================================================

from pathlib import Path

# ── General ───────────────────────────────────────────────────────────────────
OUTDIR  = config["general"]["outdir"].rstrip("/")
THREADS = config["general"]["threads"]
DB_DIR  = config["general"]["db_dir"]

# ── Per-module output dirs ────────────────────────────────────────────────────
MOD1_OUT = f"{OUTDIR}/mod1_QC"
MOD2_OUT = f"{OUTDIR}/mod2_assembly"
MOD3_OUT = f"{OUTDIR}/mod3_taxonomy"
MOD4_OUT = f"{OUTDIR}/mod4_binning"
MOD5_OUT = f"{OUTDIR}/mod5_annotation"

# ── Database paths ────────────────────────────────────────────────────────────
HOST_IDX      = config["mod1_QC"]["host_idx"]
EUKDETECT_CFG = config["mod1_QC"]["eukdetect_cfg"]
PHIX_REF      = config["mod1_QC"].get("phix_ref", "phix")

KRAKEN2_DB    = config["mod3_taxonomy"]["kraken2_db"]
KAIJU_DB      = config["mod3_taxonomy"]["kaiju_db"]
KAIJU_NODES   = config["mod3_taxonomy"]["kaiju_nodes"]

EUKCC_DB      = config["mod4_binning"]["eukcc_db"]
BUSCO_LINEAGE = config["mod4_binning"]["busco_lineage"]
AUGUSTUS_CFG  = config["mod4_binning"]["augustus_config"]

DBCAN_DB      = config["mod5_annotation"]["dbcan_db"]
ANTISMASH_DB  = config["mod5_annotation"]["antismash_db"]

# ── Sample helpers ────────────────────────────────────────────────────────────
def get_r1(wildcards):
    return samples_df.loc[wildcards.sample, "r1"]

def get_r2(wildcards):
    return samples_df.loc[wildcards.sample, "r2"]

def get_clean_r1(wildcards):
    return f"{MOD1_OUT}/04_nohost/{wildcards.sample}_R1.fq.gz"

def get_clean_r2(wildcards):
    return f"{MOD1_OUT}/04_nohost/{wildcards.sample}_R2.fq.gz"

def get_all_clean_r1(wildcards=None):
    return expand(f"{MOD1_OUT}/04_nohost/{{sample}}_R1.fq.gz", sample=SAMPLES)

def get_all_clean_r2(wildcards=None):
    return expand(f"{MOD1_OUT}/04_nohost/{{sample}}_R2.fq.gz", sample=SAMPLES)

# ── Inter-module guards ───────────────────────────────────────────────────────
def check_mod1_done():
    s = f"{MOD1_OUT}/99_multiqc_final/multiqc_report.html"
    if not Path(s).exists():
        raise ValueError(f"Module 1 not complete. Run 'HyphaeS run mod1' first.\nExpected: {s}")

def check_mod2_done():
    s = f"{MOD2_OUT}/15_coverm/coverage_table.tsv"
    if not Path(s).exists():
        raise ValueError(f"Module 2 not complete. Run 'HyphaeS run mod2' first.\nExpected: {s}")

def check_mod3_done():
    s = f"{MOD3_OUT}/16d_diversity/alpha_diversity.tsv"
    if not Path(s).exists():
        raise ValueError(f"Module 3 not complete. Run 'HyphaeS run mod3' first.\nExpected: {s}")

def check_mod4_done():
    s = f"{MOD4_OUT}/19_mag_annotation/annotation_done.flag"
    if not Path(s).exists():
        raise ValueError(f"Module 4 not complete. Run 'HyphaeS run mod4' first.\nExpected: {s}")

# Java/WSL fix
import subprocess as _sp
JAVA_BIN = _sp.run(["bash", "-c", "which java"], capture_output=True, text=True, env={**__import__("os").environ, "PATH": "/usr/bin:/usr/local/bin:" + __import__("os").environ.get("PATH","")}).stdout.strip()
JAVA_LIB_FIX = ("export PATH=" + __import__("os").path.dirname(JAVA_BIN) + ":$PATH") if JAVA_BIN else ""

# Plot script path
import pathlib
PLOT_SCRIPT = pathlib.Path(workflow.basedir) / "scripts" / "mod1_QC_plots.py"
PLOT_SCRIPT_MOD2 = pathlib.Path(workflow.basedir) / "scripts" / "mod2_assembly_plots.py"
