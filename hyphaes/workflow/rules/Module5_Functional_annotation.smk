# =============================================================================
# Module 5 — Functional Analysis
# [20a] FUNGuild  — guild assignment (saprotroph/mycorrhizal/pathogen)
# [20b] CAZyme profile prediction — novel fungi roles from dbCAN output
# [21]  Final report — MAGs, taxonomy, CAZymes, BGCs, guilds
# =============================================================================

onstart:
    check_mod4_done()

# ── Target ────────────────────────────────────────────────────────────────────
rule mod5_all:
    input:
        f"{MOD5_OUT}/20a_funguild/funguild_done.flag",
        f"{MOD5_OUT}/20b_cazyme/cazyme_profiles.tsv",
        f"{MOD5_OUT}/21_final_report/HyphaeS_final_report.html",

# ─────────────────────────────────────────────────────────────────────────────
# [20a] FUNGuild — ecological guild assignment
# Matches taxonomy to FUNGuild database
# Assigns: saprotroph / mycorrhizal / pathogen / endophyte etc.
# Input: Bracken species table (read-level) + MAG taxonomy (genome-level)
# ─────────────────────────────────────────────────────────────────────────────
rule funguild_reads:
    input:
        bracken = expand(f"{MOD3_OUT}/16c_bracken/{{sample}}.bracken", sample=SAMPLES),
    output:
        tsv  = f"{MOD5_OUT}/20a_funguild/read_level_guilds.tsv",
    params:
        outdir = f"{MOD5_OUT}/20a_funguild",
        db     = config["mod5_annotation"].get("funguild_db", "FUNGuild.db"),
    log: f"{MOD5_OUT}/logs/20a_funguild/read_level.log"
    conda: "../envs/mod5_annotation.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        # Merge all Bracken outputs into one OTU table
        python3 -c "
import pandas as pd, glob
dfs = []
for f in glob.glob('{MOD3_OUT}/16c_bracken/*.bracken'):
    sample = f.split('/')[-1].replace('.bracken','')
    df = pd.read_csv(f, sep='\t')
    df = df[['name','new_est_reads']].rename(columns={{'new_est_reads': sample}})
    dfs.append(df.set_index('name'))
merged = pd.concat(dfs, axis=1).fillna(0)
merged.index.name = 'OTU_ID'
merged.to_csv('{params.outdir}/otu_table.tsv', sep='\t')
" &> {log}

        # Run FUNGuild
        python3 -c "
import subprocess, os
subprocess.run([
    'python3', '-m', 'funguild',
    '--in', '{params.outdir}/otu_table.tsv',
    '--db', '{params.db}',
    '--out', '{output.tsv}',
    '--classifier', 'taxonomy'
], check=True)
" &>> {log}
        """

rule funguild_mags:
    input:
        mag_tax = f"{MOD4_OUT}/18_mag_taxonomy/mag_taxonomy.tsv",
    output:
        tsv = f"{MOD5_OUT}/20a_funguild/mag_level_guilds.tsv",
    params:
        outdir = f"{MOD5_OUT}/20a_funguild",
        db     = config["mod5_annotation"].get("funguild_db", "FUNGuild.db"),
    log: f"{MOD5_OUT}/logs/20a_funguild/mag_level.log"
    conda: "../envs/mod5_annotation.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        python3 -c "
import subprocess
subprocess.run([
    'python3', '-m', 'funguild',
    '--in', '{input.mag_tax}',
    '--db', '{params.db}',
    '--out', '{output.tsv}',
    '--classifier', 'taxonomy'
], check=True)
" &> {log}
        """

rule funguild_done:
    input:
        reads = f"{MOD5_OUT}/20a_funguild/read_level_guilds.tsv",
        mags  = f"{MOD5_OUT}/20a_funguild/mag_level_guilds.tsv",
    output:
        f"{MOD5_OUT}/20a_funguild/funguild_done.flag"
    shell:
        "touch {output}"

# ─────────────────────────────────────────────────────────────────────────────
# [20b] CAZyme PROFILE PREDICTION
# Combines dbCAN output (from mod4 MAG annotation) with taxonomy
# Assigns functional roles to novel fungi based on CAZyme families
# CAZyme families: GH (glycoside hydrolase), CE, PL, AA, CBM
# ─────────────────────────────────────────────────────────────────────────────
rule cazyme_profiles:
    input:
        dbcan_done = f"{MOD4_OUT}/19_mag_annotation/dbcan_done.flag",
        guilds     = f"{MOD5_OUT}/20a_funguild/mag_level_guilds.tsv",
        mag_tax    = f"{MOD4_OUT}/18_mag_taxonomy/mag_taxonomy.tsv",
    output:
        profiles = f"{MOD5_OUT}/20b_cazyme/cazyme_profiles.tsv",
        heatmap  = f"{MOD5_OUT}/20b_cazyme/cazyme_heatmap.pdf",
    params:
        dbcan_dir = f"{MOD4_OUT}/19_mag_annotation/dbcan",
        outdir    = f"{MOD5_OUT}/20b_cazyme",
    log: f"{MOD5_OUT}/logs/20b_cazyme/cazyme_profiles.log"
    conda: "../envs/mod5_annotation.yaml"
    script:
        "../scripts/cazyme_profiles.py"

# ─────────────────────────────────────────────────────────────────────────────
# [21] FINAL REPORT
# Integrates all results: MAGs, taxonomy, CAZymes, BGCs, guilds
# Generates HTML report
# ─────────────────────────────────────────────────────────────────────────────
rule final_report:
    input:
        # Module 1
        retention   = expand(f"{MOD1_OUT}/09_retention/{{sample}}_retention.tsv",
                             sample=SAMPLES),
        multiqc     = f"{MOD1_OUT}/99_multiqc_final/multiqc_report.html",
        # Module 2
        quast       = f"{MOD2_OUT}/13_quast/report.tsv",
        # Module 3
        diversity   = f"{MOD3_OUT}/16d_diversity/alpha_diversity.tsv",
        bracken_all = expand(f"{MOD3_OUT}/16c_bracken/{{sample}}.bracken",
                             sample=SAMPLES),
        # Module 4
        mag_quality = f"{MOD4_OUT}/17e_quality/eukcc_done.flag",
        mag_tax     = f"{MOD4_OUT}/18_mag_taxonomy/mag_taxonomy.tsv",
        antismash   = f"{MOD4_OUT}/19_mag_annotation/antismash_done.flag",
        # Module 5
        guilds_reads= f"{MOD5_OUT}/20a_funguild/read_level_guilds.tsv",
        guilds_mags = f"{MOD5_OUT}/20a_funguild/mag_level_guilds.tsv",
        cazymes     = f"{MOD5_OUT}/20b_cazyme/cazyme_profiles.tsv",
        cazyme_map  = f"{MOD5_OUT}/20b_cazyme/cazyme_heatmap.pdf",
    output:
        html  = f"{MOD5_OUT}/21_final_report/HyphaeS_final_report.html",
        tsv   = f"{MOD5_OUT}/21_final_report/HyphaeS_summary_table.tsv",
    params:
        outdir  = f"{MOD5_OUT}/21_final_report",
        samples = SAMPLES,
    log: f"{MOD5_OUT}/logs/21_report/final_report.log"
    conda: "../envs/mod5_annotation.yaml"
    script:
        "../scripts/final_report.py"
