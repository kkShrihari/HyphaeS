# =============================================================================
# Module 4 — Binning
# [17a] EukRep + Tiara    — eukaryote contig filter, fungi contigs only
# [17b] Ensemble binning  — MetaBAT2 + MaxBin2 + CONCOCT   (mode: ensemble)
#        OR HyphaeBin     — novel fungi binning tool         (mode: hyperbin)
# [17c] DAS_Tool          — picks best non-overlapping bins  (ensemble only)
# [17d] HyphaeBin         — PLACEHOLDER (pass) to be implemented
# [17e] EukCC + BUSCO + GUNC — MAG quality, completeness, chimerism
# [18]  MAG taxonomy      — Kaiju + EPA-ng + IQ-TREE
# [19]  MAG annotation    — Augustus + eggNOG + dbCAN + antiSMASH
# =============================================================================
# Config option:
#   binning_mode: "ensemble"   → MetaBAT2+MaxBin2+CONCOCT → DAS_Tool
#   binning_mode: "hyperbin"   → HyphaeBin only
# =============================================================================

onstart:
    check_mod2_done()

BINNING_MODE = config["mod4_binning"].get("binning_mode", "ensemble")
BINS_DIR     = f"{MOD4_OUT}/17c_dastool/bins" if BINNING_MODE == "ensemble" \
               else f"{MOD4_OUT}/17d_hyperbin/bins"

# ── Target ────────────────────────────────────────────────────────────────────
rule mod4_all:
    input:
        f"{MOD4_OUT}/17e_quality/eukcc_done.flag",
        f"{MOD4_OUT}/18_mag_taxonomy/mag_taxonomy.tsv",
        f"{MOD4_OUT}/19_mag_annotation/annotation_done.flag",

# ─────────────────────────────────────────────────────────────────────────────
# [17a] EukRep + Tiara — eukaryote contig filter
# EukRep: ML-based euk vs prok classification
# Tiara: deep-learning euk classifier, confirms fungi signal
# ─────────────────────────────────────────────────────────────────────────────
rule eukrep:
    input:
        contigs = f"{MOD2_OUT}/12_filtered/contigs_filtered.fasta"
    output:
        euk     = f"{MOD4_OUT}/17a_eukrep/euk_contigs.fasta",
        prok    = f"{MOD4_OUT}/17a_eukrep/prok_contigs.fasta",
    params:
        min_len = config["mod4_binning"].get("eukrep_min_len", 1000),
    log: f"{MOD4_OUT}/logs/17a_eukrep/eukrep.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {MOD4_OUT}/17a_eukrep

        EukRep \
            -i {input.contigs} \
            -o {output.euk} \
            --prokarya {output.prok} \
            --min {params.min_len} \
            &> {log}
        """

rule tiara:
    input:
        euk = f"{MOD4_OUT}/17a_eukrep/euk_contigs.fasta"
    output:
        tsv  = f"{MOD4_OUT}/17a_eukrep/tiara_classifications.tsv",
        fung = f"{MOD4_OUT}/17a_eukrep/fungi_contigs.fasta",
    params:
        min_len = config["mod4_binning"].get("tiara_min_len", 1000),
    threads: THREADS
    log: f"{MOD4_OUT}/logs/17a_eukrep/tiara.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail

        tiara \
            -i {input.euk} \
            -o {output.tsv} \
            --min_len {params.min_len} \
            --threads {threads} \
            &> {log}

        # Extract contigs classified as fungi or eukarya
        python3 -c "
import pandas as pd
from Bio import SeqIO

df = pd.read_csv('{output.tsv}', sep='\t')
keep = set(df[df['class_fst_stage'].isin(['fungi','eukarya'])]['sequence_id'])
recs = [r for r in SeqIO.parse('{input.euk}', 'fasta') if r.id in keep]
SeqIO.write(recs, '{output.fung}', 'fasta')
print(f'Fungi/euk contigs kept: {{len(recs)}}')
" &>> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [17b] ENSEMBLE BINNING — MetaBAT2 + MaxBin2 + CONCOCT
# Runs only when binning_mode: ensemble
# ─────────────────────────────────────────────────────────────────────────────
if BINNING_MODE == "ensemble":

    rule metabat2:
        input:
            contigs  = f"{MOD4_OUT}/17a_eukrep/fungi_contigs.fasta",
            coverage = f"{MOD2_OUT}/15_coverm/coverage_table.tsv",
        output:
            done = f"{MOD4_OUT}/17b_ensemble/metabat2/binning_done.flag",
        params:
            outdir  = f"{MOD4_OUT}/17b_ensemble/metabat2",
            min_len = config["mod4_binning"].get("metabat2_min_len", 1500),
        threads: THREADS
        log: f"{MOD4_OUT}/logs/17b_ensemble/metabat2.log"
        conda: "../envs/mod4_binning.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir}

            # Convert CoverM table to MetaBAT2 depth format
            python3 -c "
import pandas as pd
df = pd.read_csv('{input.coverage}', sep='\t')
df.to_csv('{params.outdir}/depth.txt', sep='\t', index=False)
"
            metabat2 \
                -i {input.contigs} \
                -a {params.outdir}/depth.txt \
                -o {params.outdir}/bin \
                --minContig {params.min_len} \
                --numThreads {threads} \
                &> {log}

            touch {output.done}
            """

    rule maxbin2:
        input:
            contigs  = f"{MOD4_OUT}/17a_eukrep/fungi_contigs.fasta",
            coverage = f"{MOD2_OUT}/15_coverm/coverage_table.tsv",
        output:
            done = f"{MOD4_OUT}/17b_ensemble/maxbin2/binning_done.flag",
        params:
            outdir = f"{MOD4_OUT}/17b_ensemble/maxbin2",
        threads: THREADS
        log: f"{MOD4_OUT}/logs/17b_ensemble/maxbin2.log"
        conda: "../envs/mod4_binning.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir}

            # Extract coverage values for MaxBin2
            python3 -c "
import pandas as pd
df = pd.read_csv('{input.coverage}', sep='\t')
cov_cols = [c for c in df.columns if c != 'Contig']
for col in cov_cols:
    df[['Contig', col]].to_csv(
        '{params.outdir}/' + col.replace(' ','_') + '_cov.txt',
        sep='\t', index=False, header=False
    )
"
            COV_FILES=$(ls {params.outdir}/*_cov.txt | tr '\n' ',')
            COV_FILES=${{COV_FILES%,}}

            run_MaxBin.pl \
                -contig {input.contigs} \
                -abund_list $COV_FILES \
                -out {params.outdir}/bin \
                -thread {threads} \
                &> {log}

            touch {output.done}
            """

    rule concoct:
        input:
            contigs  = f"{MOD4_OUT}/17a_eukrep/fungi_contigs.fasta",
            bams     = expand(f"{MOD2_OUT}/14_mapping/{{sample}}.sorted.bam", sample=SAMPLES),
        output:
            done = f"{MOD4_OUT}/17b_ensemble/concoct/binning_done.flag",
        params:
            outdir   = f"{MOD4_OUT}/17b_ensemble/concoct",
            chunk    = config["mod4_binning"].get("concoct_chunk", 10000),
            clusters = config["mod4_binning"].get("concoct_clusters", 400),
        threads: THREADS
        log: f"{MOD4_OUT}/logs/17b_ensemble/concoct.log"
        conda: "../envs/mod4_binning.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir}

            # Cut contigs into chunks for CONCOCT
            cut_up_fasta.py {input.contigs} \
                -c {params.chunk} -o 0 \
                --merge_last \
                -b {params.outdir}/contigs_10k.bed \
                > {params.outdir}/contigs_10k.fasta

            # Generate coverage table
            concoct_coverage_table.py \
                {params.outdir}/contigs_10k.bed \
                {input.bams} \
                > {params.outdir}/coverage_table.tsv

            # Run CONCOCT
            concoct \
                --composition_file {params.outdir}/contigs_10k.fasta \
                --coverage_file {params.outdir}/coverage_table.tsv \
                --clusters {params.clusters} \
                --threads {threads} \
                --basename {params.outdir}/concoct \
                &> {log}

            # Merge sub-contig clustering
            merge_cutup_clustering.py \
                {params.outdir}/concoct_clustering_gt1000.csv \
                > {params.outdir}/clustering_merged.csv

            # Extract bins as fasta
            extract_fasta_bins.py \
                {input.contigs} \
                {params.outdir}/clustering_merged.csv \
                --output_path {params.outdir}/bins/

            touch {output.done}
            """

    # ─────────────────────────────────────────────────────────────────────────
    # [17c] DAS_Tool — picks best non-overlapping bins from ensemble
    # ─────────────────────────────────────────────────────────────────────────
    rule dastool:
        input:
            contigs  = f"{MOD4_OUT}/17a_eukrep/fungi_contigs.fasta",
            metabat2 = f"{MOD4_OUT}/17b_ensemble/metabat2/binning_done.flag",
            maxbin2  = f"{MOD4_OUT}/17b_ensemble/maxbin2/binning_done.flag",
            concoct  = f"{MOD4_OUT}/17b_ensemble/concoct/binning_done.flag",
        output:
            done     = f"{MOD4_OUT}/17c_dastool/dastool_done.flag",
            summary  = f"{MOD4_OUT}/17c_dastool/DASTool_summary.txt",
        params:
            outdir       = f"{MOD4_OUT}/17c_dastool",
            score_thresh = config["mod4_binning"].get("dastool_score_threshold", 0.5),
            metabat2_dir = f"{MOD4_OUT}/17b_ensemble/metabat2",
            maxbin2_dir  = f"{MOD4_OUT}/17b_ensemble/maxbin2",
            concoct_dir  = f"{MOD4_OUT}/17b_ensemble/concoct/bins",
        threads: THREADS
        log: f"{MOD4_OUT}/logs/17c_dastool/dastool.log"
        conda: "../envs/mod4_binning.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir}

            # Generate scaffold-to-bin files for each binner
            Fasta_to_Scaffolds2Bin.sh -i {params.metabat2_dir} -e fa \
                > {params.outdir}/metabat2_s2b.tsv
            Fasta_to_Scaffolds2Bin.sh -i {params.maxbin2_dir} -e fasta \
                > {params.outdir}/maxbin2_s2b.tsv
            Fasta_to_Scaffolds2Bin.sh -i {params.concoct_dir} -e fa \
                > {params.outdir}/concoct_s2b.tsv

            DAS_Tool \
                -i {params.outdir}/metabat2_s2b.tsv,{params.outdir}/maxbin2_s2b.tsv,{params.outdir}/concoct_s2b.tsv \
                -l metabat2,maxbin2,concoct \
                -c {input.contigs} \
                -o {params.outdir}/DASTool \
                --score_threshold {params.score_thresh} \
                --threads {threads} \
                --write_bins \
                &> {log}

            touch {output.done}
            """

# ─────────────────────────────────────────────────────────────────────────────
# [17b/17c] HYPERBIN MODE — HyphaeBin only
# PLACEHOLDER — to be implemented
# ─────────────────────────────────────────────────────────────────────────────
if BINNING_MODE == "hyperbin":

    rule hyperbin:
        input:
            contigs  = f"{MOD4_OUT}/17a_eukrep/fungi_contigs.fasta",
            coverage = f"{MOD2_OUT}/15_coverm/coverage_table.tsv",
        output:
            done     = f"{MOD4_OUT}/17d_hyperbin/binning_done.flag",
            bins_dir = directory(f"{MOD4_OUT}/17d_hyperbin/bins"),
        params:
            outdir   = f"{MOD4_OUT}/17d_hyperbin",
        threads: THREADS
        log: f"{MOD4_OUT}/logs/17d_hyperbin/hyperbin.log"
        conda: "../envs/mod4_binning.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p {params.outdir}/bins

            # ================================================================
            # HyphaeBin — PLACEHOLDER
            # Algorithm to be implemented
            # Inputs available:
            #   - Fungi contigs : {input.contigs}
            #   - Coverage table: {input.coverage}
            #   - Output bins   : {params.outdir}/bins/
            # ================================================================
            echo "HyphaeBin placeholder — not yet implemented" | tee {log}
            echo "Place bin FASTA files in {params.outdir}/bins/"

            # Remove this pass block and implement HyphaeBin here
            pass

            touch {output.done}
            """

# ─────────────────────────────────────────────────────────────────────────────
# Helper: get bins dir based on mode
# ─────────────────────────────────────────────────────────────────────────────
def get_bins_done(wildcards=None):
    if BINNING_MODE == "ensemble":
        return f"{MOD4_OUT}/17c_dastool/dastool_done.flag"
    else:
        return f"{MOD4_OUT}/17d_hyperbin/binning_done.flag"

def get_bins_dir(wildcards=None):
    if BINNING_MODE == "ensemble":
        return f"{MOD4_OUT}/17c_dastool/DASTool_bins"
    else:
        return f"{MOD4_OUT}/17d_hyperbin/bins"

# ─────────────────────────────────────────────────────────────────────────────
# [17e] EukCC + BUSCO + GUNC — MAG quality
# EukCC : fungi-specific completeness/contamination
# BUSCO : universal single-copy orthologs
# GUNC  : chimerism detection
# ─────────────────────────────────────────────────────────────────────────────
rule eukcc:
    input:
        bins_done = get_bins_done,
    output:
        f"{MOD4_OUT}/17e_quality/eukcc/eukcc_summary.tsv"
    params:
        bins_dir = get_bins_dir,
        outdir   = f"{MOD4_OUT}/17e_quality/eukcc",
        db       = EUKCC_DB,
    threads: THREADS
    log: f"{MOD4_OUT}/logs/17e_quality/eukcc.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        eukcc folder \
            --db {params.db} \
            --out {params.outdir} \
            --threads {threads} \
            {params.bins_dir} \
            &> {log}
        """

rule busco:
    input:
        bins_done = get_bins_done,
    output:
        f"{MOD4_OUT}/17e_quality/busco/busco_summary.tsv"
    params:
        bins_dir = get_bins_dir,
        outdir   = f"{MOD4_OUT}/17e_quality/busco",
        lineage  = BUSCO_LINEAGE,
    threads: THREADS
    log: f"{MOD4_OUT}/logs/17e_quality/busco.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        for BIN in {params.bins_dir}/*.fa; do
            NAME=$(basename $BIN .fa)
            busco \
                -i $BIN \
                -o {params.outdir}/$NAME \
                -l {params.lineage} \
                -m genome \
                --cpu {threads} \
                --augustus \
                --augustus_species botrytis_cinerea \
                -f \
                &>> {log}
        done

        # Merge summaries
        python3 -c "
import os, glob, pandas as pd
rows = []
for f in glob.glob('{params.outdir}/*/short_summary*.txt'):
    name = os.path.basename(os.path.dirname(f))
    with open(f) as fh:
        for line in fh:
            if line.strip().startswith('C:'):
                rows.append({{'bin': name, 'busco_summary': line.strip()}})
pd.DataFrame(rows).to_csv('{output}', sep='\t', index=False)
"
        """

rule gunc:
    input:
        bins_done = get_bins_done,
    output:
        f"{MOD4_OUT}/17e_quality/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv"
    params:
        bins_dir = get_bins_dir,
        outdir   = f"{MOD4_OUT}/17e_quality/gunc",
    threads: THREADS
    log: f"{MOD4_OUT}/logs/17e_quality/gunc.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        gunc run \
            --input_dir {params.bins_dir} \
            --out_dir {params.outdir} \
            --threads {threads} \
            &> {log}
        """

rule quality_summary:
    input:
        eukcc = f"{MOD4_OUT}/17e_quality/eukcc/eukcc_summary.tsv",
        busco = f"{MOD4_OUT}/17e_quality/busco/busco_summary.tsv",
        gunc  = f"{MOD4_OUT}/17e_quality/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv",
    output:
        f"{MOD4_OUT}/17e_quality/eukcc_done.flag",
    params:
        min_completeness   = config["mod4_binning"]["min_completeness"],
        max_contamination  = config["mod4_binning"]["max_contamination"],
        outdir             = f"{MOD4_OUT}/17e_quality",
    log: f"{MOD4_OUT}/logs/17e_quality/summary.log"
    run:
        import pandas as pd
        eukcc = pd.read_csv(input.eukcc, sep="\t")
        hq    = eukcc[
            (eukcc["completeness"] >= params.min_completeness) &
            (eukcc["contamination"] <= params.max_contamination)
        ]
        hq.to_csv(f"{params.outdir}/HQ_MAGs.tsv", sep="\t", index=False)
        with open(output[0], "w") as f:
            f.write(f"HQ MAGs: {len(hq)}\n")

# ─────────────────────────────────────────────────────────────────────────────
# [18] MAG TAXONOMY — Kaiju + EPA-ng + IQ-TREE
# Kaiju: protein-level MAG taxonomy
# EPA-ng: phylogenetic placement into reference tree
# IQ-TREE: phylogenetic tree construction
# ─────────────────────────────────────────────────────────────────────────────
rule mag_taxonomy_kaiju:
    input:
        quality_flag = f"{MOD4_OUT}/17e_quality/eukcc_done.flag",
    output:
        f"{MOD4_OUT}/18_mag_taxonomy/kaiju_mag.tsv"
    params:
        bins_dir = get_bins_dir,
        db       = KAIJU_DB,
        nodes    = KAIJU_NODES,
        outdir   = f"{MOD4_OUT}/18_mag_taxonomy",
    threads: THREADS
    log: f"{MOD4_OUT}/logs/18_mag_taxonomy/kaiju.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        for BIN in {params.bins_dir}/*.fa; do
            NAME=$(basename $BIN .fa)
            kaiju \
                -t {params.nodes} \
                -f {params.db} \
                -i $BIN \
                -o {params.outdir}/${{NAME}}_kaiju.out \
                -z {threads} &>> {log}

            kaiju2table \
                -t {params.nodes} \
                -n $(dirname {params.nodes})/names.dmp \
                -r genus \
                -o {params.outdir}/${{NAME}}_kaiju_table.tsv \
                {params.outdir}/${{NAME}}_kaiju.out &>> {log}
        done

        # Merge all MAG taxonomy
        python3 -c "
import glob, pandas as pd
dfs = []
for f in glob.glob('{params.outdir}/*_kaiju_table.tsv'):
    df = pd.read_csv(f, sep='\t')
    df['bin'] = f.replace('_kaiju_table.tsv','')
    dfs.append(df)
pd.concat(dfs).to_csv('{output}', sep='\t', index=False)
" &>> {log}
        """

rule mag_phylogeny:
    input:
        f"{MOD4_OUT}/18_mag_taxonomy/kaiju_mag.tsv"
    output:
        tree = f"{MOD4_OUT}/18_mag_taxonomy/phylogeny/iqtree.treefile",
        tsv  = f"{MOD4_OUT}/18_mag_taxonomy/mag_taxonomy.tsv",
    params:
        bins_dir = get_bins_dir,
        outdir   = f"{MOD4_OUT}/18_mag_taxonomy/phylogeny",
        ref_tree = config["mod4_binning"].get("reference_tree", ""),
        ref_msa  = config["mod4_binning"].get("reference_msa", ""),
    threads: THREADS
    log: f"{MOD4_OUT}/logs/18_mag_taxonomy/phylogeny.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        # EPA-ng — place MAG sequences into reference tree
        if [ -n "{params.ref_tree}" ] && [ -f "{params.ref_tree}" ]; then
            epa-ng \
                --tree {params.ref_tree} \
                --ref-msa {params.ref_msa} \
                --query {params.bins_dir}/*.fa \
                --outdir {params.outdir} \
                --threads {threads} \
                &>> {log}
        fi

        # IQ-TREE — build phylogenetic tree
        iqtree2 \
            -s {params.bins_dir} \
            --prefix {params.outdir}/iqtree \
            -T {threads} \
            -m TEST \
            -B 1000 \
            &>> {log}

        # Copy kaiju taxonomy as final MAG taxonomy
        cp {input} {output.tsv}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [19] MAG ANNOTATION — Augustus + eggNOG + dbCAN + antiSMASH
# Augustus : gene prediction (fungi-trained models)
# eggNOG   : functional annotation
# dbCAN    : CAZyme annotation
# antiSMASH: biosynthetic gene clusters
# ─────────────────────────────────────────────────────────────────────────────
rule augustus:
    input:
        quality_flag = f"{MOD4_OUT}/17e_quality/eukcc_done.flag",
    output:
        done = f"{MOD4_OUT}/19_mag_annotation/augustus_done.flag"
    params:
        bins_dir   = get_bins_dir,
        outdir     = f"{MOD4_OUT}/19_mag_annotation/augustus",
        species    = config["mod4_binning"].get("augustus_species", "botrytis_cinerea"),
        config_dir = AUGUSTUS_CFG,
    threads: THREADS
    log: f"{MOD4_OUT}/logs/19_annotation/augustus.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        export AUGUSTUS_CONFIG_PATH={params.config_dir}

        for BIN in {params.bins_dir}/*.fa; do
            NAME=$(basename $BIN .fa)
            augustus \
                --species={params.species} \
                --AUGUSTUS_CONFIG_PATH={params.config_dir} \
                --outfile={params.outdir}/${{NAME}}.gff \
                --errfile={params.outdir}/${{NAME}}.err \
                $BIN &>> {log}

            # Convert GFF to protein FASTA
            getAnnoFasta.pl {params.outdir}/${{NAME}}.gff --seqfile=$BIN &>> {log}
        done

        touch {output.done}
        """

rule eggnog_mags:
    input:
        f"{MOD4_OUT}/19_mag_annotation/augustus_done.flag"
    output:
        done = f"{MOD4_OUT}/19_mag_annotation/eggnog_done.flag"
    params:
        augustus_dir = f"{MOD4_OUT}/19_mag_annotation/augustus",
        outdir       = f"{MOD4_OUT}/19_mag_annotation/eggnog",
        db           = config["mod5_annotation"]["eggnog_db"],
    threads: THREADS
    log: f"{MOD4_OUT}/logs/19_annotation/eggnog.log"
    conda: "../envs/mod5_annotation.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        for PROT in {params.augustus_dir}/*.aa; do
            NAME=$(basename $PROT .aa)
            emapper.py \
                -i $PROT \
                -o {params.outdir}/$NAME \
                --data_dir {params.db} \
                --cpu {threads} \
                --override \
                &>> {log}
        done

        touch {output.done}
        """

rule dbcan_mags:
    input:
        f"{MOD4_OUT}/19_mag_annotation/augustus_done.flag"
    output:
        done = f"{MOD4_OUT}/19_mag_annotation/dbcan_done.flag"
    params:
        augustus_dir = f"{MOD4_OUT}/19_mag_annotation/augustus",
        outdir       = f"{MOD4_OUT}/19_mag_annotation/dbcan",
        db           = DBCAN_DB,
    threads: THREADS
    log: f"{MOD4_OUT}/logs/19_annotation/dbcan.log"
    conda: "../envs/mod5_annotation.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        for PROT in {params.augustus_dir}/*.aa; do
            NAME=$(basename $PROT .aa)
            run_dbcan.py \
                $PROT protein \
                --out_dir {params.outdir}/$NAME \
                --db_dir {params.db} \
                --cpu {threads} \
                &>> {log}
        done

        touch {output.done}
        """

rule antismash:
    input:
        f"{MOD4_OUT}/19_mag_annotation/augustus_done.flag"
    output:
        done = f"{MOD4_OUT}/19_mag_annotation/antismash_done.flag"
    params:
        bins_dir    = get_bins_dir,
        outdir      = f"{MOD4_OUT}/19_mag_annotation/antismash",
        db          = ANTISMASH_DB,
        taxon       = config["mod4_binning"].get("antismash_taxon", "fungi"),
    threads: THREADS
    log: f"{MOD4_OUT}/logs/19_annotation/antismash.log"
    conda: "../envs/mod4_binning.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        for BIN in {params.bins_dir}/*.fa; do
            NAME=$(basename $BIN .fa)
            antismash \
                $BIN \
                --output-dir {params.outdir}/$NAME \
                --taxon {params.taxon} \
                --databases {params.db} \
                --cpus {threads} \
                --genefinding-tool glimmerhmm \
                &>> {log}
        done

        touch {output.done}
        """

rule annotation_done:
    input:
        eggnog   = f"{MOD4_OUT}/19_mag_annotation/eggnog_done.flag",
        dbcan    = f"{MOD4_OUT}/19_mag_annotation/dbcan_done.flag",
        antismash= f"{MOD4_OUT}/19_mag_annotation/antismash_done.flag",
    output:
        f"{MOD4_OUT}/19_mag_annotation/annotation_done.flag"
    shell:
        "touch {output}"
