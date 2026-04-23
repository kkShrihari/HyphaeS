# =============================================================================
# Module 2 — Assembly
# [11] Co-assembly — metaSPAdes (all samples pooled)
# [12] Filter contigs — min 1000 bp
# [13] Assembly QC — QUAST
# [14] Read mapping — Bowtie2 (per-sample BAMs)
# [15] Coverage calculation — CoverM
# =============================================================================

onstart:
    check_mod1_done()

# ── Target ────────────────────────────────────────────────────────────────────
rule mod2_all:
    input:
        f"{MOD2_OUT}/13_quast/report.html",
        f"{MOD2_OUT}/15_coverm/coverage_table.tsv",
        f"{MOD2_OUT}/plots/05_assembly_dashboard.png",

# ─────────────────────────────────────────────────────────────────────────────
# [11] CO-ASSEMBLY — metaSPAdes
# All samples pooled — better fungi contig recovery than per-sample assembly
# memory set in config (gmbs17 has ~535 GB)
# ─────────────────────────────────────────────────────────────────────────────
rule co_assembly:
    input:
        r1 = get_all_clean_r1,
        r2 = get_all_clean_r2,
    output:
        contigs   = f"{MOD2_OUT}/11_spades/contigs.fasta",
        scaffolds = f"{MOD2_OUT}/11_spades/scaffolds.fasta",
    params:
        outdir  = f"{MOD2_OUT}/11_spades",
        memory  = config["mod2_assembly"]["spades_memory_gb"],
        kmer    = config["mod2_assembly"].get("kmer", "21,33,55,77"),
    threads: THREADS
    log: f"{MOD2_OUT}/logs/11_spades/co_assembly.log"
    conda: "../envs/mod2_assembly.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        # Merge all sample reads for co-assembly
        R1={params.outdir}/merged_R1.fq.gz
        R2={params.outdir}/merged_R2.fq.gz
        cat {input.r1} > $R1
        cat {input.r2} > $R2

        spades.py \
            --meta \
            -1 $R1 -2 $R2 \
            -o {params.outdir} \
            -k {params.kmer} \
            --memory {params.memory} \
            --threads {threads} \
            &> {log}

        rm -f $R1 $R2
        """

# ─────────────────────────────────────────────────────────────────────────────
# [12] FILTER CONTIGS — min 1000 bp
# Lower threshold than standard (2500 bp) to keep more fungi contigs
# ─────────────────────────────────────────────────────────────────────────────
rule filter_contigs:
    input:
        f"{MOD2_OUT}/11_spades/contigs.fasta"
    output:
        f"{MOD2_OUT}/12_filtered/contigs_filtered.fasta"
    params:
        min_len = config["mod2_assembly"]["min_contig_len"],
    log: f"{MOD2_OUT}/logs/12_filter/filter.log"
    conda: "../envs/mod2_assembly.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {MOD2_OUT}/12_filtered

        seqkit seq \
            --min-len {params.min_len} \
            --out-file {output} \
            {input} \
            &> {log}

        echo "Contigs kept: $(grep -c '>' {output})" >> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [13] ASSEMBLY QC — QUAST
# Reports N50, total length, % reads mapped
# ─────────────────────────────────────────────────────────────────────────────
rule assembly_qc:
    input:
        contigs = f"{MOD2_OUT}/12_filtered/contigs_filtered.fasta",
        r1      = get_all_clean_r1,
        r2      = get_all_clean_r2,
    output:
        html = f"{MOD2_OUT}/13_quast/report.html",
        tsv  = f"{MOD2_OUT}/13_quast/report.tsv",
    params:
        outdir = f"{MOD2_OUT}/13_quast",
        min_len = config["mod2_assembly"]["min_contig_len"],
    threads: THREADS
    log: f"{MOD2_OUT}/logs/13_quast/quast.log"
    conda: "../envs/mod2_assembly.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        R1={params.outdir}/merged_R1.fq.gz
        R2={params.outdir}/merged_R2.fq.gz
        cat {input.r1} > $R1
        cat {input.r2} > $R2

        quast.py \
            {input.contigs} \
            --pe1 $R1 --pe2 $R2 \
            --output-dir {params.outdir} \
            --threads {threads} \
            --min-contig {params.min_len} \
            --no-icarus \
            &> {log}

        rm -f $R1 $R2
        """

# ─────────────────────────────────────────────────────────────────────────────
# [14] READ MAPPING — Bowtie2
# Map each sample back to co-assembly — per-sample BAM files
# BAMs are required by CoverM for per-contig coverage
# ─────────────────────────────────────────────────────────────────────────────
rule build_assembly_index:
    input:
        f"{MOD2_OUT}/12_filtered/contigs_filtered.fasta"
    output:
        multiext(f"{MOD2_OUT}/14_mapping/assembly_idx",
                 ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                 ".rev.1.bt2", ".rev.2.bt2")
    params:
        prefix = f"{MOD2_OUT}/14_mapping/assembly_idx"
    threads: THREADS
    log: f"{MOD2_OUT}/logs/14_mapping/bowtie2_build.log"
    conda: "../envs/mod2_assembly.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {MOD2_OUT}/14_mapping
        bowtie2-build --threads {threads} {input} {params.prefix} &> {log}
        """

rule map_reads:
    input:
        r1  = get_clean_r1,
        r2  = get_clean_r2,
        idx = multiext(f"{MOD2_OUT}/14_mapping/assembly_idx",
                       ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2")
    output:
        bam = f"{MOD2_OUT}/14_mapping/{{sample}}.sorted.bam",
        bai = f"{MOD2_OUT}/14_mapping/{{sample}}.sorted.bam.bai",
    params:
        prefix = f"{MOD2_OUT}/14_mapping/assembly_idx"
    threads: THREADS
    log: f"{MOD2_OUT}/logs/14_mapping/{{sample}}.log"
    conda: "../envs/mod2_assembly.yaml"
    shell:
        """
        set -euo pipefail

        bowtie2 \
            -x {params.prefix} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            --sensitive \
            2>> {log} \
        | samtools view -b -F 4 \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [15] COVERAGE CALCULATION — CoverM
# Coverage per contig per sample — used as binning signal
# ─────────────────────────────────────────────────────────────────────────────
rule coverm:
    input:
        bams    = expand(f"{MOD2_OUT}/14_mapping/{{sample}}.sorted.bam", sample=SAMPLES),
        contigs = f"{MOD2_OUT}/12_filtered/contigs_filtered.fasta",
    output:
        f"{MOD2_OUT}/15_coverm/coverage_table.tsv"
    params:
        outdir = f"{MOD2_OUT}/15_coverm",
        method = config["mod2_assembly"].get("coverm_method", "mean"),
    threads: THREADS
    log: f"{MOD2_OUT}/logs/15_coverm/coverm.log"
    conda: "../envs/mod2_assembly.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        coverm contig \
            --bam-files {input.bams} \
            --methods {params.method} \
            --threads {threads} \
            --output-file {output} \
            &> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [16] ASSEMBLY PLOTS
# ─────────────────────────────────────────────────────────────────────────────
rule mod2_plots:
    input:
        contigs  = f"{MOD2_OUT}/12_filtered/contigs_filtered.fasta",
        coverage = f"{MOD2_OUT}/15_coverm/coverage_table.tsv",
        bams     = expand(f"{MOD2_OUT}/14_mapping/{{sample}}.sorted.bam", sample=SAMPLES),
    output:
        f"{MOD2_OUT}/plots/05_assembly_dashboard.png"
    conda: "../envs/mod2_assembly.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {MOD2_OUT}/plots
        python {PLOT_SCRIPT_MOD2} \
            --results {MOD2_OUT} \
            --outdir {MOD2_OUT}/plots
        """
