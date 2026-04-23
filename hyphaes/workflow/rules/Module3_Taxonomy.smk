# =============================================================================
# Module 3 — Taxonomy (read level)
# [16a] Kraken2 — PlusPF DB (community overview, bacteria/fungi ratio)
# [16b] Kaiju — nr_euk DB (protein-level, detects novel fungi)
# [16c] Bracken — species-level abundance re-estimation
# [16d] Diversity analysis — alpha/beta/PCoA via phyloseq/R
# =============================================================================

onstart:
    check_mod1_done()

# ── Target ────────────────────────────────────────────────────────────────────
rule mod3_all:
    input:
        expand(f"{MOD3_OUT}/16a_kraken2/{{sample}}.report", sample=SAMPLES),
        expand(f"{MOD3_OUT}/16b_kaiju/{{sample}}_kaiju.tsv", sample=SAMPLES),
        expand(f"{MOD3_OUT}/16c_bracken/{{sample}}.bracken", sample=SAMPLES),
        f"{MOD3_OUT}/16d_diversity/alpha_diversity.tsv",
        f"{MOD3_OUT}/16d_diversity/beta_diversity_pcoa.pdf",

# ─────────────────────────────────────────────────────────────────────────────
# [16a] KRAKEN2 — PlusPF database
# Community overview + bacteria/fungi ratio
# PlusPF = standard + protozoa + fungi database
# ─────────────────────────────────────────────────────────────────────────────
rule kraken2:
    input:
        r1 = get_clean_r1,
        r2 = get_clean_r2,
    output:
        report  = f"{MOD3_OUT}/16a_kraken2/{{sample}}.report",
        out     = f"{MOD3_OUT}/16a_kraken2/{{sample}}.kraken2",
    params:
        db         = KRAKEN2_DB,
        confidence = config["mod3_taxonomy"].get("kraken2_confidence", 0.1),
    threads: THREADS
    log: f"{MOD3_OUT}/logs/16a_kraken2/{{sample}}.log"
    conda: "../envs/mod3_taxonomy.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {MOD3_OUT}/16a_kraken2

        kraken2 \
            --db {params.db} \
            --paired \
            --threads {threads} \
            --confidence {params.confidence} \
            --report {output.report} \
            --output {output.out} \
            --gzip-compressed \
            {input.r1} {input.r2} \
            &> {log}
        """

rule kraken2_multiqc:
    input:
        expand(f"{MOD3_OUT}/16a_kraken2/{{sample}}.report", sample=SAMPLES)
    output:
        f"{MOD3_OUT}/16a_kraken2/multiqc_report.html"
    log: f"{MOD3_OUT}/logs/16a_kraken2/multiqc.log"
    conda: "../envs/mod3_taxonomy.yaml"
    shell:
        """
        multiqc {MOD3_OUT}/16a_kraken2 \
            -o {MOD3_OUT}/16a_kraken2 \
            --force &> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [16b] KAIJU — nr_euk database
# Protein-level classification — detects novel fungi missed by Kraken2
# nr_euk = NCBI nr restricted to eukaryotes
# ─────────────────────────────────────────────────────────────────────────────
rule kaiju:
    input:
        r1 = get_clean_r1,
        r2 = get_clean_r2,
    output:
        raw = f"{MOD3_OUT}/16b_kaiju/{{sample}}_kaiju.out",
        tsv = f"{MOD3_OUT}/16b_kaiju/{{sample}}_kaiju.tsv",
    params:
        db    = KAIJU_DB,
        nodes = KAIJU_NODES,
        evalue= config["mod3_taxonomy"].get("kaiju_evalue", "0.01"),
    threads: THREADS
    log: f"{MOD3_OUT}/logs/16b_kaiju/{{sample}}.log"
    conda: "../envs/mod3_taxonomy.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {MOD3_OUT}/16b_kaiju

        # Run Kaiju
        kaiju \
            -t {params.nodes} \
            -f {params.db} \
            -i <(zcat {input.r1}) \
            -j <(zcat {input.r2}) \
            -o {output.raw} \
            -z {threads} \
            -E {params.evalue} \
            &> {log}

        # Add taxonomy names to output
        kaiju-addTaxonNames \
            -t {params.nodes} \
            -n $(dirname {params.nodes})/names.dmp \
            -i {output.raw} \
            -o {output.tsv} \
            -r superkingdom,phylum,class,order,family,genus,species \
            &>> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [16c] BRACKEN — species-level abundance re-estimation
# Re-estimates species abundances from Kraken2 reports
# Must match k-mer length used to build Kraken2 DB (default 150)
# ─────────────────────────────────────────────────────────────────────────────
rule bracken:
    input:
        report = f"{MOD3_OUT}/16a_kraken2/{{sample}}.report",
    output:
        bracken = f"{MOD3_OUT}/16c_bracken/{{sample}}.bracken",
        report  = f"{MOD3_OUT}/16c_bracken/{{sample}}_bracken.report",
    params:
        db        = KRAKEN2_DB,
        read_len  = config["mod3_taxonomy"].get("bracken_read_len", 150),
        level     = config["mod3_taxonomy"].get("bracken_level", "S"),
        threshold = config["mod3_taxonomy"].get("bracken_threshold", 10),
    log: f"{MOD3_OUT}/logs/16c_bracken/{{sample}}.log"
    conda: "../envs/mod3_taxonomy.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p {MOD3_OUT}/16c_bracken

        bracken \
            -d {params.db} \
            -i {input.report} \
            -o {output.bracken} \
            -w {output.report} \
            -r {params.read_len} \
            -l {params.level} \
            -t {params.threshold} \
            &> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [16d] DIVERSITY ANALYSIS — phyloseq/R
# Alpha diversity (Shannon, Simpson, Chao1)
# Beta diversity (Bray-Curtis PCoA)
# ─────────────────────────────────────────────────────────────────────────────
rule diversity_analysis:
    input:
        bracken = expand(f"{MOD3_OUT}/16c_bracken/{{sample}}.bracken", sample=SAMPLES),
    output:
        alpha   = f"{MOD3_OUT}/16d_diversity/alpha_diversity.tsv",
        beta    = f"{MOD3_OUT}/16d_diversity/beta_diversity_pcoa.pdf",
        combined= f"{MOD3_OUT}/16d_diversity/abundance_table.tsv",
    params:
        indir   = f"{MOD3_OUT}/16c_bracken",
        outdir  = f"{MOD3_OUT}/16d_diversity",
        samples = ",".join(SAMPLES),
    log: f"{MOD3_OUT}/logs/16d_diversity/diversity.log"
    conda: "../envs/mod3_taxonomy.yaml"
    script:
        "../scripts/diversity_analysis.R"
