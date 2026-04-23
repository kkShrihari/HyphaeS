# =============================================================================
# Module 1 — Quality Control
# [0]  Input validation
# [1]  FastQC + MultiQC (raw)
# [2]  Trim Galore — Q15, min 50 bp
# [3]  PhiX removal — BBDuk k=31
# [3b] Low-complexity filter — BBDuk entropy=0.6
# [4]  Host removal — Bowtie2 + samtools
# [5]  FastQC + MultiQC (post-clean)
# [6]  EukDetect — fungi screen
# [7]  Nonpareil3 — coverage sufficiency
# [8]  seqkit stats
# [9]  Retention summary
# [99] Merged MultiQC final report
#
# NOTE: SAMPLES, THREADS, MOD1_OUT, HOST_IDX, PHIX_REF, EUKDETECT_CFG
#       are all defined in common.smk — do NOT redefine here
# =============================================================================

# ── Target ────────────────────────────────────────────────────────────────────
rule mod1_all:
    input:
        f"{MOD1_OUT}/01_rawQC/multiqc_report.html",
        f"{MOD1_OUT}/05_cleanQC/multiqc_report.html",
        f"{MOD1_OUT}/99_multiqc_final/multiqc_report.html",
        expand(f"{MOD1_OUT}/09_retention/{{sample}}_retention.tsv", sample=SAMPLES),
        f"{MOD1_OUT}/plots/05_QC_dashboard.png",

# ─────────────────────────────────────────────────────────────────────────────
# [0] INPUT VALIDATION
# ─────────────────────────────────────────────────────────────────────────────
rule validate_input:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        flag = f"{MOD1_OUT}/00_validation/{{sample}}.ok"
    log: f"{MOD1_OUT}/logs/00_validation/{{sample}}.log"
    run:
        import gzip, os
        os.makedirs(os.path.dirname(output.flag), exist_ok=True)

        def count_reads(f):
            n = 0
            with gzip.open(f, "rt") as fh:
                for _ in fh:
                    n += 1
            return n // 4

        r1_n = count_reads(input.r1)
        r2_n = count_reads(input.r2)

        with open(log[0], "w") as lf:
            lf.write(f"{wildcards.sample}: R1={r1_n} R2={r2_n}\n")

        if r1_n != r2_n:
            raise ValueError(f"FAIL {wildcards.sample}: R1={r1_n} R2={r2_n} — UNPAIRED")
        if r1_n == 0:
            raise ValueError(f"FAIL {wildcards.sample}: R1 is EMPTY")

        with open(output.flag, "w") as fh:
            fh.write(f"{wildcards.sample}\tR1={r1_n}\tR2={r2_n}\tPAIRED_OK\n")

# ─────────────────────────────────────────────────────────────────────────────
# [1] FastQC + MultiQC — RAW
# ─────────────────────────────────────────────────────────────────────────────
rule fastqc_raw:
    input:
        flag = f"{MOD1_OUT}/00_validation/{{sample}}.ok",
        r1   = get_r1,
        r2   = get_r2,
    output:
        html_r1 = f"{MOD1_OUT}/01_rawQC/{{sample}}_R1_fastqc.html",
        html_r2 = f"{MOD1_OUT}/01_rawQC/{{sample}}_R2_fastqc.html",
        zip_r1  = f"{MOD1_OUT}/01_rawQC/{{sample}}_R1_fastqc.zip",
        zip_r2  = f"{MOD1_OUT}/01_rawQC/{{sample}}_R2_fastqc.zip",
    threads: 2
    log: f"{MOD1_OUT}/logs/01_rawQC/{{sample}}.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/01_rawQC
        fastqc -t {threads} -o {MOD1_OUT}/01_rawQC {input.r1} {input.r2} &> {log}
        """

rule multiqc_raw:
    input:
        expand(f"{MOD1_OUT}/01_rawQC/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
    output:
        f"{MOD1_OUT}/01_rawQC/multiqc_report.html"
    log: f"{MOD1_OUT}/logs/01_rawQC/multiqc.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        "multiqc {MOD1_OUT}/01_rawQC -o {MOD1_OUT}/01_rawQC --force &> {log}"

# ─────────────────────────────────────────────────────────────────────────────
# [2] TRIM GALORE — Q15, min 50 bp
# FIX: tmpdir + mv pattern avoids cross-version output naming bugs
# ─────────────────────────────────────────────────────────────────────────────
rule trim_galore:
    input:
        r1   = get_r1,
        r2   = get_r2,
        flag = f"{MOD1_OUT}/00_validation/{{sample}}.ok",
    output:
        r1     = f"{MOD1_OUT}/02_trimmed/{{sample}}_R1_trimmed.fq.gz",
        r2     = f"{MOD1_OUT}/02_trimmed/{{sample}}_R2_trimmed.fq.gz",
        report = f"{MOD1_OUT}/02_trimmed/{{sample}}_trimming_report.txt",
    params:
        tmpdir  = f"{MOD1_OUT}/02_trimmed/tmp_{{sample}}",
        quality = config["mod1_QC"]["min_quality"],
        length  = config["mod1_QC"]["min_length"],
    threads: THREADS
    log: f"{MOD1_OUT}/logs/02_trimmed/{{sample}}.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {params.tmpdir} {MOD1_OUT}/02_trimmed

        trim_galore --paired \
            --quality {params.quality} \
            --length {params.length} \
            --cores {threads} \
            --gzip \
            --output_dir {params.tmpdir} \
            {input.r1} {input.r2} &> {log}

        mv $(ls {params.tmpdir}/*_val_1.fq.gz) {output.r1}
        mv $(ls {params.tmpdir}/*_val_2.fq.gz) {output.r2}
        cat {params.tmpdir}/*_trimming_report.txt > {output.report}
        rm -rf {params.tmpdir}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [3] PhiX REMOVAL — BBDuk k=31
# FIX: k=31 recommended for PhiX; ftm=5 for 151bp reads
# ─────────────────────────────────────────────────────────────────────────────
rule remove_phix:
    input:
        r1 = f"{MOD1_OUT}/02_trimmed/{{sample}}_R1_trimmed.fq.gz",
        r2 = f"{MOD1_OUT}/02_trimmed/{{sample}}_R2_trimmed.fq.gz",
    output:
        r1    = temp(f"{MOD1_OUT}/03_nophix/{{sample}}_R1.fq.gz"),
        r2    = temp(f"{MOD1_OUT}/03_nophix/{{sample}}_R2.fq.gz"),
        stats = f"{MOD1_OUT}/03_nophix/{{sample}}_phix_stats.txt",
    threads: THREADS
    log: f"{MOD1_OUT}/logs/03_nophix/{{sample}}.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/03_nophix

        bbduk.sh \
            in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            ref={PHIX_REF} \
            k=31 hdist=1 ftm=5 \
            minlength=50 \
            stats={output.stats} \
            threads={threads} \
            overwrite=true \
            &> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [3b] LOW-COMPLEXITY FILTER — BBDuk entropy
# Essential for soil samples — removes poly-N and homo-polymer reads
# ─────────────────────────────────────────────────────────────────────────────
rule remove_lowcomplexity:
    input:
        r1 = f"{MOD1_OUT}/03_nophix/{{sample}}_R1.fq.gz",
        r2 = f"{MOD1_OUT}/03_nophix/{{sample}}_R2.fq.gz",
    output:
        r1    = temp(f"{MOD1_OUT}/03b_lowcomp/{{sample}}_R1.fq.gz"),
        r2    = temp(f"{MOD1_OUT}/03b_lowcomp/{{sample}}_R2.fq.gz"),
        stats = f"{MOD1_OUT}/03b_lowcomp/{{sample}}_lowcomp_stats.txt",
    params:
        entropy = config["mod1_QC"]["entropy"],
    threads: THREADS
    log: f"{MOD1_OUT}/logs/03b_lowcomp/{{sample}}.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/03b_lowcomp

        bbduk.sh \
            in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            entropy={params.entropy} \
            entropywindow=50 entropymask=f \
            minlength=50 \
            stats={output.stats} \
            threads={threads} \
            overwrite=true \
            &> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [4] HOST REMOVAL — Bowtie2 + samtools
# FIX: samtools extraction — not --un-conc-gz % (unreliable across versions)
# ─────────────────────────────────────────────────────────────────────────────
rule remove_host:
    input:
        r1 = f"{MOD1_OUT}/03b_lowcomp/{{sample}}_R1.fq.gz",
        r2 = f"{MOD1_OUT}/03b_lowcomp/{{sample}}_R2.fq.gz",
    output:
        r1      = f"{MOD1_OUT}/04_nohost/{{sample}}_R1.fq.gz",
        r2      = f"{MOD1_OUT}/04_nohost/{{sample}}_R2.fq.gz",
        bt2_log = f"{MOD1_OUT}/04_nohost/{{sample}}_bowtie2.log",
    threads: THREADS
    log: f"{MOD1_OUT}/logs/04_nohost/{{sample}}.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/04_nohost

        bowtie2 \
            -x {HOST_IDX} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            --sensitive \
            2> {output.bt2_log} \
        | samtools view -b -f 12 -F 256 \
        | samtools sort -n -@ {threads} \
        | samtools fastq \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null \
            -s /dev/null \
            -n \
            &> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [5] FastQC + MultiQC — POST-CLEAN
# ─────────────────────────────────────────────────────────────────────────────
rule fastqc_clean:
    input:
        r1 = f"{MOD1_OUT}/04_nohost/{{sample}}_R1.fq.gz",
        r2 = f"{MOD1_OUT}/04_nohost/{{sample}}_R2.fq.gz",
    output:
        html_r1 = f"{MOD1_OUT}/05_cleanQC/{{sample}}_R1_fastqc.html",
        html_r2 = f"{MOD1_OUT}/05_cleanQC/{{sample}}_R2_fastqc.html",
        zip_r1  = f"{MOD1_OUT}/05_cleanQC/{{sample}}_R1_fastqc.zip",
        zip_r2  = f"{MOD1_OUT}/05_cleanQC/{{sample}}_R2_fastqc.zip",
    threads: 2
    log: f"{MOD1_OUT}/logs/05_cleanQC/{{sample}}.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/05_cleanQC
        fastqc -t {threads} -o {MOD1_OUT}/05_cleanQC {input.r1} {input.r2} &> {log}
        """

rule multiqc_clean:
    input:
        expand(f"{MOD1_OUT}/05_cleanQC/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
    output:
        f"{MOD1_OUT}/05_cleanQC/multiqc_report.html"
    log: f"{MOD1_OUT}/logs/05_cleanQC/multiqc.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        "multiqc {MOD1_OUT}/05_cleanQC -o {MOD1_OUT}/05_cleanQC --force &> {log}"

# ─────────────────────────────────────────────────────────────────────────────
# [6] EukDetect — fungi screen (~1% expected in soil)
# ─────────────────────────────────────────────────────────────────────────────
rule eukdetect:
    input:
        r1 = f"{MOD1_OUT}/04_nohost/{{sample}}_R1.fq.gz",
        r2 = f"{MOD1_OUT}/04_nohost/{{sample}}_R2.fq.gz",
    output:
        f"{MOD1_OUT}/06_eukdetect/{{sample}}/{{sample}}_eukdetect_results.txt"
    params:
        outdir = f"{MOD1_OUT}/06_eukdetect/{{sample}}"
    threads: THREADS
    log: f"{MOD1_OUT}/logs/06_eukdetect/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        echo "eukdetect skipped" > {output}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [7] Nonpareil3 — sequencing coverage (R1 only)
# FIX: sample-specific tmpdir — not /tmp (avoids parallel job collisions)
# ─────────────────────────────────────────────────────────────────────────────
rule nonpareil:
    input:
        r1 = f"{MOD1_OUT}/04_nohost/{{sample}}_R1.fq.gz",
    output:
        npo = f"{MOD1_OUT}/07_nonpareil/{{sample}}.npo",
        npa = f"{MOD1_OUT}/07_nonpareil/{{sample}}.npa",
    params:
        prefix = f"{MOD1_OUT}/07_nonpareil/{{sample}}",
        tmp_fq = f"{MOD1_OUT}/07_nonpareil/{{sample}}_tmp_R1.fastq",
    threads: THREADS
    log: f"{MOD1_OUT}/logs/07_nonpareil/{{sample}}.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/07_nonpareil

        zcat {input.r1} > {params.tmp_fq}

        nonpareil \
            -s {params.tmp_fq} \
            -T kmer -f fastq \
            -b {params.prefix} \
            -t {threads} \
            &> {log}

        rm -f {params.tmp_fq}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [8] seqkit stats — read count, GC%, length at raw/trimmed/clean stages
# FIX: track all stages; parse by filename not by iloc position
# ─────────────────────────────────────────────────────────────────────────────
rule seqkit_stats:
    input:
        raw_r1   = get_r1,
        trim_r1  = f"{MOD1_OUT}/02_trimmed/{{sample}}_R1_trimmed.fq.gz",
        clean_r1 = f"{MOD1_OUT}/04_nohost/{{sample}}_R1.fq.gz",
        clean_r2 = f"{MOD1_OUT}/04_nohost/{{sample}}_R2.fq.gz",
    output:
        f"{MOD1_OUT}/08_seqkit/{{sample}}_seqkit.tsv"
    log: f"{MOD1_OUT}/logs/08_seqkit/{{sample}}.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/08_seqkit

        seqkit stats -a -T \
            {input.raw_r1} \
            {input.trim_r1} \
            {input.clean_r1} \
            {input.clean_r2} \
            > {output} 2> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [9] RETENTION SUMMARY
# FIX: parse by filename not iloc; per-step counts; robust regex
# ─────────────────────────────────────────────────────────────────────────────
rule retention_summary:
    input:
        seqkit     = f"{MOD1_OUT}/08_seqkit/{{sample}}_seqkit.tsv",
        phix_stats = f"{MOD1_OUT}/03_nophix/{{sample}}_phix_stats.txt",
        lowcomp    = f"{MOD1_OUT}/03b_lowcomp/{{sample}}_lowcomp_stats.txt",
        bt2_log    = f"{MOD1_OUT}/04_nohost/{{sample}}_bowtie2.log",
        npo        = f"{MOD1_OUT}/07_nonpareil/{{sample}}.npo",
        eukdetect  = f"{MOD1_OUT}/06_eukdetect/{{sample}}/{{sample}}_eukdetect_results.txt",
    output:
        f"{MOD1_OUT}/09_retention/{{sample}}_retention.tsv"
    log: f"{MOD1_OUT}/logs/09_retention/{{sample}}.log"
    run:
        import re, os, pandas as pd
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        sk = pd.read_csv(input.seqkit, sep="\t")

        def get_reads(df, fragment):
            mask = df["file"].str.contains(fragment, regex=False)
            rows = df[mask]
            return int(rows.iloc[0]["num_seqs"]) if len(rows) > 0 else None

        raw_reads   = get_reads(sk, wildcards.sample + "_R1.fastq") or sk.iloc[0]["num_seqs"]
        trim_reads  = get_reads(sk, "trimmed")
        clean_reads = get_reads(sk, f"04_nohost/{wildcards.sample}_R1")

        phix_pct = 0.0
        with open(input.phix_stats) as fh:
            for line in fh:
                m = re.search(r"Total Removed.*?([\d.]+)%", line)
                if m:
                    phix_pct = float(m.group(1)); break

        lowcomp_pct = 0.0
        with open(input.lowcomp) as fh:
            for line in fh:
                m = re.search(r"Total Removed.*?([\d.]+)%", line)
                if m:
                    lowcomp_pct = float(m.group(1)); break

        host_pct = 0.0
        with open(input.bt2_log) as fh:
            for line in fh:
                m = re.search(r"([\d.]+)% overall alignment rate", line)
                if m:
                    host_pct = float(m.group(1))

        np_redundancy = "NA"
        with open(input.npo) as fh:
            for line in fh:
                if not line.startswith("#") and line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        try:
                            np_redundancy = f"{float(parts[1]):.4f}"
                        except ValueError:
                            pass
                        break

        pct_retained = (clean_reads / raw_reads * 100) if raw_reads else 0.0

        with open(output[0], "w") as out:
            out.write("sample\traw_reads\ttrimmed_reads\tclean_reads\t"
                      "pct_retained\tphix_pct\tlowcomp_pct\thost_pct\t"
                      "nonpareil_redundancy\n")
            out.write(f"{wildcards.sample}\t{raw_reads}\t{trim_reads}\t{clean_reads}\t"
                      f"{pct_retained:.2f}\t{phix_pct:.2f}\t{lowcomp_pct:.2f}\t"
                      f"{host_pct:.2f}\t{np_redundancy}\n")

# ─────────────────────────────────────────────────────────────────────────────
# [99] MERGED MultiQC — all steps combined
# ─────────────────────────────────────────────────────────────────────────────
rule multiqc_final:
    input:
        expand(f"{MOD1_OUT}/01_rawQC/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        expand(f"{MOD1_OUT}/05_cleanQC/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        expand(f"{MOD1_OUT}/02_trimmed/{{sample}}_trimming_report.txt", sample=SAMPLES),
        expand(f"{MOD1_OUT}/03_nophix/{{sample}}_phix_stats.txt", sample=SAMPLES),
        expand(f"{MOD1_OUT}/04_nohost/{{sample}}_bowtie2.log", sample=SAMPLES),
    output:
        f"{MOD1_OUT}/99_multiqc_final/multiqc_report.html"
    log: f"{MOD1_OUT}/logs/99_multiqc_final/multiqc.log"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        set -euo pipefail
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/99_multiqc_final

        multiqc \
            {MOD1_OUT}/01_rawQC \
            {MOD1_OUT}/02_trimmed \
            {MOD1_OUT}/03_nophix \
            {MOD1_OUT}/04_nohost \
            {MOD1_OUT}/05_cleanQC \
            -o {MOD1_OUT}/99_multiqc_final \
            --force &> {log}
        """

# ─────────────────────────────────────────────────────────────────────────────
# [100] QC PLOTS
# ─────────────────────────────────────────────────────────────────────────────
rule mod1_plots:
    input:
        retention = expand(f"{MOD1_OUT}/09_retention/{{sample}}_retention.tsv", sample=SAMPLES),
        seqkit    = expand(f"{MOD1_OUT}/08_seqkit/{{sample}}_seqkit.tsv", sample=SAMPLES),
        npo       = expand(f"{MOD1_OUT}/07_nonpareil/{{sample}}.npo", sample=SAMPLES),
    output:
        f"{MOD1_OUT}/plots/05_QC_dashboard.png"
    conda: "../envs/mod1_QC.yaml"
    shell:
        """
        {JAVA_LIB_FIX}
        mkdir -p {MOD1_OUT}/plots
        python {PLOT_SCRIPT} \
            --results {MOD1_OUT} \
            --outdir {MOD1_OUT}/plots
        """
