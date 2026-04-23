"""
HyphaeS CLI — Hyphal Metagenomics Suite
"""

import os
import sys
import glob
import shutil
import subprocess
from pathlib import Path
from datetime import datetime

import click
import yaml
import pandas as pd

PKG_DIR      = Path(__file__).parent
TEMPLATES    = PKG_DIR / "templates"
SNAKEFILE    = PKG_DIR / "workflow" / "Snakefile"
SETUP_SCRIPT = PKG_DIR / "workflow" / "scripts" / "download_databases.sh"

MODULES = ["mod1", "mod2", "mod3", "mod4", "mod5"]

MODULE_NAMES = {
    "mod1": "Quality Control",
    "mod2": "Assembly",
    "mod3": "Taxonomy",
    "mod4": "Binning",
    "mod5": "Functional Annotation",
}

MODULE_REQUIRED_INPUTS = {
    "mod1": [],
    "mod2": [
        "results/mod1_QC/04_nohost/{sample}_R1.fq.gz",
        "results/mod1_QC/04_nohost/{sample}_R2.fq.gz",
        "results/mod1_QC/09_retention/{sample}_retention.tsv",
    ],
    "mod3": ["results/mod2_assembly/{sample}/contigs.fasta"],
    "mod4": ["results/mod3_taxonomy/{sample}/taxonomy.tsv"],
    "mod5": ["results/mod4_binning/{sample}/dastool_bins/"],
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def load_config():
    cfg = Path("config.yaml")
    if not cfg.exists():
        click.echo("ERROR: config.yaml not found. Run 'HyphaeS init' first.", err=True)
        sys.exit(1)
    with open(cfg) as f:
        return yaml.safe_load(f)

def load_samples():
    tsv = Path("samples.tsv")
    if not tsv.exists():
        click.echo("ERROR: samples.tsv not found. Run 'HyphaeS init' first.", err=True)
        sys.exit(1)
    return pd.read_csv(tsv, sep="\t")

def get_log_path(label):
    log_dir = Path("logs")
    log_dir.mkdir(exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    return log_dir / f"{label}_{ts}.log"

def stream_command(cmd, log_path):
    with open(log_path, "w") as lf:
        lf.write(f"Command : {' '.join(str(c) for c in cmd)}\n")
        lf.write(f"Started : {datetime.now()}\n")
        lf.write("=" * 60 + "\n\n")
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        for line in process.stdout:
            sys.stdout.write(line)
            lf.write(line)
        process.wait()
        lf.write(f"\n\nFinished  : {datetime.now()}\n")
        lf.write(f"Exit code : {process.returncode}\n")
    return process.returncode

def get_frontend():
    """Use mamba if available, else conda."""
    r = subprocess.run(["which", "mamba"], capture_output=True, text=True)
    return "mamba" if r.returncode == 0 else "conda"

def get_conda_envs_dir():
    frontend = get_frontend()
    r = subprocess.run([frontend, "info", "--base"], capture_output=True, text=True)
    if r.returncode == 0:
        return str(Path(r.stdout.strip()) / "envs")
    return None

def get_snakemake():
    """Find the correct snakemake — same Python env as HyphaeS."""
    import sys
    from pathlib import Path
    # Use snakemake from same Python as HyphaeS
    snake = Path(sys.executable).parent / "snakemake"
    if snake.exists():
        return str(snake)
    return "snakemake"

def build_snakemake_cmd(module, cores, conda_prefix, conda_frontend, extra_flags):
    cmd = [
        get_snakemake(),
        "--snakefile",      str(SNAKEFILE),
        "--configfile",     "config.yaml",
        "--config",         f"active_module={module}",
        "--cores",          str(cores),
        "--use-conda",
        "--conda-frontend", conda_frontend,
        "--rerun-incomplete",
        "--latency-wait",   "60",
        "--keep-going",
        "--printshellcmds",
    ]
    if conda_prefix:
        cmd += ["--conda-prefix", conda_prefix]
    cmd += extra_flags
    return cmd

# ── CLI GROUP ─────────────────────────────────────────────────────────────────
@click.group()
@click.version_option("0.1.0", prog_name="HyphaeS")
def main():
    """
    HyphaeS — Hyphal Metagenomics Suite

    Shotgun soil metagenomics for low-abundance fungi (~1%).
    Paired-end short reads · Snakemake · 5 modules.

    \b
    First time setup:
      HyphaeS init --reads /path/to/fastq/
      HyphaeS setup
      HyphaeS install-envs

    \b
    Run analysis:
      HyphaeS check mod1
      HyphaeS run mod1 --dry
      HyphaeS run mod1
    """
    pass

# ─────────────────────────────────────────────────────────────────────────────
# INIT
# ─────────────────────────────────────────────────────────────────────────────
@main.command()
@click.option("--force", is_flag=True,
              help="Overwrite existing config.yaml and samples.tsv")
@click.option("--reads", default=None, metavar="DIR",
              help="Folder containing fastq.gz files — auto-generates samples.tsv")
def init(force, reads):
    """Initialize a new HyphaeS project in the current folder."""
    click.echo("\nHyphaeS — Initializing project...\n")

    cfg_dest = Path("config.yaml")
    if cfg_dest.exists() and not force:
        click.echo("  SKIP : config.yaml already exists (use --force to overwrite)")
    else:
        shutil.copy(TEMPLATES / "config.yaml", cfg_dest)
        click.echo("  Created : config.yaml")

    tsv_dest = Path("samples.tsv")

    if reads:
        reads_path = Path(reads)
        if not reads_path.exists():
            click.echo(f"  ERROR: --reads folder not found: {reads}", err=True)
            sys.exit(1)

        r1_files = []
        for pat in ["*_R1.fastq.gz", "*_R1.fq.gz", "*_1.fastq.gz", "*_1.fq.gz"]:
            r1_files.extend(glob.glob(str(reads_path / pat)))
        r1_files = sorted(set(r1_files))

        if not r1_files:
            click.echo(f"  ERROR: No R1 fastq.gz files found in: {reads}", err=True)
            sys.exit(1)

        rows = []
        for r1 in r1_files:
            r1p = Path(r1)
            r2p = sample = None
            for r1_suf, r2_suf in [
                ("_R1.fastq.gz", "_R2.fastq.gz"),
                ("_R1.fq.gz",    "_R2.fq.gz"),
                ("_1.fastq.gz",  "_2.fastq.gz"),
                ("_1.fq.gz",     "_2.fq.gz"),
            ]:
                if r1p.name.endswith(r1_suf):
                    r2p    = r1p.parent / r1p.name.replace(r1_suf, r2_suf)
                    sample = r1p.name.replace(r1_suf, "")
                    break
            if r2p and r2p.exists():
                rows.append({"sample": sample, "r1": str(r1p.resolve()), "r2": str(r2p.resolve())})

        if not rows:
            click.echo("  ERROR: No valid R1/R2 pairs found.", err=True)
            sys.exit(1)

        if tsv_dest.exists() and not force:
            click.echo("  SKIP : samples.tsv already exists (use --force to overwrite)")
        else:
            with open(tsv_dest, "w") as f:
                f.write("sample\tr1\tr2\n")
                for row in rows:
                    f.write(f"{row['sample']}\t{row['r1']}\t{row['r2']}\n")
            click.echo(f"  Created : samples.tsv — {len(rows)} samples detected:")
            for row in rows:
                click.echo(f"    {row['sample']}")
    else:
        if tsv_dest.exists() and not force:
            click.echo("  SKIP : samples.tsv already exists (use --force to overwrite)")
        else:
            shutil.copy(TEMPLATES / "samples.tsv", tsv_dest)
            click.echo("  Created : samples.tsv (template — edit with your paths)")

    click.echo("\nNext steps:")
    click.echo("  1. HyphaeS setup        — download databases")
    click.echo("  2. HyphaeS install-envs — install conda environments")
    click.echo("  3. HyphaeS run mod1 --dry\n")

# ─────────────────────────────────────────────────────────────────────────────
# SETUP
# ─────────────────────────────────────────────────────────────────────────────
@main.command()
@click.option("--db-dir",         default=None)
@click.option("--skip-grch38",    is_flag=True)
@click.option("--skip-eukdetect", is_flag=True)
def setup(db_dir, skip_grch38, skip_eukdetect):
    """Download and configure all required databases."""
    db_dir   = db_dir or str(Path("databases").resolve())
    log_path = get_log_path("db_setup")
    click.echo(f"\nHyphaeS — Setting up databases...\n  To: {db_dir}\n")
    env = os.environ.copy()
    if skip_grch38:    env["SKIP_GRCH38"]    = "1"
    if skip_eukdetect: env["SKIP_EUKDETECT"] = "1"
    rc = stream_command(["bash", str(SETUP_SCRIPT), db_dir], log_path)
    if rc != 0:
        click.echo(f"\nERROR: Setup failed. Check {log_path}", err=True)
        sys.exit(rc)
    click.echo(f"\nSetup complete.\n")

# ─────────────────────────────────────────────────────────────────────────────
# INSTALL-ENVS
# Uses snakemake --conda-create-envs-only to install all conda environments
# defined in workflow/envs/*.yaml — same approach as ATLAS
# ─────────────────────────────────────────────────────────────────────────────
@main.command("install-envs")
@click.option("--module",
              type=click.Choice(MODULES + ["all"]),
              default="all", show_default=True,
              help="Which module environment to install")
@click.option("--conda-prefix", default=None,
              help="Where to store conda environments")
@click.option("--conda-frontend",
              type=click.Choice(["conda", "mamba"]),
              default="mamba", show_default=True,
              help="Conda frontend to use")
def install_envs(module, conda_prefix, conda_frontend):
    """Install conda environments for all pipeline tools.

    \b
    Same approach as ATLAS — snakemake reads the yaml files
    in workflow/envs/ and creates environments automatically.

    \b
    Examples:
      HyphaeS install-envs
      HyphaeS install-envs --module mod1
    """
    # Auto-detect frontend
    if conda_frontend == "auto":
        conda_frontend = get_frontend()

    log_path = get_log_path("install_envs")
    click.echo(f"\nHyphaeS — Installing conda environments...")
    click.echo(f"  Module         : {module}")
    click.echo(f"  Conda frontend : {conda_frontend}")
    click.echo(f"  Log            : {log_path}\n")

    # Fix mamba prefix bug — create envs dir and pre-create prefix folders
    Path.home().joinpath(".local/share/mamba/envs").mkdir(parents=True, exist_ok=True)

    # Clean any empty stale conda env folders
    conda_dir = Path(".snakemake/conda")
    if conda_dir.exists():
        for item in conda_dir.iterdir():
            if item.is_dir() and not any(item.iterdir()):
                shutil.rmtree(item)
        # Pre-create prefix folder for each yaml to avoid mamba bug
        # for yaml_file in conda_dir.glob("*.yaml"):
        #     prefix = yaml_file.with_suffix("")
        #     prefix.mkdir(exist_ok=True)

    cmd = [
        get_snakemake(),
        "--snakefile",             str(SNAKEFILE),
        "--configfile",            "config.yaml",
        "--cores",                 "1",
        "--use-conda",
        "--conda-create-envs-only",
        "--conda-frontend",        conda_frontend,
    ]
    if conda_prefix:
        cmd += ["--conda-prefix", conda_prefix]
    if module != "all":
        cmd += ["--config", f"active_module={module}"]

    rc = stream_command(cmd, log_path)
    if rc != 0:
        click.echo(f"\nERROR: Environment install failed. Check {log_path}", err=True)
        sys.exit(rc)
    click.echo(f"\nAll environments installed. Log: {log_path}\n")

# ─────────────────────────────────────────────────────────────────────────────
# CHECK
# ─────────────────────────────────────────────────────────────────────────────
@main.command()
@click.argument("module", type=click.Choice(MODULES))
@click.option("--verbose", "-v", is_flag=True)
def check(module, verbose):
    """Validate inputs before running a module."""
    click.echo(f"\nHyphaeS — Checking {module}...\n")
    failed = False

    if not Path("config.yaml").exists():
        click.echo("  FAIL : config.yaml missing")
        sys.exit(1)
    click.echo("  OK   : config.yaml found")

    df       = load_samples()
    col_fail = False
    for col in ["sample", "r1", "r2"]:
        if col not in df.columns:
            click.echo(f"  FAIL : samples.tsv missing column: '{col}'")
            col_fail = True
            failed   = True
    if not col_fail:
        click.echo(f"  OK   : samples.tsv valid ({len(df)} samples)")

    missing_files = []
    for _, row in df.iterrows():
        for col in ["r1", "r2"]:
            if not Path(str(row[col])).exists():
                missing_files.append(f"    {row['sample']} {col.upper()}: {row[col]}")
    if missing_files:
        click.echo(f"  FAIL : Missing fastq files:")
        for m in missing_files:
            click.echo(m)
        failed = True
    else:
        click.echo("  OK   : All fastq paths exist")

    required = MODULE_REQUIRED_INPUTS.get(module, [])
    if required:
        samples      = df["sample"].tolist()
        prev_mod     = f"mod{int(module[-1]) - 1}"
        missing_outs = []
        for pattern in required:
            for s in samples:
                p = Path(pattern.format(sample=s))
                if not p.exists():
                    missing_outs.append(f"    {p}")
        if missing_outs:
            click.echo(f"  FAIL : Missing {prev_mod} outputs:")
            for m in missing_outs:
                click.echo(m)
            failed = True
        else:
            click.echo(f"  OK   : {prev_mod} outputs present")

    if failed:
        click.echo(f"\nCheck FAILED.\n")
        sys.exit(1)
    click.echo(f"\nAll checks passed — ready to run {module}.\n")

# ─────────────────────────────────────────────────────────────────────────────
# RUN
# ─────────────────────────────────────────────────────────────────────────────
@main.command()
@click.argument("module", type=click.Choice(MODULES + ["all"]))
@click.option("--dry",    "-n", is_flag=True)
@click.option("--force",  "-f", is_flag=True)
@click.option("--cores",  "-c", default=None, type=int)
@click.option("--until",        default=None, metavar="RULE")
@click.option("--touch",        is_flag=True)
@click.option("--reason",       is_flag=True)
@click.option("--verbose", "-v",is_flag=True)
@click.option("--quiet",   "-q",is_flag=True)
@click.option("--conda-prefix", default=None)
@click.option("--conda-frontend",
              type=click.Choice(["conda", "mamba"]),
              default="conda", show_default=True)
@click.option("--resources",    default=None, metavar="KEY=VAL")
@click.option("--skip-check",   is_flag=True)
def run(module, dry, force, cores, until, touch, reason,
        verbose, quiet, conda_prefix, conda_frontend, resources, skip_check):
    """Run a pipeline module or all modules sequentially."""
    if cores is None:
        cfg   = load_config()
        cores = cfg.get("general", {}).get("threads", 16)

    if not conda_prefix:
        conda_prefix = get_conda_envs_dir()

    extra = []
    if dry:       extra.append("--dry-run")
    if force:     extra.append("--forceall")
    if touch:     extra.append("--touch")
    if reason:    extra.append("--reason")
    if verbose:   extra.append("--verbose")
    if quiet:     extra += ["--quiet", "all"]
    if until:     extra += ["--until", until]
    if resources: extra += ["--resources", resources]

    run_modules = MODULES if module == "all" else [module]

    for mod in run_modules:
        click.echo(f"\n{'='*60}")
        click.echo(f"  HyphaeS — {mod}: {MODULE_NAMES[mod]}")
        click.echo(f"  Started : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        click.echo(f"  Cores   : {cores}")
        if dry:   click.echo("  Mode    : DRY RUN")
        if force: click.echo("  Mode    : FORCE RERUN ALL")
        click.echo(f"{'='*60}\n")

        if not dry and not skip_check:
            ctx = click.get_current_context()
            ctx.invoke(check, module=mod, verbose=False)

        log_path = get_log_path(f"pipeline_{mod}")
        click.echo(f"  Log: {log_path}\n")

        cmd = build_snakemake_cmd(mod, cores, conda_prefix, conda_frontend, extra)
        rc  = stream_command(cmd, log_path)

        if rc != 0:
            click.echo(f"\nERROR: {mod} failed. Check {log_path}", err=True)
            sys.exit(rc)

        click.echo(f"\n  {mod} complete : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        click.echo(f"  Log saved to  : {log_path}\n")

    click.echo("\nHyphaeS pipeline complete.\n")

# ─────────────────────────────────────────────────────────────────────────────
# INFO
# ─────────────────────────────────────────────────────────────────────────────
@main.command()
def info():
    """Show current project status."""
    click.echo("\nHyphaeS — Project Status\n")
    cfg_ok = Path("config.yaml").exists()
    smp_ok = Path("samples.tsv").exists()
    click.echo(f"  config.yaml : {'✓ found' if cfg_ok else '✗ missing'}")
    click.echo(f"  samples.tsv : {'✓ found' if smp_ok else '✗ missing'}")
    if smp_ok:
        df = pd.read_csv("samples.tsv", sep="\t")
        click.echo(f"  Samples     : {len(df)}")
    click.echo("")
    for mod, path in {
        "mod1": "results/mod1_QC/99_multiqc_final/multiqc_report.html",
        "mod2": "results/mod2_assembly",
        "mod3": "results/mod3_taxonomy",
        "mod4": "results/mod4_binning",
        "mod5": "results/mod5_annotation",
    }.items():
        done = Path(path).exists()
        click.echo(f"  {mod}  {MODULE_NAMES[mod]:<28}: {'✓ complete' if done else '○ not run'}")
    click.echo("")

# ─────────────────────────────────────────────────────────────────────────────
# CLEAN
# ─────────────────────────────────────────────────────────────────────────────
@main.command()
@click.option("--module", "modules", type=click.Choice(MODULES + ["all"]), default="all")
@click.option("--logs",  is_flag=True)
@click.option("--yes", "-y", is_flag=True)
def clean(modules, logs, yes):
    """Remove pipeline output files."""
    targets = MODULES if modules == "all" else [modules]
    dirs_to_remove = []
    for mod in targets:
        dirs_to_remove.extend(glob.glob(f"results/{mod}*"))
        if logs:
            dirs_to_remove.extend(glob.glob(f"logs/pipeline_{mod}*"))

    if not dirs_to_remove:
        click.echo("\nNothing to clean.\n")
        return

    click.echo("\nWill remove:")
    for d in dirs_to_remove:
        click.echo(f"  {d}")

    if not yes:
        click.confirm("\nProceed?", abort=True)

    for d in dirs_to_remove:
        p = Path(d)
        if p.is_dir():  shutil.rmtree(p)
        elif p.is_file(): p.unlink()
        click.echo(f"  Removed: {d}")

    click.echo("\nClean complete.\n")

if __name__ == "__main__":
    main()
