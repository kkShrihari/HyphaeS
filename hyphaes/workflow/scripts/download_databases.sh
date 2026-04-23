#!/bin/bash
# =============================================================================
# HyphaeS — DATABASE SETUP (all 5 modules)
# Usage:
#   bash download_databases.sh
#   bash download_databases.sh /opt/databases
# =============================================================================

set -euo pipefail

BASE_DB="${1:-${HOME}/databases}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_CONFIG="${SCRIPT_DIR}/../../templates/config.yaml"
LOG_FILE="${BASE_DB}/db_setup.log"

mkdir -p "${BASE_DB}"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================================"
echo "  HyphaeS Database Setup — All Modules"
echo "  Installing to : ${BASE_DB}"
echo "  Started       : $(date)"
echo "============================================================"

DONE_DIR="${BASE_DB}/.done"
mkdir -p "${DONE_DIR}"
is_done()  { [ -f "${DONE_DIR}/$1.done" ]; }
mark_done(){ touch "${DONE_DIR}/$1.done"; }

# =============================================================================
# [1/8] Human GRCh38 — Bowtie2 pre-built index (~4 GB) — Mod1
# =============================================================================
echo ""; echo "[1/8] Human GRCh38 Bowtie2 index..."
HOST_DIR="${BASE_DB}/host"
HOST_IDX="${HOST_DIR}/GRCh38_noalt_as/GRCh38_noalt_as"
mkdir -p "${HOST_DIR}"

if is_done "grch38"; then echo "  Skipping."; else
    wget -c --progress=bar https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip \
        -O "${HOST_DIR}/GRCh38_noalt_as.zip"
    unzip "${HOST_DIR}/GRCh38_noalt_as.zip" -d "${HOST_DIR}/"
    rm "${HOST_DIR}/GRCh38_noalt_as.zip"
    mark_done "grch38"; echo "  Done."
fi

# =============================================================================
# [2/8] EukDetect DB v2 (~5 GB) — Mod1
# =============================================================================
echo ""; echo "[2/8] EukDetect database v2..."
EUK_DIR="${BASE_DB}/eukdetect"
EUK_DB_DIR="${EUK_DIR}/eukdetect_database_v2"
EUK_FNA="${EUK_DB_DIR}/ncbi_eukprot_met_arch_markers.fna"
mkdir -p "${EUK_DIR}"

if is_done "eukdetect"; then echo "  Skipping."; else
    wget -c --progress=bar "https://ndownloader.figshare.com/files/34885596" \
        -O "${EUK_DIR}/eukdetect_database_v2.tar.gz"
    tar -xzf "${EUK_DIR}/eukdetect_database_v2.tar.gz" -C "${EUK_DIR}/"
    rm "${EUK_DIR}/eukdetect_database_v2.tar.gz"
    bwa index "${EUK_FNA}"
    python -c "from ete3 import NCBITaxa; NCBITaxa().update_taxonomy_database()"
    cp ~/.etetoolkit/taxa.sqlite "${EUK_DB_DIR}/"
    mark_done "eukdetect"; echo "  Done."
fi

EUKDETECT_INSTALL=$(python -c "import eukdetect,os; print(os.path.dirname(eukdetect.__file__))" 2>/dev/null || echo "")
EUK_CONFIG="${EUK_DIR}/eukdetect_config.yml"
cat > "${EUK_CONFIG}" << EOF
database_dir: "${EUK_DB_DIR}"
database_prefix: "ncbi_eukprot_met_arch_markers.fna"
eukdetect_dir: "${EUKDETECT_INSTALL}"
samples:
EOF

# =============================================================================
# [3/8] Kraken2 PlusPF (~100 GB) — Mod3
# =============================================================================
echo ""; echo "[3/8] Kraken2 PlusPF database (~100 GB)..."
K2_DIR="${BASE_DB}/kraken2/pluspf"
mkdir -p "${K2_DIR}"

if is_done "kraken2"; then echo "  Skipping."; else
    wget -c --progress=bar \
        "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz" \
        -O "${K2_DIR}/k2_pluspf.tar.gz"
    tar -xzf "${K2_DIR}/k2_pluspf.tar.gz" -C "${K2_DIR}/"
    rm "${K2_DIR}/k2_pluspf.tar.gz"
    bracken-build -d "${K2_DIR}" -t 16 -l 150
    mark_done "kraken2"; echo "  Done."
fi

# =============================================================================
# [4/8] Kaiju nr_euk (~50 GB) — Mod3
# =============================================================================
echo ""; echo "[4/8] Kaiju nr_euk database (~50 GB)..."
KAIJU_DIR="${BASE_DB}/kaiju"
mkdir -p "${KAIJU_DIR}"

if is_done "kaiju"; then echo "  Skipping."; else
    wget -c --progress=bar \
        "https://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2023-05-10.tgz" \
        -O "${KAIJU_DIR}/kaiju_db_nr_euk.tgz"
    tar -xzf "${KAIJU_DIR}/kaiju_db_nr_euk.tgz" -C "${KAIJU_DIR}/"
    rm "${KAIJU_DIR}/kaiju_db_nr_euk.tgz"
    mark_done "kaiju"; echo "  Done."
fi

# =============================================================================
# [5/8] EukCC DB (~15 GB) — Mod4
# =============================================================================
echo ""; echo "[5/8] EukCC database (~15 GB)..."
EUKCC_DIR="${BASE_DB}/eukcc"
mkdir -p "${EUKCC_DIR}"

if is_done "eukcc"; then echo "  Skipping."; else
    wget -c --progress=bar \
        "https://openstack.cebitec.uni-bielefeld.de/s/pQ8F6RGKbvXkfhW/download" \
        -O "${EUKCC_DIR}/eukcc_db.tar.gz"
    tar -xzf "${EUKCC_DIR}/eukcc_db.tar.gz" -C "${EUKCC_DIR}/"
    rm "${EUKCC_DIR}/eukcc_db.tar.gz"
    mark_done "eukcc"; echo "  Done."
fi

# =============================================================================
# [6/8] dbCAN DB (~1 GB) — Mod4/5
# =============================================================================
echo ""; echo "[6/8] dbCAN database..."
DBCAN_DIR="${BASE_DB}/dbcan"
mkdir -p "${DBCAN_DIR}"

if is_done "dbcan"; then echo "  Skipping."; else
    for f in CAZyDB.07262023.fa dbCAN-HMMdb-V12.txt tcdb.fa tf-1.hmm tf-2.hmm stp.hmm; do
        wget -c --progress=bar \
            "https://bcb.unl.edu/dbCAN2/download/Databases/${f}" \
            -O "${DBCAN_DIR}/${f}"
    done
    diamond makedb --in "${DBCAN_DIR}/CAZyDB.07262023.fa" -d "${DBCAN_DIR}/CAZyDB"
    hmmpress "${DBCAN_DIR}/dbCAN-HMMdb-V12.txt"
    mark_done "dbcan"; echo "  Done."
fi

# =============================================================================
# [7/8] eggNOG-mapper DB (~50 GB) — Mod4/5
# =============================================================================
echo ""; echo "[7/8] eggNOG-mapper database (~50 GB)..."
EGGNOG_DIR="${BASE_DB}/eggnog"
mkdir -p "${EGGNOG_DIR}"

if is_done "eggnog"; then echo "  Skipping."; else
    download_eggnog_data.py --data_dir "${EGGNOG_DIR}" -y
    mark_done "eggnog"; echo "  Done."
fi

# =============================================================================
# [8/8] FUNGuild DB — Mod5
# =============================================================================
echo ""; echo "[8/8] FUNGuild database..."
FUNGUILD_DIR="${BASE_DB}/funguild"
mkdir -p "${FUNGUILD_DIR}"

if is_done "funguild"; then echo "  Skipping."; else
    wget -c --progress=bar \
        "https://github.com/UMNFuN/FUNGuild/raw/master/FUNGuild.db" \
        -O "${FUNGUILD_DIR}/FUNGuild.db"
    mark_done "funguild"; echo "  Done."
fi

# =============================================================================
# AUTO-PATCH config.yaml
# =============================================================================
echo ""; echo "Auto-updating config.yaml..."

KAIJU_FMI=$(ls ${KAIJU_DIR}/*.fmi 2>/dev/null | head -1 || echo "")
AUGUSTUS_CFG=$(python -c "import subprocess; r=subprocess.run(['which','augustus'],capture_output=True,text=True); import os; p=os.path.dirname(r.stdout.strip()); print(os.path.join(os.path.dirname(p),'config'))" 2>/dev/null || echo "")

sed -i "s|host_idx:.*|host_idx: \"${HOST_IDX}\"|"                "${PIPELINE_CONFIG}"
sed -i "s|eukdetect_cfg:.*|eukdetect_cfg: \"${EUK_CONFIG}\"|"    "${PIPELINE_CONFIG}"
sed -i "s|kraken2_db:.*|kraken2_db: \"${K2_DIR}\"|"              "${PIPELINE_CONFIG}"
sed -i "s|kaiju_db:.*|kaiju_db: \"${KAIJU_FMI}\"|"               "${PIPELINE_CONFIG}"
sed -i "s|kaiju_nodes:.*|kaiju_nodes: \"${KAIJU_DIR}/nodes.dmp\"|" "${PIPELINE_CONFIG}"
sed -i "s|eukcc_db:.*|eukcc_db: \"${EUKCC_DIR}\"|"               "${PIPELINE_CONFIG}"
sed -i "s|augustus_config:.*|augustus_config: \"${AUGUSTUS_CFG}\"|" "${PIPELINE_CONFIG}"
sed -i "s|dbcan_db:.*|dbcan_db: \"${DBCAN_DIR}\"|"               "${PIPELINE_CONFIG}"
sed -i "s|eggnog_db:.*|eggnog_db: \"${EGGNOG_DIR}\"|"            "${PIPELINE_CONFIG}"
sed -i "s|funguild_db:.*|funguild_db: \"${FUNGUILD_DIR}/FUNGuild.db\"|" "${PIPELINE_CONFIG}"

echo "  Done."

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo "  ALL DONE — $(date)"
echo "============================================================"
echo "  Database              | Size    | Used by"
echo "  ----------------------|---------|--------"
echo "  GRCh38 Bowtie2 index  | ~4 GB   | mod1"
echo "  EukDetect DB v2       | ~5 GB   | mod1"
echo "  Kraken2 PlusPF        | ~100 GB | mod3"
echo "  Kaiju nr_euk          | ~50 GB  | mod3, mod4"
echo "  EukCC DB              | ~15 GB  | mod4"
echo "  dbCAN DB              | ~1 GB   | mod4, mod5"
echo "  eggNOG-mapper DB      | ~50 GB  | mod4, mod5"
echo "  FUNGuild DB           | <1 MB   | mod5"
echo "  ----------------------|---------|--------"
echo "  Total                 | ~225 GB |"
echo ""
echo "  Log: ${LOG_FILE}"
echo "============================================================"
