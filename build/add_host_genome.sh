#!/usr/bin/env bash
# add_host_genome.sh
# ------------------
# Generic framework for extending any GTDB-derived Kraken2/Bracken database
# with a eukaryotic host genome node.
#
# Background:
#   GTDB taxonomy is prokaryote-only. Any Kraken2/Bracken database built from
#   GTDB-derived MAGs therefore lacks a host taxonomy entry. When host-associated
#   samples are analysed, residual host reads that evade upstream host depletion
#   are either left unclassified or systematically mis-assigned to microbial taxa.
#   This script resolves that structural gap by appending a single custom host
#   node under root (taxid 1) and rebuilding the database index.
#
# Prerequisites:
#   - Kraken2 (tested with v2.1.x)
#   - Bracken (bracken-build, tested with v3.0.x)
#   - Python 3
#   - taxonomy/add_host_taxon.py (in the same repository)
#
# Usage:
#   bash add_host_genome.sh \
#     --src-db       /path/to/any_gtdb_kraken2_db \
#     --dst-db       /path/to/output_db_with_host \
#     --host-fasta   /path/to/host_genome.fna \
#     --host-name    "Rattus norvegicus" \
#     --host-taxid   3735 \
#     [--threads     8] \
#     [--read-lengths "50 100 150 200 250"]
#
# Taxid selection:
#   Choose a sequential integer absent from the source database's nodes.dmp.
#   Check the highest existing taxid with:
#     awk -F'\t' '{print $1}' /path/to/db/taxonomy/nodes.dmp | sort -n | tail -5
#   Use a value at least 1 higher than the current maximum (or any unused value).
#   Do NOT reuse the NCBI taxid for the same organism, as this can cause
#   lineage conflicts in downstream NCBI-based tools.
#
# Examples:
#   Mouse   : --host-name "Mus musculus"      --host-taxid 3734
#   Rat     : --host-name "Rattus norvegicus"  --host-taxid 3735
#   Hamster : --host-name "Mesocricetus auratus" --host-taxid 3736
#
# Reference:
#   Chen J-Y et al. (2026). Extending the MGnify mouse-gut Kraken2 database
#   with host taxonomy greatly reduces systematic false-positive classifications.
#   Microbial Genomics (Technical Resource), in preparation.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAXONOMY_SCRIPT="${SCRIPT_DIR}/../taxonomy/add_host_taxon.py"

# ---------- defaults ----------
SRC_DB=""
DST_DB=""
HOST_FASTA=""
HOST_NAME=""
HOST_TAXID=""
THREADS=8
READ_LENGTHS="50 100 150 200 250"

usage() {
  cat <<EOF
Usage: $0 \\
  --src-db       <path>   Source GTDB-derived Kraken2 database directory
  --dst-db       <path>   Destination directory for the extended database
  --host-fasta   <fasta>  Host genome FASTA (uncompressed or .gz)
  --host-name    <str>    Scientific name (e.g. "Rattus norvegicus")
  --host-taxid   <int>    Custom taxid not already present in nodes.dmp
  [--threads     <int>]   CPU threads (default: 8)
  [--read-lengths <str>]  Bracken lengths, space-separated (default: "50 100 150 200 250")
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --src-db)       SRC_DB="$2";        shift 2 ;;
    --dst-db)       DST_DB="$2";        shift 2 ;;
    --host-fasta)   HOST_FASTA="$2";    shift 2 ;;
    --host-name)    HOST_NAME="$2";     shift 2 ;;
    --host-taxid)   HOST_TAXID="$2";    shift 2 ;;
    --threads)      THREADS="$2";       shift 2 ;;
    --read-lengths) READ_LENGTHS="$2";  shift 2 ;;
    -h|--help)      usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

[[ -z "$SRC_DB" || -z "$DST_DB" || -z "$HOST_FASTA" || -z "$HOST_NAME" || -z "$HOST_TAXID" ]] && {
  echo "[ERROR] Missing required arguments." >&2
  usage
}

[[ ! -d "$SRC_DB" ]]     && { echo "[ERROR] Source database not found: $SRC_DB"; exit 1; }
[[ ! -f "$HOST_FASTA" ]] && { echo "[ERROR] Host FASTA not found: $HOST_FASTA"; exit 1; }
[[ ! -f "$TAXONOMY_SCRIPT" ]] && { echo "[ERROR] Script not found: $TAXONOMY_SCRIPT"; exit 1; }

echo "============================================================="
echo " add_host_genome.sh — GTDB Kraken2 database host-node extension"
echo "============================================================="
echo " Source DB   : $SRC_DB"
echo " Target DB   : $DST_DB"
echo " Host FASTA  : $HOST_FASTA"
echo " Host name   : $HOST_NAME"
echo " Host taxid  : $HOST_TAXID"
echo " Threads     : $THREADS"
echo " Read lengths: $READ_LENGTHS"
echo "============================================================="

# ---------- Step 1: copy source ----------
echo ""
echo "[1/5] Copying source database to $DST_DB ..."
mkdir -p "$DST_DB"
cp -r "$SRC_DB/." "$DST_DB/"
# Remove compiled index files; they are regenerated in step 5
rm -f "$DST_DB/hash.k2d" "$DST_DB/taxo.k2d" "$DST_DB/opts.k2d" \
      "$DST_DB/seqid2taxid.map" "$DST_DB/database.kraken"
echo "      Done."

# ---------- Step 2: modify taxonomy ----------
echo ""
echo "[2/5] Appending host node (taxid ${HOST_TAXID}: ${HOST_NAME}) to taxonomy files..."
python3 "$TAXONOMY_SCRIPT" \
  --names        "${DST_DB}/taxonomy/names.dmp" \
  --nodes        "${DST_DB}/taxonomy/nodes.dmp" \
  --species-name "$HOST_NAME" \
  --taxid        "$HOST_TAXID"
echo "      Taxonomy modification complete."

# ---------- Step 3: prepare host FASTA ----------
echo ""
echo "[3/5] Tagging FASTA headers with kraken:taxid|${HOST_TAXID} ..."
HOST_K2="${DST_DB}/library/host/host_k2.fa"
mkdir -p "${DST_DB}/library/host"

# Support gzipped FASTA
if [[ "$HOST_FASTA" == *.gz ]]; then
  zcat "$HOST_FASTA" | sed "s/^>\(\S*\)/>\1|kraken:taxid|${HOST_TAXID}/" > "$HOST_K2"
else
  sed "s/^>\(\S*\)/>\1|kraken:taxid|${HOST_TAXID}/" "$HOST_FASTA" > "$HOST_K2"
fi
echo "      Tagged FASTA: $HOST_K2"

# ---------- Step 4: add to Kraken2 library ----------
echo ""
echo "[4/5] Adding host genome to Kraken2 library..."
kraken2-build --add-to-library "$HOST_K2" --db "$DST_DB"
echo "      Done."

# ---------- Step 5a: build Kraken2 database ----------
echo ""
echo "[5a/5] Building Kraken2 database..."
kraken2-build --build --db "$DST_DB" --threads "$THREADS"
echo "       Kraken2 build complete."

# ---------- Step 5b: build Bracken ----------
echo ""
echo "[5b/5] Building Bracken k-mer distributions (lengths: ${READ_LENGTHS})..."
for len in $READ_LENGTHS; do
  echo "       bracken-build -l ${len}..."
  bracken-build -d "$DST_DB" -k 35 -l "$len" -t "$THREADS"
done
echo "       Bracken build complete."

# ---------- Verification ----------
echo ""
echo "============================================================="
echo " Build complete. Verifying host k-mer coverage..."
echo "============================================================="
kraken2-inspect --db "$DST_DB" --report-zero-counts \
  | grep -i "${HOST_NAME}\|taxid ${HOST_TAXID}\|^${HOST_TAXID}\b" \
  | head -5

echo ""
echo "Host node in taxonomy:"
echo "  names.dmp : $(grep "^${HOST_TAXID}" "${DST_DB}/taxonomy/names.dmp")"
echo "  nodes.dmp : $(grep "^${HOST_TAXID}" "${DST_DB}/taxonomy/nodes.dmp")"
echo ""
echo "Full inspect: kraken2-inspect --db $DST_DB --report-zero-counts"
