#!/usr/bin/env python3
"""
add_host_taxon.py
-----------------
Adds a custom Mus musculus taxonomy node (taxid 3734) to the MGnify mouse-gut
Kraken2 database v1.0 taxonomy files (names.dmp and nodes.dmp), producing the
corrected v1.1 taxonomy.

Usage:
    python add_host_taxon.py \
        --names  /path/to/v1.0/taxonomy/names.dmp \
        --nodes  /path/to/v1.0/taxonomy/nodes.dmp \
        --out-names  /path/to/v1.1/taxonomy/names.dmp \
        --out-nodes  /path/to/v1.1/taxonomy/nodes.dmp

If --out-names / --out-nodes are omitted, the input files are modified in-place
(a .bak backup is created automatically).

Note:
    taxid 3734 is a custom sequential identifier within the MGnify taxonomy tree.
    It is numerically distinct from NCBI Taxonomy ID 10090 (Mus musculus).
    The node is placed directly under root (taxid 1) because the GTDB-derived
    taxonomy does not include eukaryotic superkingdom nodes.

Reference:
    Chen J-Y et al. (2026). Extending the MGnify mouse-gut Kraken2 database
    with host taxonomy eliminates systematic false-positive classifications.
    GigaScience (Technical Note), in preparation [DOI: pending]
"""

import argparse
import shutil
import sys
from pathlib import Path

TAXID = "3734"
PARENT_TAXID = "1"
RANK = "species"

# Exact lines to append (tab-separated, NCBI dmp format)
NAMES_LINE = "\t|\t".join([TAXID, "Mus musculus", "", "scientific name"]) + "\t|"
NODES_LINE = "\t|\t".join([
    TAXID,        # taxid
    PARENT_TAXID, # parent taxid
    RANK,         # rank
    "",           # embl code
    "0",          # division id
    "1",          # inherited div flag
    "11",         # genetic code id
    "1",          # inherited GC flag
    "0",          # mitochondrial genetic code id
    "1",          # inherited MGC flag
    "1",          # GenBank hidden flag
    "0",          # hidden subtree root flag
]) + "\t|"


def check_not_already_present(filepath: Path, taxid: str) -> None:
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(taxid + "\t"):
                print(f"[WARNING] taxid {taxid} already present in {filepath}. Skipping append.")
                return
    return None


def append_line(src: Path, dst: Path, line: str, label: str) -> None:
    if dst != src:
        shutil.copy2(src, dst)
    # Check if already present
    with open(dst, "r") as f:
        content = f.read()
    if content.startswith(TAXID + "\t") or ("\n" + TAXID + "\t") in content:
        print(f"[WARNING] taxid {TAXID} already present in {dst.name}. No change made.")
        return
    # Ensure file ends with newline before appending
    if content and not content.endswith("\n"):
        with open(dst, "a") as f:
            f.write("\n")
    with open(dst, "a") as f:
        f.write(line + "\n")
    print(f"[OK] Appended {label} entry to {dst}")


def main():
    parser = argparse.ArgumentParser(
        description="Add Mus musculus taxid 3734 to MGnify v1.0 taxonomy files."
    )
    parser.add_argument("--names", required=True, help="Path to input names.dmp")
    parser.add_argument("--nodes", required=True, help="Path to input nodes.dmp")
    parser.add_argument("--out-names", default=None, help="Path to output names.dmp (default: in-place)")
    parser.add_argument("--out-nodes", default=None, help="Path to output nodes.dmp (default: in-place)")
    args = parser.parse_args()

    names_src = Path(args.names)
    nodes_src = Path(args.nodes)

    if not names_src.exists():
        sys.exit(f"[ERROR] names.dmp not found: {names_src}")
    if not nodes_src.exists():
        sys.exit(f"[ERROR] nodes.dmp not found: {nodes_src}")

    names_dst = Path(args.out_names) if args.out_names else names_src
    nodes_dst = Path(args.out_nodes) if args.out_nodes else nodes_src

    # Create parent directories for output if needed
    names_dst.parent.mkdir(parents=True, exist_ok=True)
    nodes_dst.parent.mkdir(parents=True, exist_ok=True)

    # Backup if modifying in-place
    if names_dst == names_src:
        shutil.copy2(names_src, str(names_src) + ".bak")
        print(f"[OK] Backup created: {names_src}.bak")
    if nodes_dst == nodes_src:
        shutil.copy2(nodes_src, str(nodes_src) + ".bak")
        print(f"[OK] Backup created: {nodes_src}.bak")

    append_line(names_src, names_dst, NAMES_LINE, "names.dmp")
    append_line(nodes_src, nodes_dst, NODES_LINE, "nodes.dmp")

    print(f"\nDone. Verify with:")
    print(f"  grep '^{TAXID}' {names_dst}")
    print(f"  grep '^{TAXID}' {nodes_dst}")


if __name__ == "__main__":
    main()
