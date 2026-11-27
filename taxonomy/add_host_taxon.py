#!/usr/bin/env python3
"""
add_host_taxon.py
-----------------
Appends a custom host taxonomy node to the taxonomy files of any GTDB-derived
Kraken2 database, enabling correct host-read classification in host-associated
metagenomic workflows.

Background:
    GTDB taxonomy covers prokaryotes only. Any Kraken2/Bracken database built
    from GTDB-derived MAGs therefore lacks eukaryotic host lineages. Without an
    explicit host node, residual host reads that evade upstream host-depletion
    are either left unclassified or mis-assigned to microbial taxa. This script
    resolves that structural gap by appending a single host node directly under
    root (taxid 1), which is necessary because GTDB does not include eukaryotic
    superkingdom nodes.

Usage (mouse — reproduces v1.1 build):
    python add_host_taxon.py \\
        --names  /path/to/taxonomy/names.dmp \\
        --nodes  /path/to/taxonomy/nodes.dmp

Usage (any host):
    python add_host_taxon.py \\
        --names        /path/to/taxonomy/names.dmp \\
        --nodes        /path/to/taxonomy/nodes.dmp \\
        --species-name "Rattus norvegicus" \\
        --taxid        3735

    Choosing a safe taxid:
        Use a sequential integer that is not already present in the database's
        nodes.dmp (check with: awk -F'\\t' '{print $1}' nodes.dmp | sort -n | tail).
        Avoid NCBI taxids for the same organism (Mus musculus = 10090,
        Rattus norvegicus = 10116) to prevent downstream NCBI-lineage conflicts.

Options:
    --names         Path to input names.dmp (modified in-place; .bak created)
    --nodes         Path to input nodes.dmp (modified in-place; .bak created)
    --out-names     Optional separate output path for names.dmp
    --out-nodes     Optional separate output path for nodes.dmp
    --species-name  Scientific name of the host (default: "Mus musculus")
    --taxid         Custom taxid to assign (default: 3734)
    --parent-taxid  Parent taxid in the GTDB tree (default: 1 = root)
    --rank          Taxonomic rank string (default: "species")

Reference:
    Chen J-Y (2026). Extending the MGnify mouse-gut Kraken2 database
    with host taxonomy greatly reduces systematic false-positive classifications.
    Microbial Genomics (Technical Resource), in preparation.
"""

import argparse
import shutil
import sys
from pathlib import Path


def build_names_line(taxid: str, species_name: str) -> str:
    return "\t|\t".join([taxid, species_name, "", "scientific name"]) + "\t|"


def build_nodes_line(taxid: str, parent_taxid: str, rank: str) -> str:
    return "\t|\t".join([
        taxid,        # taxid
        parent_taxid, # parent taxid
        rank,         # rank
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


def is_taxid_present(filepath: Path, taxid: str) -> bool:
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(taxid + "\t"):
                return True
    return False


def append_line(src: Path, dst: Path, line: str, label: str) -> None:
    if dst != src:
        shutil.copy2(src, dst)

    if is_taxid_present(dst, line.split("\t")[0]):
        print(f"[WARNING] taxid already present in {dst.name}. No change made.")
        return

    with open(dst, "r") as f:
        content = f.read()

    if content and not content.endswith("\n"):
        with open(dst, "a") as f:
            f.write("\n")

    with open(dst, "a") as f:
        f.write(line + "\n")

    print(f"[OK] Appended {label} entry to {dst}")


def main():
    parser = argparse.ArgumentParser(
        description="Add a custom host taxonomy node to a GTDB-derived Kraken2 database.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--names",        required=True,            help="Path to names.dmp")
    parser.add_argument("--nodes",        required=True,            help="Path to nodes.dmp")
    parser.add_argument("--out-names",    default=None,             help="Output names.dmp (default: in-place)")
    parser.add_argument("--out-nodes",    default=None,             help="Output nodes.dmp (default: in-place)")
    parser.add_argument("--species-name", default="Mus musculus",   help='Species name (default: "Mus musculus")')
    parser.add_argument("--taxid",        default="3734",           help="Custom taxid (default: 3734)")
    parser.add_argument("--parent-taxid", default="1",              help="Parent taxid (default: 1 = root)")
    parser.add_argument("--rank",         default="species",        help='Taxonomic rank (default: "species")')
    args = parser.parse_args()

    names_src = Path(args.names)
    nodes_src = Path(args.nodes)

    if not names_src.exists():
        sys.exit(f"[ERROR] names.dmp not found: {names_src}")
    if not nodes_src.exists():
        sys.exit(f"[ERROR] nodes.dmp not found: {nodes_src}")

    names_dst = Path(args.out_names) if args.out_names else names_src
    nodes_dst = Path(args.out_nodes) if args.out_nodes else nodes_src

    names_dst.parent.mkdir(parents=True, exist_ok=True)
    nodes_dst.parent.mkdir(parents=True, exist_ok=True)

    # Backup if in-place
    if names_dst == names_src:
        shutil.copy2(names_src, str(names_src) + ".bak")
        print(f"[OK] Backup: {names_src}.bak")
    if nodes_dst == nodes_src:
        shutil.copy2(nodes_src, str(nodes_src) + ".bak")
        print(f"[OK] Backup: {nodes_src}.bak")

    names_line = build_names_line(args.taxid, args.species_name)
    nodes_line = build_nodes_line(args.taxid, args.parent_taxid, args.rank)

    append_line(names_src, names_dst, names_line, "names.dmp")
    append_line(nodes_src, nodes_dst, nodes_line, "nodes.dmp")

    print(f"\nDone. Verify with:")
    print(f"  grep '^{args.taxid}' {names_dst}")
    print(f"  grep '^{args.taxid}' {nodes_dst}")


if __name__ == "__main__":
    main()
