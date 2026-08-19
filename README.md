# MGnify Mouse-Gut Kraken2/Bracken Database v1.1

This repository provides scripts and taxonomy files to build a corrected version (v1.1) of the [MGnify mouse-gut Kraken2/Bracken pre-built database](https://www.ebi.ac.uk/metagenomics/), which extends the original v1.0 release with a *Mus musculus* host taxonomy node.

> **Associated publication:** Chen J-Y et al. (2026). A missing host taxonomy node in the MGnify mouse-gut Kraken2 database causes systematic false-positive metagenomic classifications: identification, root-cause analysis, and a validated corrected database. *PeerJ* [DOI: pending]

---

## Background

The MGnify mouse-gut Kraken2 database v1.0 is built on a taxonomy derived from the [Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org/) release 214, which covers exclusively prokaryotic lineages. As a result, *Mus musculus* has no representation in the v1.0 taxonomy tree.

When residual host reads persist in sequencing data after minimap2 host depletion (as routinely occurs), these reads have no valid taxonomic destination and are systematically mis-assigned to the most k-mer-similar bacterial nodes. This creates false-positive bacterial classifications that are particularly pronounced in germ-free (GF) negative-control animals, but are present in all mouse gut metagenomic datasets.

**v1.1 adds a custom *Mus musculus* taxonomy node (taxid 3734) and indexes the GRCm39 reference genome assembly**, so residual host reads are correctly assigned rather than distributed across bacterial taxa.

> **Note on the custom taxid:** The taxid 3734 used in this database is a custom sequential identifier within the MGnify-derived taxonomy tree. It is **numerically distinct** from NCBI Taxonomy ID 10090 (the standard NCBI identifier for *Mus musculus*). Users integrating v1.1 output with NCBI-based tools should remap this taxid accordingly.

---

## Repository structure

```
├── README.md
├── LICENSE
├── taxonomy/
│   ├── add_host_taxon.py        # Script to add taxid 3734 to v1.0 taxonomy files
│   └── host_taxonomy_patch.txt  # Exact lines added to names.dmp and nodes.dmp
├── build/
│   └── build_v1.1.sh            # Complete database build script
└── verify/
    └── kraken2_inspect_summary.md  # Post-build verification results
```

---

## How to build v1.1

### Prerequisites

- MGnify mouse-gut Kraken2 database v1.0 (downloaded from [MGnify FTP](https://ftp.ebi.ac.uk/pub/databases/metagenomics/))
- Kraken2 (tested with v2.17.1)
- Bracken (tested with v2.x)
- GRCm39 *Mus musculus* reference genome: `GCF_000001635.27_GRCm39_genomic.fna` (available from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/))
- Python 3 (for the taxonomy modification script)

### Step 1 — Modify the taxonomy files

Navigate to your v1.0 database directory and run:

```bash
python taxonomy/add_host_taxon.py \
  --names  /path/to/mgnify_mouse_gut_v1.0/taxonomy/names.dmp \
  --nodes  /path/to/mgnify_mouse_gut_v1.0/taxonomy/nodes.dmp \
  --out-names  /path/to/mgnify_mouse_gut_v1.1/taxonomy/names.dmp \
  --out-nodes  /path/to/mgnify_mouse_gut_v1.1/taxonomy/nodes.dmp
```

Or manually append the lines documented in `taxonomy/host_taxonomy_patch.txt`.

### Step 2 — Add the GRCm39 genome to the library

```bash
# Prepare FASTA headers for Kraken2 (assign taxid 3734 to all sequences)
# The add_host_taxon.py script handles header formatting automatically.
# Alternatively, use sed:
sed 's/^>\(\S*\)/>\1|kraken:taxid|3734/' GCF_000001635.27_GRCm39_genomic.fna > GRCm39_k2.fa

kraken2-build --add-to-library GRCm39_k2.fa \
  --db /path/to/mgnify_mouse_gut_v1.1
```

### Step 3 — Build the database

```bash
# K-mer parameters (k=35, minimizer length l=31) are inherited from v1.0 opts.k2d
kraken2-build --build --db /path/to/mgnify_mouse_gut_v1.1 --threads 8
```

### Step 4 — Build Bracken databases

```bash
for len in 50 100 150 200 250; do
  bracken-build -d /path/to/mgnify_mouse_gut_v1.1 -k 35 -l $len -t 8
done
```

### Step 5 — Verify

```bash
kraken2-inspect --db /path/to/mgnify_mouse_gut_v1.1 \
  --report-zero-counts > kraken2_inspect_v1.1.txt

grep "Mus musculus" kraken2_inspect_v1.1.txt
```

Expected: taxid 3734 should account for ~23.55% of all k-mers.

See `verify/kraken2_inspect_summary.md` for our build verification results.

---

## Expected behaviour after correction

| Condition | v1.0 result | v1.1 result |
|-----------|-------------|-------------|
| GF negative-control sample | *Mus musculus* reads: 0% (taxid absent); false-positive bacterial reads present | *Mus musculus* reads: ~98.75% of species-classifiable reads; bacterial FP reads: ~0 |
| SPF mouse samples | Host reads dispersed across bacterial taxa | Host reads correctly assigned to *Mus musculus* (taxid 3734) |

---

## Display artefact in EPI2ME

Because the *Mus musculus* node (taxid 3734) is placed directly under root (taxid 1) — an unavoidable consequence of GTDB not including eukaryotic superkingdom nodes — EPI2ME abundance tables will display the label **"Mus musculus_Incertae_sedis"** in the superkingdom and kingdom fields. This is a display artefact only and does not affect Kraken2 read-level classification or species-level abundance estimates.

---

## Citation

If you use this corrected database or these scripts, please cite:

> Chen J-Y et al. (2026). A missing host taxonomy node in the MGnify mouse-gut Kraken2 database causes systematic false-positive metagenomic classifications: identification, root-cause analysis, and a validated corrected database. *PeerJ* [DOI: pending]

Please also cite the original MGnify database:

> Richardson L et al. (2023). MGnify: the microbiome sequence data analysis resource in 2023. *Nucleic Acids Research* 51(D1):D753–D759. doi:10.1093/nar/gkac1080

---

## License

Scripts in this repository are released under the MIT License. The underlying MGnify database v1.0 retains its original license from EMBL-EBI.
