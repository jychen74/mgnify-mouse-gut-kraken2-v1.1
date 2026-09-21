# Test package - Level 1 (no database required)

Reproduces the read-level numbers behind the three-layer validation of the
germ-free (GF) mouse negative control (Q3HM run, barcode02) from the shipped
EPI2ME wf-metagenomics outputs. It does **not** re-run Kraken2, so it needs no
database and no raw fastq files; it verifies the counting and the figure code.
Level 2 (re-classification with the v1.1 database) and Level 3 (full rebuild)
are described in the top-level README.

## Contents

| Path | What it is |
|------|------------|
| `data/*.assignments.tsv` | Unmodified EPI2ME `reads_assignments` files (Kraken2 read-level calls; columns: status, read ID, taxid, length, k-mer string, lineage) |
| `data/*.bracken.report` | Unmodified EPI2ME Bracken reports (`new_est_reads`, threshold 10) |
| `data/MD5SUMS` | Checksums of the files above (`cd data && md5sum -c MD5SUMS`) |
| `provenance/*.versions_all.txt` | Tool versions reported by the workflow for each layer |
| `provenance/*.params.json` | Selected workflow parameters (instance ID, input, database, `kraken2_confidence`, Bracken threshold) |
| `expected/expected_counts.tsv` | Expected value for every metric checked |
| `check_level1.R` | Recomputes every metric and compares it with the expected value (base R only) |
| `fig4_readlevel_sankey.R` | Regenerates the read-level flow figure (Layer II to Layer III) |

Layers: **I** = raw `fastq_pass` + database v1.0; **II** = `no_host_fastq`
(after minimap2 host depletion) + v1.0; **III** = the same `no_host_fastq`
reads + v1.1.

## Run

```bash
cd tests/level1
Rscript check_level1.R          # expect: 27 / 27 metrics PASS, exit status 0
Rscript fig4_readlevel_sankey.R # writes output/fig4_readlevel_sankey.{pdf,tiff,png} + data CSVs
```

Verified with R 4.3.3 on Linux (fresh copy of this folder; all files pass `md5sum -c`).
`fig4_readlevel_sankey.R` needs the R packages `ggplot2`, `dplyr`, `ggalluvial`,
`ggpubr`, `ggrepel`, `scales` (and `Cairo` support for the PDF export).

## Notes on what the numbers are

- Kraken2 read-level counts (`*.assignments.tsv`) and Bracken-estimated counts
  (`*.bracken.report`) are different quantities. Bracken redistributes reads
  and applies a 10-read threshold, so for Layer III 1,214 reads are assigned to
  taxid 3734 by Kraken2, whereas Bracken estimates 1,265 of 1,280 reads
  (98.83%).
- Layer II and Layer III contain the identical 6,184 reads; only the database differs.
- `kraken2_confidence` was 0 in all three runs.
