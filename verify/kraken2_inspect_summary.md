# Kraken2 database v1.1 — verification results

**Build date:** August 2026  
**Kraken2 version used for build:** 2.17.1  
**K-mer parameters:** k = 35, minimizer length l = 31 (inherited from v1.0 `opts.k2d`)

## Host k-mer coverage

After build, `kraken2-inspect` confirmed that the custom *Mus musculus* node (taxid 3734) accounts for **23.55%** of all k-mers in the completed v1.1 database, consistent with successful indexing of the GRCm39 reference genome assembly (GCF_000001635.27; 61 sequences, ~2.6 Gb).

## Three-layer validation (GF negative-control, exp_12 barcode02)

| Layer | Input | Database | *Mus musculus* reads | Bacterial FP reads |
|-------|-------|----------|---------------------|-------------------|
| ① Before depletion | fastq_pass | v1.0 | 0% (taxid absent) | 6,552 |
| ② After depletion  | no_host_fastq | v1.0 | 0% (taxid absent) | 294 |
| ③ After depletion  | no_host_fastq | **v1.1** | **98.75%** | ~0 |

Classification performed with Kraken2 v2.1.3 (EPI2ME wf-metagenomics v2.14.2).  
Layer ③ uses identical input to Layer ②; only the database version differs.

## Cross-batch validation (exp_14, Q2HM cohort, 6 SPF mice, independent facility)

All six barcodes yielded detectable *Mus musculus* reads (taxid 3734) at the species level with v1.1, ranging from 0.22% to 3.53% of species-level-classifiable reads. No *Mus musculus* reads were detectable with v1.0.

See the associated publication for full details.
