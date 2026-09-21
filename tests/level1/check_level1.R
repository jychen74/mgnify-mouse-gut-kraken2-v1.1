#!/usr/bin/env Rscript
# =============================================================================
# Level 1 check - reproduces the read-level numbers of the GF-mouse (barcode02)
# three-layer validation from the shipped EPI2ME outputs. NO database needed.
# Base R only. Usage (from any directory):  Rscript check_level1.R
# Exit status 0 = all metrics match tests/level1/expected/expected_counts.tsv
# =============================================================================
args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- if (length(args_file)) dirname(normalizePath(args_file)) else getwd()
dat  <- function(f) file.path(here, "data", f)

read_assign <- function(f) {
  stopifnot(file.exists(f))
  read.delim(f, header = FALSE, quote = "", comment.char = "", fill = TRUE,
             stringsAsFactors = FALSE,
             col.names = c("status", "read_id", "taxid", "len", "kmers", "lineage"),
             colClasses = c("character", "character", "integer",
                            "integer", "character", "character"))
}
read_bracken <- function(f) {
  stopifnot(file.exists(f))
  read.delim(f, quote = "", comment.char = "", stringsAsFactors = FALSE)
}

I   <- read_assign(dat("layerI_rawfastq_v1.0.assignments.tsv"))
II  <- read_assign(dat("layerII_nohost_v1.0.assignments.tsv"))
III <- read_assign(dat("layerIII_nohost_v1.1.assignments.tsv"))
bI   <- read_bracken(dat("layerI_rawfastq_v1.0.bracken.report"))
bII  <- read_bracken(dat("layerII_nohost_v1.0.bracken.report"))
bIII <- read_bracken(dat("layerIII_nohost_v1.1.bracken.report"))

cls10 <- function(x) ifelse(x$status == "U", "U", ifelse(x$taxid == 1, "root", "bac"))
cls11 <- function(x) ifelse(x$status == "U", "U",
                     ifelse(x$taxid == 3734, "mus", ifelse(x$taxid == 1, "root", "bac")))

m <- merge(data.frame(read_id = II$read_id,  a = cls10(II),  stringsAsFactors = FALSE),
           data.frame(read_id = III$read_id, b = cls11(III), stringsAsFactors = FALSE),
           by = "read_id")
flow <- function(a, b) sum(m$a == a & m$b == b)
mus_bracken <- bIII$new_est_reads[bIII$taxonomy_id == 3734]

obs <- c(
  layerI_total_reads                    = nrow(I),
  layerI_classified                     = sum(I$status == "C"),
  layerI_unclassified                   = sum(I$status == "U"),
  layerI_removed_by_minimap2            = sum(!I$read_id %in% II$read_id),
  layerII_total_reads                   = nrow(II),
  layerII_reads_found_in_layerI         = sum(II$read_id %in% I$read_id),
  layerII_classified_v1.0               = sum(II$status == "C"),
  layerII_unclassified_v1.0             = sum(II$status == "U"),
  layerIII_total_reads                  = nrow(III),
  layerIII_same_read_set_as_layerII     = as.integer(setequal(II$read_id, III$read_id)),
  layerIII_unclassified_v1.1            = sum(III$status == "U"),
  layerIII_taxid3734_reads              = sum(cls11(III) == "mus"),
  layerIII_bacterial_reads              = sum(cls11(III) == "bac"),
  layerIII_root_reads                   = sum(cls11(III) == "root"),
  layerIII_classified_total             = sum(III$status == "C"),
  flow_bacteria_to_taxid3734            = flow("bac", "mus"),
  flow_unclassified_to_taxid3734        = flow("U", "mus"),
  flow_bacteria_to_bacteria             = flow("bac", "bac"),
  flow_bacteria_to_root                 = flow("bac", "root"),
  flow_unclassified_to_bacteria_or_root = flow("U", "bac") + flow("U", "root"),
  flow_unclassified_to_unclassified     = flow("U", "U"),
  bracken_layerI_total_est_reads        = sum(bI$new_est_reads),
  bracken_layerII_total_est_reads       = sum(bII$new_est_reads),
  bracken_layerIII_total_est_reads      = sum(bIII$new_est_reads),
  bracken_layerIII_mus_est_reads        = mus_bracken,
  bracken_layerIII_bacterial_est_reads  = sum(bIII$new_est_reads) - mus_bracken,
  bracken_layerIII_mus_pct_2dp          = round(100 * mus_bracken / sum(bIII$new_est_reads), 2)
)

exp <- read.delim(file.path(here, "expected", "expected_counts.tsv"),
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
stopifnot(setequal(exp$metric, names(obs)))
res <- data.frame(metric   = exp$metric,
                  expected = exp$expected,
                  observed = unname(obs[exp$metric]),
                  stringsAsFactors = FALSE)
res$status <- ifelse(abs(res$expected - res$observed) < 1e-9, "PASS", "FAIL")

print(res, row.names = FALSE)
cat(sprintf("\n%d / %d metrics PASS\n", sum(res$status == "PASS"), nrow(res)))
if (any(res$status != "PASS")) quit(status = 1)
