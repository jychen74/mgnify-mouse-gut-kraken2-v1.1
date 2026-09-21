# =============================================================================
# Fig 4 — Read-level flow of GF mouse (barcode02) reads across Layers I → III
# Paper: "Extending the MGnify mouse-gut Kraken2/Bracken database..."
# Input : Kraken2 read-level assignments shipped in tests/level1/data (EPI2ME wf-metagenomics reads_assignments)
#   Layer I  : raw fastq_pass      + v1.0   (01M2MQW75DGVM9P5HHSHR3SM43)
#   Layer II : no_host_fastq       + v1.0   (01KZR0HVM2BXXEVEZ3S1552WVG)
#   Layer III: no_host_fastq       + v1.1   (01KZWY3BS86J6DDSEKPYSRMZK0)
# Panel A: Layer I classification -> fate after minimap2 depletion
# Panel B: the SAME 6,184 retained reads, Layer II (v1.0) -> Layer III (v1.1)
# All counts are Kraken2 read-level assignments (NOT Bracken re-estimates).
# Panel B: 3 reads (Unclassified -> Bacteria/Root) are drawn dark grey, not in legend.
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggalluvial)
library(ggpubr)     # theme_pubr(), ggarrange()
library(ggrepel)    # installed as a ggpubr dependency
library(scales)

# ── 0. Paths (run from tests/level1/, or set MGNIFY_L1 to that folder) ───────
here    <- Sys.getenv("MGNIFY_L1", ".")
stopifnot(dir.exists(file.path(here, "data")))
out_dir <- file.path(here, "output")
dir.create(out_dir, showWarnings = FALSE)

runs <- c(I   = "layerI_rawfastq_v1.0",
          II  = "layerII_nohost_v1.0",
          III = "layerIII_nohost_v1.1")

# ── 1. Read assignments ───────────────────────────────────────────────────────
read_assign <- function(run_name) {
  f <- file.path(here, "data", paste0(run_name, ".assignments.tsv"))
  stopifnot(file.exists(f))
  d <- read.delim(f, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  fill = TRUE, stringsAsFactors = FALSE,
                  col.names = c("status", "read_id", "taxid", "len", "kmers", "lineage"),
                  colClasses = c("character", "character", "integer",
                                 "integer", "character", "character"))
  d[, c("status", "read_id", "taxid")]
}
I   <- read_assign(runs["I"])
II  <- read_assign(runs["II"])
III <- read_assign(runs["III"])

# ── 2. Sanity checks (stop if the input is not what the paper describes) ─────
stopifnot(nrow(I) == 34830, nrow(II) == 6184, nrow(III) == 6184,
          all(II$read_id %in% I$read_id), setequal(II$read_id, III$read_id))

cls_v10 <- function(d) ifelse(d$status == "U", "Unclassified",
                       ifelse(d$taxid == 1, "Root (ambiguous)", "Bacteria"))
cls_v11 <- function(d) ifelse(d$status == "U", "Unclassified",
                       ifelse(d$taxid == 3734, "Mus musculus",
                       ifelse(d$taxid == 1, "Root (ambiguous)", "Bacteria")))

# ── 3. Panel A data: Layer I class -> minimap2 fate ──────────────────────────
A <- I %>%
  mutate(`Layer I\n(raw, v1.0)` = cls_v10(I),
         `minimap2\nhost depletion` = ifelse(read_id %in% II$read_id,
                                             "Retained", "Removed")) %>%
  count(`Layer I\n(raw, v1.0)`, `minimap2\nhost depletion`, name = "n")

# ── 4. Panel B data: Layer II -> Layer III on identical reads ────────────────
B <- inner_join(
  data.frame(read_id = II$read_id,  II  = cls_v10(II),  stringsAsFactors = FALSE),
  data.frame(read_id = III$read_id, III = cls_v11(III), stringsAsFactors = FALSE),
  by = "read_id") %>%
  mutate(flow = case_when(
    II == "Bacteria"     & III == "Mus musculus"     ~ "Bacteria -> Mus musculus",
    II == "Unclassified" & III == "Mus musculus"     ~ "Unclassified -> Mus musculus",
    II == "Bacteria"     & III == "Bacteria"         ~ "Bacteria -> Bacteria (residual)",
    II == "Bacteria"     & III == "Root (ambiguous)" ~ "Bacteria -> Root (ambiguous)",
    II == "Unclassified" & III == "Unclassified"     ~ "Unclassified -> Unclassified",
    TRUE                                             ~ "Other"),
    II  = factor(II,  levels = c("Bacteria", "Unclassified")),
    III = factor(III, levels = c("Mus musculus", "Bacteria",
                                 "Root (ambiguous)", "Unclassified"))) %>%
  count(II, III, flow, name = "n")

stopifnot(sum(B$n) == 6184,
          B$n[B$flow == "Bacteria -> Mus musculus"] +
          B$n[B$flow == "Unclassified -> Mus musculus"] == 1214)   # Kraken2 Mus reads

# Save the numbers behind the figure (same habit as fig2_data.csv)
write.csv(A, file.path(out_dir, "fig4_panelA_data.csv"), row.names = FALSE)
write.csv(B, file.path(out_dir, "fig4_panelB_data.csv"), row.names = FALSE)

# ── 5. Colours (bacteria = orange, Mus = green: consistent with Fig 1) ───────
col_bac <- "#D55E00"; col_mus <- "#009E73"; col_mus_l <- "#8FD6C0"
col_root <- "#7B6FB0"; col_un <- "#B8B8B8"

# ASCII "->" on purpose: avoids missing-glyph problems on some fonts/devices
flow_cols <- c("Bacteria -> Mus musculus"        = col_mus,
               "Unclassified -> Mus musculus"    = col_mus_l,
               "Bacteria -> Bacteria (residual)" = col_bac,
               "Bacteria -> Root (ambiguous)"    = col_root,
               "Unclassified -> Unclassified"    = col_un,
               "Other"                           = "#666666")

# Stratum labels: left of axis 1, right of axis 2 (so they never sit on the flows)
strip_sfx <- function(x) sub("\\.1$", "", as.character(x))
lab_layer <- function(axis_idx, nudge, hjust, ...) {
  geom_text_repel(
    stat = "stratum", ...,
    aes(label = ifelse(after_stat(x) == axis_idx,
                       paste0(strip_sfx(after_stat(stratum)), "\n",
                              comma(after_stat(count))), NA)),
    size = 4.4, direction = "y", nudge_x = nudge, hjust = hjust,
    segment.size = 0.25, min.segment.length = 0.1, lineheight = 0.9,
    na.rm = TRUE)   # NA labels are intentional (other axis) -> silence the warning
}
x_exp <- expansion(add = c(0.9, 0.9))

# ── 6. Panel A ────────────────────────────────────────────────────────────────
pA <- ggplot(A, aes(axis1 = `Layer I\n(raw, v1.0)`,
                    axis2 = `minimap2\nhost depletion`, y = n)) +
  geom_flow(aes(fill = `Layer I\n(raw, v1.0)`), width = 0.14, alpha = 0.8) +
  geom_stratum(width = 0.14, fill = "white", colour = "grey30", linewidth = 0.4) +
  lab_layer(1, -0.27, 1) + lab_layer(2, 0.27, 0) +
  scale_x_discrete(limits = c("Layer I\n(raw, v1.0)", "minimap2\nhost depletion"),
                   expand = x_exp) +
  scale_y_continuous(labels = comma, expand = c(0, 0)) +
  scale_fill_manual(values = c("Bacteria" = col_bac, "Unclassified" = col_un,
                               "Root (ambiguous)" = col_root), name = NULL) +
  labs(y = "Reads (Kraken2)", x = NULL,
       title = "Host depletion, GF barcode02") +
  theme_pubr(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12))

# ── 7. Panel B ────────────────────────────────────────────────────────────────
pB <- ggplot(B, aes(axis1 = II, axis2 = III, y = n)) +
  geom_flow(aes(fill = flow), width = 0.14, alpha = 0.85, discern = TRUE) +
  geom_stratum(width = 0.14, fill = "white", colour = "grey30", linewidth = 0.4,
               discern = TRUE) +
  lab_layer(1, -0.27, 1, discern = TRUE) + lab_layer(2, 0.27, 0, discern = TRUE) +
  scale_x_discrete(limits = c("Layer II\n(depleted, v1.0)", "Layer III\n(depleted, v1.1)"),
                   expand = x_exp) +
  scale_y_continuous(labels = comma, expand = c(0, 0)) +
  scale_fill_manual(values = flow_cols, name = NULL, breaks = names(flow_cols)[1:5]) +
  labs(y = "Reads (Kraken2)", x = NULL,
       title = "Same 6,184 reads, v1.0 vs v1.1") +
  theme_pubr(base_size = 14) +
  theme(legend.position = "bottom", legend.text = element_text(size = 11),
        legend.key.size = unit(0.45, "cm"), legend.spacing.x = unit(0.6, "cm"),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12)) +
  guides(fill = guide_legend(ncol = 1))

# ── 8. Combine & export ───────────────────────────────────────────────────────
fig4 <- ggarrange(pA, pB, ncol = 2, widths = c(1, 1.25), labels = c("A", "B"))
print(fig4)

fw <- 30; fh <- 16
ggsave(file.path(out_dir, "fig4_readlevel_sankey.pdf"),  fig4, width = fw, height = fh,
       units = "cm", device = cairo_pdf)
ggsave(file.path(out_dir, "fig4_readlevel_sankey.tiff"), fig4, width = fw, height = fh,
       units = "cm", dpi = 600, compression = "lzw")
ggsave(file.path(out_dir, "fig4_readlevel_sankey.png"),  fig4, width = fw, height = fh,
       units = "cm", dpi = 300, bg = "white")
message("Done. Exported fig4_readlevel_sankey (pdf/tiff/png) + fig4_panelA/B_data.csv")
