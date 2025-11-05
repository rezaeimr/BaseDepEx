## ============================================================
## Part 6 — Scatterplots of log2FC with MYC categories
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

## -------------------- Parameters --------------------
analysis_level <- "gene"   # or "isoform"
padj_cutoff <- 0.01

## -------------------- Paths --------------------
base_root <- getwd()
de_tbl   <- file.path(base_root, "results/1_DEAnalysis/tables/shrunk")
cat_tbl  <- file.path(base_root, "results/4_categories/tables")
fig_out  <- file.path(base_root, "results/6_scatterplots/figures")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

## -------------------- Helper functions --------------------
pick_id_col <- function(df) {
  cand <- intersect(c("isoform_id", "gene_id", "feature_id"), names(df))
  if (length(cand) == 0) stop("No ID column found in table.")
  cand[1]
}

load_tt <- function(name) {
  fp <- file.path(de_tbl, paste0("tT_", name, ".tsv"))
  if (!file.exists(fp)) stop("Missing DE table: ", fp)
  read.delim(fp, check.names = FALSE)
}

load_ids <- function(tag, drug, analysis_level) {
  fn <- file.path(cat_tbl, paste0(tag, "_", drug, ".tsv"))
  if (!file.exists(fn)) return(character(0))
  df <- read.delim(fn, check.names = FALSE)
  
  if (analysis_level == "gene" && "gene_id" %in% names(df)) {
    ids <- unique(na.omit(df$gene_id))
  } else if (analysis_level == "isoform" && "isoform_id" %in% names(df)) {
    ids <- unique(na.omit(df$isoform_id))
  } else {
    ids <- character(0)
  }
  ids
}

## -------------------- Drug list --------------------
drug_list <- c("11j", "KVS", "11j_PlaB", "KVS_PlaB", "PlaB")

## -------------------- Colors --------------------
myc_colors <- c(
  "MYC_enhanced_up"     = "#D62728",
  "MYC_enhanced_down"   = "#1F77B4",
  "MYC_suppressed_up"   = "#FF7F0E",
  "MYC_suppressed_down" = "#9467BD",
  "Switched_increase"   = "#2CA02C",
  "Switched_decrease"   = "#17BECF"
)

## -------------------- Main loop --------------------
for (drug in drug_list) {
  message("Plotting: ", drug)
  
  # Load DE tables
  t_no  <- load_tt(drug)
  t_yes <- load_tt(paste0(drug, "_OHT"))
  
  id_no  <- pick_id_col(t_no)
  id_yes <- pick_id_col(t_yes)
  t_no  <- rename(t_no,  feature_id = !!id_no)
  t_yes <- rename(t_yes, feature_id = !!id_yes)
  
  # Merge both tables
  df <- full_join(
    t_no  %>% dplyr::select(feature_id, log2FoldChange, padj) %>%
      rename(log2FC_noOHT = log2FoldChange, padj_noOHT = padj),
    t_yes %>% dplyr::select(feature_id, log2FoldChange, padj) %>%
      rename(log2FC_withOHT = log2FoldChange, padj_withOHT = padj),
    by = "feature_id"
  )
  
  df$category <- "Other"
  
  # Load MYC category IDs depending on analysis level
  enh_up   <- load_ids("MYC_enhanced_up", drug, analysis_level)
  enh_down <- load_ids("MYC_enhanced_down", drug, analysis_level)
  sup_up   <- load_ids("MYC_suppressed_up", drug, analysis_level)
  sup_down <- load_ids("MYC_suppressed_down", drug, analysis_level)
  swi_inc  <- load_ids("Switched_increase", drug, analysis_level)
  swi_dec  <- load_ids("Switched_decrease", drug, analysis_level)
  
  df$category[df$feature_id %in% enh_up]   <- "MYC_enhanced_up"
  df$category[df$feature_id %in% enh_down] <- "MYC_enhanced_down"
  df$category[df$feature_id %in% sup_up]   <- "MYC_suppressed_up"
  df$category[df$feature_id %in% sup_down] <- "MYC_suppressed_down"
  df$category[df$feature_id %in% swi_inc]  <- "Switched_increase"
  df$category[df$feature_id %in% swi_dec]  <- "Switched_decrease"
  
  label_positions <- df %>%
    filter(category != "Other") %>%
    group_by(category) %>%
    summarise(
      x = median(log2FC_noOHT, na.rm = TRUE),
      y = median(log2FC_withOHT, na.rm = TRUE),
      .groups = "drop"
    )
  
  ## -------------------- Plot --------------------
  p <- ggplot(df, aes(x = log2FC_noOHT, y = log2FC_withOHT)) +
    geom_point(data = subset(df, category == "Other"),
               color = "gray85", alpha = 0.6, size = 1.1) +
    geom_point(data = subset(df, category != "Other"),
               aes(color = category), alpha = 0.85, size = 1.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = myc_colors) +
    labs(title = paste0("Drug: ", drug),
         x = expression("log"[2]*"FC (Drug vs DMSO, OHT OFF)"),
         y = expression("log"[2]*"FC (Drug vs DMSO, OHT ON)"),
         color = "Category") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
  
  ggsave(file.path(fig_out, paste0("scatter_", drug, "_cat.png")),
         p, width = 6, height = 5, dpi = 300)
  
  message("  Saved: scatter_", drug, "_cat.png")
}

message("Scatterplots with category labels complete: ", fig_out)
