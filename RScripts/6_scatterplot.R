## ============================================================
## 6_scatterplots.R — Scatterplots of log2FC with MYC categories
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

set.seed(1)

## -------------------- Directories --------------------
project_root <- getwd()

de_tbl  <- file.path(project_root, "results/1_DEAnalysis/tables/shrunk")
cat_tbl <- file.path(project_root, "results/4_categories/tables")
fig_out <- file.path(project_root, "results/6_scatterplots/figures")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

## -------------------- Drugs --------------------
drug_list <- c("11j", "KVS", "11j_PlaB", "KVS_PlaB", "PlaB")

## -------------------- Helper functions --------------------
pick_id_col <- function(df) {
  cand <- intersect(c("gene_id", "ensembl_id", "transcript_id", "feature_id", "isoform_id"), names(df))
  if (length(cand) == 0) stop("No ID column found in table.")
  cand[1]
}

load_tt <- function(name) {
  fp <- file.path(de_tbl, paste0("tT_", name, ".tsv"))
  if (!file.exists(fp)) stop("Missing DE table: ", fp)
  read.delim(fp, check.names = FALSE)
}

load_ids <- function(tag, drug) {
  fn <- file.path(cat_tbl, paste0(tag, "_", drug, ".tsv"))
  if (!file.exists(fn)) return(character(0))
  df <- read.delim(fn, check.names = FALSE)
  cand <- intersect(c("gene_id", "ensembl_id", "feature_id", "isoform_id"), names(df))
  if (length(cand) == 0) return(character(0))
  unique(df[[cand[1]]])
}

## -------------------- Colors --------------------
myc_colors <- c(
  "MYC_enhanced_up"     = "#E41A1C",
  "MYC_enhanced_down"   = "#377EB8",
  "MYC_suppressed_up"   = "#FF7F00",
  "MYC_suppressed_down" = "#4DAF4A",
  "Switched_increase"   = "#984EA3",
  "Switched_decrease"   = "#A65628",
  "Other"               = "gray85"
)

## ============================================================
## Main plotting loop
## ============================================================

for (drug in drug_list) {
  message(">>> Plotting: ", drug)
  
  # --- Load DE tables (drug and drug + OHT)
  t_no  <- load_tt(drug)
  t_yes <- load_tt(paste0(drug, "_OHT"))
  
  id_no  <- pick_id_col(t_no)
  id_yes <- pick_id_col(t_yes)
  
  t_no  <- rename(t_no,  feature_id = !!id_no)
  t_yes <- rename(t_yes, feature_id = !!id_yes)
  
  df <- full_join(
    t_no  %>% dplyr::select(feature_id, log2FoldChange, padj) %>%
      rename(log2FC_noOHT = log2FoldChange, padj_noOHT = padj),
    t_yes %>% dplyr::select(feature_id, log2FoldChange, padj) %>%
      rename(log2FC_withOHT = log2FoldChange, padj_withOHT = padj),
    by = "feature_id"
  )
  
  df$category <- "Other"
  
  # --- Load MYC categories ---
  enh_up   <- load_ids("MYC_enhanced_up", drug)
  enh_down <- load_ids("MYC_enhanced_down", drug)
  sup_up   <- load_ids("MYC_suppressed_up", drug)
  sup_down <- load_ids("MYC_suppressed_down", drug)
  swi_inc  <- load_ids("Switched_increase", drug)
  swi_dec  <- load_ids("Switched_decrease", drug)
  
  # --- Assign categories ---
  df$category[df$feature_id %in% enh_up]   <- "MYC_enhanced_up"
  df$category[df$feature_id %in% enh_down] <- "MYC_enhanced_down"
  df$category[df$feature_id %in% sup_up]   <- "MYC_suppressed_up"
  df$category[df$feature_id %in% sup_down] <- "MYC_suppressed_down"
  df$category[df$feature_id %in% swi_inc]  <- "Switched_increase"
  df$category[df$feature_id %in% swi_dec]  <- "Switched_decrease"
  
  # --- Compute label positions ---
  label_positions <- df %>%
    filter(category != "Other") %>%
    group_by(category) %>%
    summarise(
      x = median(log2FC_noOHT, na.rm = TRUE),
      y = median(log2FC_withOHT, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- Scatterplot ---
  p <- ggplot(df, aes(x = log2FC_noOHT, y = log2FC_withOHT)) +
    geom_point(data = subset(df, category == "Other"),
               color = "gray85", alpha = 0.6, size = 1.1) +
    geom_point(data = subset(df, category != "Other"),
               aes(color = category), alpha = 0.9, size = 1.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = myc_colors, name = "Category") +
    labs(
      title = paste0("MYC response — ", drug),
      x = expression("log"[2]*"FC (Drug vs DMSO, OHT OFF)"),
      y = expression("log"[2]*"FC (Drug vs DMSO, OHT ON)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  ggsave(file.path(fig_out, paste0("Scatter_MYCCat_", drug, ".png")),
         p, width = 6.5, height = 5.5, dpi = 300)
  
  message("    Saved → Scatter_MYCCat_", drug, ".png")
}

message("\n>>> 6_scatterplots complete → ", fig_out, "\n")

