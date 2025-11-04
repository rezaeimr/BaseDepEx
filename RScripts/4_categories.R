## ============================================================
## 4_categories.R — UpSet plots + MYC category classification
## Works for both gene and isoform levels
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexUpset)
})

## ---------------- Parameters ----------------
analysis_level <- "gene"  # choose: "gene" or "isoform"

## ---------------- Directories ----------------
base_root <- getwd()

de_tbl   <- file.path(base_root, "results/1_DEAnalysis/tables/shrunk")
int_tbl  <- file.path(base_root, "results/2_interaction/tables")
out_tbl  <- file.path(base_root, "results/4_categories/tables")
out_fig  <- file.path(base_root, "results/4_categories/figures")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

## ---------------- Drugs ----------------
drugs <- c("11j", "KVS", "11j_PlaB", "KVS_PlaB", "PlaB")

## ---------------- Helper functions ----------------
strip_version <- function(x) sub("\\..*$", "", x)

load_id_table <- function(file, analysis_level) {
  if (!file.exists(file)) return(data.frame())
  df <- read.delim(file, check.names = FALSE)
  
  if (analysis_level == "gene") {
    if (!"gene_id" %in% names(df)) return(data.frame())
    df <- df %>% dplyr::select(gene_id) %>% distinct()
    df$gene_id <- strip_version(df$gene_id)
  } else if (analysis_level == "isoform") {
    if (!all(c("isoform_id", "gene_id") %in% names(df))) return(data.frame())
    df <- df %>% dplyr::select(isoform_id, gene_id) %>% distinct()
    df$isoform_id <- strip_version(df$isoform_id)
    df$gene_id <- strip_version(df$gene_id)
  }
  return(df)
}

## ============================================================
## PART 1 — UpSet plots
## ============================================================

for (drug in drugs) {
  message(">>> [UpSet] Processing: ", drug)
  
  # --- Load sets ---
  up_drug       <- load_id_table(file.path(de_tbl,  paste0("up_",   drug, ".tsv")), analysis_level)
  up_drug_OHT   <- load_id_table(file.path(de_tbl,  paste0("up_",   drug, "_OHT.tsv")), analysis_level)
  down_drug     <- load_id_table(file.path(de_tbl,  paste0("down_", drug, ".tsv")), analysis_level)
  down_drug_OHT <- load_id_table(file.path(de_tbl,  paste0("down_", drug, "_OHT.tsv")), analysis_level)
  up_int        <- load_id_table(file.path(int_tbl, paste0("up_int_",   drug, ".tsv")), analysis_level)
  down_int      <- load_id_table(file.path(int_tbl, paste0("down_int_", drug, ".tsv")), analysis_level)
  
  ## --- Build combined data frame explicitly ---
  if (analysis_level == "gene") {
    all_ids <- unique(na.omit(c(
      up_drug$gene_id, up_drug_OHT$gene_id, down_drug$gene_id, down_drug_OHT$gene_id,
      up_int$gene_id, down_int$gene_id
    )))
    
    if (length(all_ids) == 0) {
      message("  [skip] No genes found for ", drug)
      next
    }
    
    upset_df <- tibble(gene_id = all_ids)
    ref_col <- upset_df$gene_id
    
    upset_df[[paste0("up_", drug)]]           <- ref_col %in% up_drug$gene_id
    upset_df[[paste0("up_", drug, "_OHT")]]   <- ref_col %in% up_drug_OHT$gene_id
    upset_df[[paste0("down_", drug)]]         <- ref_col %in% down_drug$gene_id
    upset_df[[paste0("down_", drug, "_OHT")]] <- ref_col %in% down_drug_OHT$gene_id
    upset_df[[paste0("up_int_", drug)]]       <- ref_col %in% up_int$gene_id
    upset_df[[paste0("down_int_", drug)]]     <- ref_col %in% down_int$gene_id
    
  } else if (analysis_level == "isoform") {
    all_iso <- unique(na.omit(c(
      up_drug$isoform_id, up_drug_OHT$isoform_id, down_drug$isoform_id, down_drug_OHT$isoform_id,
      up_int$isoform_id, down_int$isoform_id
    )))
    
    if (length(all_iso) == 0) {
      message("  [skip] No isoforms found for ", drug)
      next
    }
    
    # Map isoform_id -> gene_id from all sets
    map_df <- bind_rows(
      up_drug, up_drug_OHT, down_drug, down_drug_OHT, up_int, down_int
    ) %>% distinct(isoform_id, gene_id)
    
    upset_df <- tibble(isoform_id = all_iso) %>%
      left_join(map_df, by = "isoform_id")
    
    ref_col <- upset_df$isoform_id
    
    upset_df[[paste0("up_", drug)]]           <- ref_col %in% up_drug$isoform_id
    upset_df[[paste0("up_", drug, "_OHT")]]   <- ref_col %in% up_drug_OHT$isoform_id
    upset_df[[paste0("down_", drug)]]         <- ref_col %in% down_drug$isoform_id
    upset_df[[paste0("down_", drug, "_OHT")]] <- ref_col %in% down_drug_OHT$isoform_id
    upset_df[[paste0("up_int_", drug)]]       <- ref_col %in% up_int$isoform_id
    upset_df[[paste0("down_int_", drug)]]     <- ref_col %in% down_int$isoform_id
  }
  
  ## --- Save table ---
  write.table(
    upset_df,
    file.path(out_tbl, paste0("upset_data_", drug, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  ## --- Draw UpSet plot ---
  sets <- setdiff(names(upset_df), c("gene_id", "isoform_id"))
  
  p <- upset(
    upset_df,
    intersect = sets,
    base_annotations = list(
      'Intersection size' = intersection_size(counts = TRUE, text = list(size = 3.5))
    ),
    set_sizes = upset_set_size(),
    width_ratio = 0.3,
    height_ratio = 0.3
  ) + theme_bw(base_size = 9)
  
  ggsave(file.path(out_fig, paste0("UpSet_", drug, ".png")), p, width = 9, height = 6, dpi = 300)
}

message("\n>>> Part 1 complete: UpSet plots and tables saved.\n")

## ============================================================
## PART 2 — MYC Category Detection
## ============================================================

for (drug in drugs) {
  message(">>> [Categories] Processing: ", drug)
  
  fp <- file.path(out_tbl, paste0("upset_data_", drug, ".tsv"))
  if (!file.exists(fp)) {
    message("  [skip] No UpSet data for ", drug)
    next
  }
  df <- read.delim(fp, check.names = FALSE)
  
  up_noOHT   <- paste0("up_", drug)
  up_OHT     <- paste0("up_", drug, "_OHT")
  down_noOHT <- paste0("down_", drug)
  down_OHT   <- paste0("down_", drug, "_OHT")
  up_int     <- paste0("up_int_", drug)
  down_int   <- paste0("down_int_", drug)
  
  # Safe accessor
  has_col <- function(x) if (x %in% names(df)) df[[x]] else rep(FALSE, nrow(df))
  U1 <- has_col(up_noOHT); U2 <- has_col(up_OHT)
  D1 <- has_col(down_noOHT); D2 <- has_col(down_OHT)
  UI <- has_col(up_int); DI <- has_col(down_int)
  
  # Classification
  category <- rep(NA_character_, nrow(df))
  category[(U1 & U2 & UI) | (!U1 & U2 & UI)] <- "MYC_enhanced_up"
  category[(D1 & D2 & DI) | (!D1 & D2 & DI)] <- "MYC_enhanced_down"
  category[(U1 & U2 & DI) | (U1 & !U2 & DI)] <- "MYC_suppressed_up"
  category[(D1 & D2 & UI) | (D1 & !D2 & UI)] <- "MYC_suppressed_down"
  category[(D1 & U2 & UI)]                   <- "Switched_increase"
  category[(U1 & D2 & DI)]                   <- "Switched_decrease"
  
  df$category <- category
  
  
  ## --- Save per category ---
  cat_list <- split(df, df$category)
  for (nm in names(cat_list)) {
    sub <- cat_list[[nm]]
    
    # If the subset is a data frame → keep gene_id/isoform_id
    if (is.data.frame(sub)) {
      subset_df <- sub[, intersect(c("gene_id", "isoform_id"), names(sub)), drop = FALSE]
    } else {
      # If it's a vector → convert explicitly to a named data frame
      colname <- if (analysis_level == "gene") "gene_id" else "isoform_id"
      subset_df <- data.frame(setNames(list(sub), colname))
    }
    
    write.table(
      subset_df,
      file.path(out_tbl, paste0(nm, "_", drug, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  
}

cat("\n>>> Part 2 complete: MYC categories saved under results/4_categories/tables\n")


