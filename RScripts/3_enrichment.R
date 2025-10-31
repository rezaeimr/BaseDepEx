## ============================================================
## 3_enrichment.R — GSEA + ORA (gene_id-based)
## Color rules + ratio labels
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
})

set.seed(1)
message(">>> 3_enrichment running — GSEA + ORA with ratio labels")

## -------------------- Directories --------------------
project_root <- "/hpcnfs/data/BA/MYC_NMD_keep/mrezaei/gene_level"
data_dir     <- file.path(project_root, "data")

de_tbl   <- file.path(project_root, "results", "1_DEAnalysis", "tables", "shrunk")
int_tbl  <- file.path(project_root, "results", "2_interaction", "tables")

res_dir  <- file.path(project_root, "results", "3_enrichment")
tbl_dir  <- file.path(res_dir, "tables")
fig_dir  <- file.path(res_dir, "figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------- Parameters --------------------
padj_cutoff <- 0.05

## -------------------- Annotation --------------------
annotation <- read.delim(file.path(data_dir, "annotation.txt"), check.names = FALSE)
annotation <- annotation %>%
  dplyr::select(gene_id, symbol) %>%
  distinct(gene_id, .keep_all = TRUE)
id2sym <- setNames(annotation$symbol, annotation$gene_id)

## ============================================================
## MSigDB
## ============================================================
message("Loading MSigDB...")
msigdb <- list(
  HALLMARK = msigdbr("Mus musculus", "H"),
  GOBP     = msigdbr("Mus musculus", "C5", "BP"),
  C2CP     = msigdbr("Mus musculus", "C2", "CP"),
  KEGG     = msigdbr("Mus musculus", "C2", "CP:KEGG")
)
gmt_list <- lapply(msigdb, function(df)
  split(df$ensembl_gene, df$gs_name)
)

strip_version <- function(x) sub("\\..*$", "", x)

## ============================================================
## Helper plotting functions
## ============================================================

plot_gsea_bar <- function(gsea_df, out_file, title = NULL) {
  df <- gsea_df %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    arrange(NES) %>%
    slice_head(n = 20) %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES),
           ratio_label = paste0(leadingEdgeCount, "/", size))
  if (nrow(df) == 0) return(NULL)
  
  p <- ggplot(df, aes(x = NES, y = pathway, fill = NES > 0)) +
    geom_col() +
    geom_text(aes(x = NES / 2, label = ratio_label),
              color = "black", size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c("TRUE" = "#D62728", "FALSE" = "#1F77B4"),
                      labels = c("FALSE" = "Suppressed (NES<0)", "TRUE" = "Activated (NES>0)"),
                      name = "Direction") +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL, title = title) +
    theme_bw(9)
  ggsave(out_file, p, width = 7, height = 5, dpi = 300)
}

plot_ora_bar <- function(res, out_file, title = NULL, direction = "up") {
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(NULL)
  df <- df %>%
    arrange(p.adjust) %>%
    mutate(RatioLabel = paste0(Count, "/", setSize))
  color <- if (direction == "up") "#D62728" else "#1F77B4"
  
  p <- ggplot(df, aes(x = Count, y = reorder(Description, Count))) +
    geom_col(fill = color) +
    geom_text(aes(label = RatioLabel),
              position = position_stack(vjust = 0.5), size = 3.2, color = "black", fontface = "bold") +
    labs(x = "Gene count", y = NULL, title = title) +
    theme_bw(9) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(out_file, p, width = 7, height = 5, dpi = 300)
}

## ============================================================
## GSEA
## ============================================================

tt_files <- c(
  list.files(de_tbl,  pattern = "^tT_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^tT_.*\\.tsv$", full.names = TRUE)
)
if (length(tt_files) == 0) stop("No tT_*.tsv files found for GSEA.")

for (fp in tt_files) {
  tag_raw <- tools::file_path_sans_ext(basename(fp))
  tag <- sub("^tT_", "", tag_raw)
  cat("\n>>> GSEA on:", tag, "\n")
  
  df <- read.delim(fp, check.names = FALSE)
  if (!all(c("gene_id", "log2FoldChange") %in% names(df))) next
  
  df2 <- df %>% filter(!is.na(gene_id) & !is.na(log2FoldChange))
  if (nrow(df2) < 50) next
  ranks <- setNames(df2$log2FoldChange, df2$gene_id)
  
  for (set_tag in names(gmt_list)) {
    cat("  -", set_tag, "\n")
    gres <- tryCatch(
      fgseaMultilevel(pathways = gmt_list[[set_tag]], stats = ranks),
      error = function(e) NULL
    )
    if (is.null(gres) || nrow(gres) == 0) next
    
    gres <- gres %>%
      arrange(padj, desc(NES)) %>%
      mutate(leadingEdgeCount = vapply(leadingEdge, length, integer(1)),
             leadingEdge_symbols = vapply(
               leadingEdge,
               function(ids) paste(na.omit(id2sym[ids]), collapse = ";"),
               character(1)
             )) %>%
      mutate(across(where(is.list),
                    ~ vapply(.x, paste, collapse = ";", character(1))))
    
    out_tsv <- file.path(tbl_dir, paste0("GSEA_", set_tag, "_", tag, ".tsv"))
    write.table(gres, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    
    out_png <- file.path(fig_dir, paste0("GSEA_", set_tag, "_", tag, ".png"))
    plot_gsea_bar(gres, out_png, paste0("GSEA - ", set_tag, " (", tag, ")"))
  }
}

cat("\n>>> GSEA complete.\n")

## ============================================================
## ORA (fixed setSize + ratio labeling)
## ============================================================

up_files <- c(
  list.files(de_tbl,  pattern = "^up_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^up_.*\\.tsv$", full.names = TRUE)
)
down_files <- c(
  list.files(de_tbl,  pattern = "^down_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^down_.*\\.tsv$", full.names = TRUE)
)
all_files <- c(up_files, down_files)
if (length(all_files) == 0) stop("No up/down tables found for ORA.")

for (fp in all_files) {
  tag_raw <- tools::file_path_sans_ext(basename(fp))
  tag <- sub("\\.tsv$", "", tag_raw)
  cat("\n>>> ORA on:", tag, "\n")
  
  df <- read.delim(fp, check.names = FALSE)
  if (!"gene_id" %in% colnames(df)) next
  genes_ens <- unique(na.omit(strip_version(df$gene_id)))
  if (length(genes_ens) < 10) {
    cat("  [info] <10 genes; skipping.\n")
    next
  }
  
  direction <- if (grepl("^up_", tag)) "up" else "down"
  bar_color <- if (direction == "up") "#D62728" else "#1F77B4"
  
  for (set_name in names(msigdb)) {
    cat("  -", set_name, "\n")
    
    # TERM2GENE and set sizes
    term2gene <- msigdb[[set_name]] %>%
      dplyr::mutate(ensembl_gene = strip_version(ensembl_gene)) %>%
      dplyr::distinct(gs_name, ensembl_gene)
    
    term_sizes <- term2gene %>%
      dplyr::count(gs_name, name = "setSize")
    
    # ORA
    res <- tryCatch(
      clusterProfiler::enricher(
        gene = genes_ens,
        TERM2GENE = term2gene,
        pAdjustMethod = "BH",
        qvalueCutoff = padj_cutoff
      ),
      error = function(e) NULL
    )
    if (is.null(res)) next
    
    df_out <- as.data.frame(res)
    if (nrow(df_out) == 0) next
    
    # Add setSize by joining on ID (ID == gs_name)
    df_out <- df_out %>%
      dplyr::left_join(term_sizes, by = c("ID" = "gs_name"))
    
    # If some terms didn’t match (shouldn’t happen), fill safely
    df_out$setSize[is.na(df_out$setSize)] <- 0L
    
    # SYMBOLS column from gene IDs
    df_out$SYMBOLS <- vapply(
      strsplit(df_out$geneID, "/"),
      function(ids) paste(na.omit(id2sym[strip_version(ids)]), collapse = ";"),
      FUN.VALUE = character(1)
    )
    
    # Ratio columns
    df_out$Ratio      <- df_out$Count / pmax(df_out$setSize, 1L)
    df_out$RatioLabel <- paste0(df_out$Count, "/", df_out$setSize)
    
    # Save FULL table
    out_tsv <- file.path(tbl_dir, paste0("ORA_", set_name, "_", tag, ".tsv"))
    write.table(df_out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Plot ONLY significant
    df_sig <- df_out %>% dplyr::filter(p.adjust < padj_cutoff)
    if (nrow(df_sig) == 0) next
    
    df_plot <- df_sig %>%
      dplyr::arrange(p.adjust)
    
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Ratio, y = reorder(Description, Ratio))) +
      ggplot2::geom_col(fill = bar_color) +
      # Put label in the middle of each bar:
      ggplot2::geom_text(ggplot2::aes(x = Ratio / 2, label = RatioLabel),
                         size = 3.2, color = "black", fontface = "bold") +
      ggplot2::labs(x = "Gene ratio", y = NULL) +
      ggplot2::theme_bw(9) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    
    out_png <- file.path(fig_dir, paste0("ORA_", set_name, "_", tag, ".png"))
    ggplot2::ggsave(out_png, p, width = 7, height = 5, dpi = 300)
  }
}

cat("\n>>> ORA complete\nTables: ", tbl_dir, "\nFigures: ", fig_dir, "\n")

