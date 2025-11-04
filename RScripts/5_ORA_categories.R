## ============================================================
## 5_ORA_Categories.R — ORA on MYC categories
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Mm.eg.db)
})

set.seed(1)
message(">>> 5_ORA_Categories running — ORA on MYC category genes")

## -------------------- Directories --------------------
project_root <- "/hpcnfs/data/BA/MYC_NMD_keep/mrezaei/gene_level"
cat_tbl <- file.path(project_root, "results/4_categories/tables")
ora_tbl <- file.path(project_root, "results/5_ORA_categories/tables")
ora_fig <- file.path(project_root, "results/5_ORA_categories/figures")
dir.create(ora_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(ora_fig, recursive = TRUE, showWarnings = FALSE)

## -------------------- Parameters --------------------
padj_cutoff <- 0.05

## -------------------- Load MSigDB gene sets (ENSEMBL IDs) --------------------
message("Loading MSigDB ...")
msigdb <- list(
  HALLMARK = msigdbr("Mus musculus", "H"),
  GOBP     = msigdbr("Mus musculus", "C5", "BP"),
  C2CP     = msigdbr("Mus musculus", "C2", "CP"),
  KEGG     = msigdbr("Mus musculus", "C2", "CP:KEGG")
)

gene_sets <- lapply(msigdb, function(df) {
  df %>% dplyr::select(gs_name, ensembl_gene) %>% distinct()
})

set_sizes <- lapply(gene_sets, function(df) {
  df %>% group_by(gs_name) %>% summarise(setSize = n(), .groups = "drop")
})

strip_version <- function(x) sub("\\..*$", "", x)

## -------------------- Helper for ORA plot --------------------
plot_ora_bar <- function(res, term_sizes_df, out_file, title_lab) {
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(invisible(NULL))
  df <- df %>%
    left_join(term_sizes_df, by = c("ID" = "gs_name")) %>%
    mutate(
      setSize = ifelse(is.na(setSize), 0L, setSize),
      RatioLabel = paste0(Count, "/", setSize)
    ) %>%
    arrange(p.adjust) %>%
    slice_head(n = 20)
  p <- ggplot(df, aes(x = Count, y = reorder(Description, Count), fill = -log10(p.adjust))) +
    geom_col() +
    geom_text(aes(label = RatioLabel),
              position = position_stack(vjust = 0.5),
              size = 3, color = "black") +
    scale_fill_gradientn(colors = c("blue", "white", "red"), name = "-log10(p.adjust)") +
    labs(x = "Count", y = NULL, title = title_lab) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(out_file, p, width = 7, height = 5, dpi = 300)
}

## ============================================================
## Run ORA for all MYC categories
## ============================================================

cat_files <- list.files(cat_tbl, pattern = "\\.tsv$", full.names = TRUE)
if (length(cat_files) == 0) stop("No category files found in ", cat_tbl)

for (fp in cat_files) {
  tag <- tools::file_path_sans_ext(basename(fp))
  message("\n>>> ORA on category: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  
  ## Detect ID column (prefer gene_id if present)
  id_col <- if ("gene_id" %in% names(df)) "gene_id" else if ("feature_id" %in% names(df)) "feature_id" else if ("isoform_id" %in% names(df)) "isoform_id" else NULL
  if (is.null(id_col)) {
    message("  [skip] No suitable ID column in ", tag)
    next
  }
  
  ## Extract unique Ensembl IDs
  genes_ens <- unique(na.omit(strip_version(df[[id_col]])))
  if (length(genes_ens) < 5) {
    message("  [info] <5 genes; skipping.")
    next
  }
  
  ## Run ORA for each gene set
  for (set_name in names(gene_sets)) {
    message("  - ", set_name)
    term2gene <- gene_sets[[set_name]] %>%
      mutate(ensembl_gene = strip_version(ensembl_gene)) %>%
      distinct(gs_name, ensembl_gene)
    
    res <- tryCatch(
      enricher(
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
    
    ## Add SYMBOLS for clarity
    ENTREZ2SYMBOL <- AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = unique(unlist(strsplit(df_out$geneID, "/"))),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    df_out$SYMBOLS <- vapply(
      strsplit(df_out$geneID, "/"),
      function(ids) paste(na.omit(ENTREZ2SYMBOL[ids]), collapse = ";"),
      FUN.VALUE = character(1)
    )
    
    ## Save ORA table
    out_tsv <- file.path(ora_tbl, paste0("ORA_", set_name, "_", tag, ".tsv"))
    write.table(df_out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    
    ## Plot ORA (significant terms)
    plot_ora_bar(
      res,
      term_sizes_df = set_sizes[[set_name]],
      out_file = file.path(ora_fig, paste0("ORA_", set_name, "_", tag, ".png")),
      title_lab = paste(set_name, "-", tag)
    )
  }
}

cat("\n>>> ORA complete.\nTables → ", ora_tbl, "\nFigures → ", ora_fig, "\n")

