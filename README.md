# BaseDepEx
A modular pipeline for Baseline-dependent expression analysis at both gene and isoform levels and downstream analyses.

This repository contains a modular RNA-seq analysis workflow for studying **drug responses under different baseline conditions**.

Currently included:
- `RScripts/1_DEAnalysis.R`: Differential expression analysis and visualization (MA, Volcano, PCA, dispersion)
- `RScripts/2_interaction.R`: Multifactorial DESeq2 for interaction of treatments with baseline

- `data/`: Example input files (gene_count.tsv, isoform_count.tsv, metadata.txt, annotation.txt)

- `results/`: Output folders for tables and figures

---

## Usage

Use RStudio or run from the repository root:
```bash
Rscript RScripts/1_DEAnalysis.R