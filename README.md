# BaseDepEx
A modular pipeline for Baseline-dependent expression analysis at both gene and isoform levels and downstream analyses.

This repository contains a modular RNA-seq analysis workflow for studying **drug responses under different baseline conditions**.

Currently included:
- `RScripts/1_DEAnalysis.R`: Differential expression analysis and visualization (MA, Volcano, PCA, dispersion)
- `data/`: Example input files (counts, metadata, annotation)
- `results/`: Output folders for tables and figures

---

## Usage

Run from the repository root:
```bash
Rscript RScripts/1_DEAnalysis.R