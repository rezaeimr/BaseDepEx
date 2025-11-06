# **BaseDepEx: Baseline-Dependent Expression Analysis Pipeline**

**BaseDepEx** is a modular R-based workflow designed for *baseline-dependent differential expression analysis*.  
It systematically quantifies how drug treatments alter transcription under different baseline conditions (e.g., OHT OFF vs OHT ON) and identifies genes whose responses depend on MYC activation or other regulatory states.

The pipeline supports both **gene-level** and **isoform-level** analyses, performs full DESeq2-based differential expression, interaction modeling, and enrichment analyses, and classifies genes into mechanistic categories reflecting **MYC-enhanced**, **MYC-suppressed**, and **switched** responses.

---

## **Pipeline Overview**

The workflow consists of six modular steps, each implemented as a standalone R script located in the `RScripts/` directory.

| Step | Script | Description |
|------|---------|-------------|
| **1️- Differential Expression** | **1_DEAnalysis.R** | Performs DE analysis using **DESeq2** for all drug treatments at both baseline (OHT OFF) and activated (OHT ON) states. Supports both gene and isoform count matrices. Produces log₂FC tables, MA plots, volcano plots, PCA, and dispersion plots. |
| **2️- Interaction Analysis** | **2_interaction.R** | Tests for *drug × baseline (OHT)* interactions to identify genes with baseline-dependent drug responses. Outputs log₂FC tables and MA/volcano plots for interaction terms. |
| **3️- Enrichment Analysis** | **3_enrichment.R** | Performs both **GSEA** and **ORA** using MSigDB gene sets (Hallmark, GO BP, KEGG, and C2:CP). Generates tables and barplots with ratio labels, distinguishing activated (red, NES > 0) and suppressed (blue, NES < 0) pathways. |
| **4️- Category Detection** | **4_categories.R** | Integrates DE and interaction results to classify genes into **MYC-enhanced**, **MYC-suppressed**, and **Switched** regulatory categories. Visualizes overlaps using **UpSet plots**. |
| **5️- ORA on MYC Categories** | **5_ORA_Categories.R** | Performs ORA separately for each MYC regulatory category to identify enriched pathways specific to each functional group. Produces both tables and publication-quality barplots. |
| **6️- Scatterplots** | **6_scatterplots.R** | Compares log₂FC values between baseline and activated conditions, highlighting MYC-enhanced, MYC-suppressed, and switched gene categories across treatments. |

---

## **Core Concepts**

- **Baseline dependence** — measures how the transcriptional response to a drug changes under baseline activation (e.g., OHT treatment).  
- **MYC regulatory categorization** — genes are grouped into:
  - *MYC-enhanced* (amplified under activation),
  - *MYC-suppressed* (weakened under activation),
  - *Switched* (direction of regulation changes).
- **Gene- and isoform-level flexibility** — define `analysis_level = "gene"` or `"isoform"` at the top of each script.  
- **Modular reproducibility** — each step produces self-contained results that seamlessly feed into downstream steps.

---

## **MYC Category Definitions**

| Category | Definition |
|-----------|-------------|
| **MYC Enhanced Up** | Upregulated in drug (OHT OFF), up in drug + OHT, and up in interaction; or up in +OHT and interaction but not DEG in OHT OFF. |
| **MYC Enhanced Down** | Down in drug (OHT OFF), down in drug + OHT, and down in interaction; or down in +OHT and interaction but not DEG in OHT OFF. |
| **MYC Suppressed Up** | Up in drug (OHT OFF), up in drug + OHT, and *down* in interaction; or up in OHT OFF and down in interaction but not DEG in +OHT. |
| **MYC Suppressed Down** | Down in drug (OHT OFF), down in drug + OHT, and *up* in interaction; or down in OHT OFF and up in interaction but not DEG in +OHT. |
| **Switched (Increase)** | Down in OHT OFF, up in OHT ON, and up in interaction. |
| **Switched (Decrease)** | Up in OHT OFF, down in OHT ON, and down in interaction. |

---

## **Directory Structure**

```
├── data
│   ├── annotation.txt
│   └── metadata_NMD.txt
├── LICENSE
├── README.md
├── requirements.txt
├── results
│   ├── 1_DEAnalysis
│   │   └── README.md
│   └── 2_interaction
│       └── README.md
└── RScripts
    ├── 1_DEAnalysis.R
    ├── 2_interaction.R
    ├── 3_enrichment.R
    ├── 4_categories.R
    ├── 5_ORA_categories.R
    └── 6_scatterplot.R
```

## Input Data Requirements

All input files used in the BaseDepEx pipeline must be placed inside the **`data/`** directory.  
This directory contains raw count matrices, metadata, and annotation tables necessary for both **gene-level** and **isoform-level** differential expression analysis. Examples of Metadata and annotation files are provided in the **`data/`** directory.

### Required Files

| File Name | Description | Required Columns | Notes |
|------------|--------------|------------------|--------|
| **`metadata_NMD.txt`** | Experimental design metadata describing each sample. | `SampleID`, `drug`, `OHT`, `group`, `replicate` | The `group` column defines DESeq2 groups (e.g., `11j`, `11j_OHT`, `DMSO`, etc.). <br> The `drug` and `OHT` columns are required for interaction analyses (`2_interaction.R`). |
| **`gene_counts.tsv`** | Gene-level count matrix (one row per gene). | `gene_id`, followed by all sample names (matching `SampleID` in metadata). | Used when `analysis_level = "gene"`. |
| **`isoform_counts.tsv`** | Isoform-level count matrix (one row per transcript/isoform). | `isoform_id`, `gene_id`, followed by all sample names. | Used when `analysis_level = "isoform"`. |
| **`annotation.txt`** | Reference annotation table mapping Ensembl IDs to symbols. | For **gene-level:** `gene_id`, `symbol`.<br> For **isoform-level:** `isoform_id`, `gene_id`, `symbol`. | Must match the IDs used in count tables (e.g., `ENSMUSG...` for mouse). |



