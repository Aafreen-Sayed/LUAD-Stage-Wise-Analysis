# Stage-Wise Evolutionary Transcriptomic Analysis of LUAD

## Project Overview
This repository contains the computational framework developed to analyze the **evolutionary progression of Lung Adenocarcinoma (LUAD)**. Unlike standard Tumor vs. Normal comparisons, this pipeline implements a **stage-wise contrast model** (Stage I vs II, II vs III, III vs IV) to isolate the molecular events driving clinical disease transitions.

## Pipeline Structure

### 1. Data Preprocessing (`01_Data_Preprocessing.R`)
- Standardizes clinical staging metadata (converting Roman numerals to numeric stages).
- Filters low-quality samples and normalizes microarray intensities.

### 2. Evolutionary DGE Analysis (`02_StageWise_DGE_Pipeline.R`)
- Implements a linear model using `limma`.
- Defines adjacent stage contrasts to identify phase-specific dysregulated genes.
- **Key Logic:** `makeContrasts(StageN - StageN-1)`

### 3. Functional Enrichment (`03_Functional_Enrichment.R`)
- Maps DEGs to biological pathways using `clusterProfiler`.
- Performs GO (Gene Ontology) and KEGG pathway enrichment analysis.
- Generates network visualizations (`cnetplot`) to link gene clusters to disease mechanisms.

## Tools Used
- **R**: v4.2+
- **Bioconductor**: limma, GEOquery, clusterProfiler, org.Hs.eg.db
- **Data Handling**: dplyr, openxlsx

---
*Note: The raw data files and manuscript findings are not included in this repository as the work is currently under review.*
