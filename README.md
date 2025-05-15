# GBA1-multiomics

Last Updated: May 2025

# Summary
The following repository encompasses all scripts developed and used in the manuscript titled "Dissecting the biological impact of GBA1 mutations using multi-omics in an isogenic setting". This study explores potential molecular mechanisms and downstream consequences of GCase reduction.

# Analyses

**Languages**: Bash, R, Python

Each analysis directory includes:

- A notebook that calls the scripts  
- A collection of scripts used in the notebook  

## Overview of Analysis Directories

| Analysis Directory       | Description |
|--------------------------|-------------|
| `00_raw_processing`             | 1) Process raw Illumina RNA data to get a final table of quantifications. 2) Clean up .adat file received from Somalogic for downstream analysis. |
| `01_PCA_scree`        | 1) Generate PCs, make PCA plots, and add PCs to analysis tables. 2) Generate scree plots of PCs to determine how many to keep for analysis. |
| `02_regressions`             | 1) Run regression with GBA1 genotypes vs outcomes (genes, proteins, metabolites). 2) Generate QQ plot to check lambda values. |
| `03_plots` | 1) Example scripts of how figures were plotted and made. |
