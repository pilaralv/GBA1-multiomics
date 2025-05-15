# GBA1-multiomics
Analysis of multiomic data for GBA1 related disease mechanisms

# Analyses

**Languages**: Bash, R, Python

Each analysis directory includes:

- A notebook that calls the scripts  
- A collection of scripts used in the notebook  

## Overview of Analysis Directories

| Analysis Directory       | Description |
|--------------------------|-------------|
| `00_ONT_RNA`             | 1) Process raw Illumina RNA data to get a final table of quantifications. 2) Clean up .adat file received from Somalogic for downstream analysis. |
| `01_Illumina_RNA`        | 1) Generate PCs, make PCA plots, and add PCs to analysis tables. 2) Generate scree plots of PCs to determine how many to keep for analysis. |
| `02_ONT_DNA`             | 1) Run regression with GBA1 genotypes vs outcomes (genes, proteins, metabolites). 2) Generate QQ plot to check lambda values. |
| `03_Branchpoint_predictions` | 1) Example scripts of how figures were plotted and made. |
