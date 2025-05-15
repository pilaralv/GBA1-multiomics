install.packages("SomaDataIO")
library(SomaDataIO)
library(dplyr)
library(stringr)

f <- "CHI-23-015_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP.adat"

my_adat <- read_adat(f)

# Replace seqid with gene name

# This file has only seq ids and gene names matching with human, all other organisms were manually removed
genes <- read.csv("seq_id_genes.csv", header=F)
genes <- genes %>%
  mutate(V1 = str_replace_all(paste0("seq.", V1), "-", "."))
# Clean up gene names
genes$V2 <- gsub("\\|", "_", genes$V2)

# Get the column names from seq$V1
columns_to_keep <- genes$V1

# Subset f to keep only columns present in columns_to_keep
subset_adat <- my_adat[, intersect(names(my_adat), columns_to_keep)]

final_subset <- cbind(my_adat[, 1:32], subset_adat)

# Replace seqs with gene name
# Iterate over the column names of final_subset
for (col_name in colnames(final_subset)) {
  # Check if the column name matches any entry in the V1 column of genes
  if (col_name %in% genes$V1) {
    # Find the corresponding V2 entry in genes
    new_col_name <- genes$V2[genes$V1 == col_name]
    # Replace the column name in final_subset
    colnames(final_subset)[colnames(final_subset) == col_name] <- new_col_name
  }
}

write.csv(final_subset, "somalogic_table_clean.csv", row.names=F)

# MANUALLY IN EXCEL UPDATE GENE NAMES AND REMOVED NON ISOGENIC LINES (NO FIBROBLASTS, LCLS, OR CALIBRATORS)
