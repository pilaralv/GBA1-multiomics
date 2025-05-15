library(tidyverse)
library(optparse)
library(Biostrings)
library(tximport)
library(readr)
library(tximportData)
library(jsonlite)

arguments <- parse_args(OptionParser(), positional_arguments = 4)
gene_map<-arguments$args[1]
quants_salmon_path<-arguments$args[2]
out_name<-arguments$args[3]
sample<-arguments$args[4]

tx2gene <- read.table(gene_map) 
tx2gene<- tx2gene[, c(1,3,2,4)]
txi <- tximport(quants_salmon_path, type="salmon", tx2gene = tx2gene,abundanceCol = TPM)
df <- data.frame(txi)
df$id <- row.names(df)
colnames(df) <- c('TPM', 'COUNTS', 'LENGTH', 'countsFromAbundance', 'ID')
# Function to append suffix to a column name
append_suffix_to_column <- function(data, column_name, suffix) {
  new_column_name <- paste(column_name, suffix, sep = "_")
  names(data)[names(data) == column_name] <- new_column_name
  data
}

# Example usage
df <- append_suffix_to_column(df, "TPM", sample)
df <- append_suffix_to_column(df, "COUNTS", sample)
df <- append_suffix_to_column(df, "LENGTH", sample)


write_csv(df, out_name)


