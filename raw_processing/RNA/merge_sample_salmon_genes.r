# Load the required package
library(dplyr)



# Function to merge files by "ID" column
merge_files_by_id <- function(file_names) {
  # Initialize an empty dataframe to store the merged data
  merged_data <- data.frame()
  
  # Loop through the file names and merge data
  for (file_name in file_names) {
    # Read the data from the file
    data <- read.csv(file_name)  # Change read.csv to appropriate function based on file format
    
    # Merge data with the existing merged_data using "ID" column
    if (nrow(merged_data) == 0) {
      # If it's the first file, directly assign it to merged_data
      merged_data <- data
    } else {
      # If not the first file, perform an inner join based on "ID" column
      merged_data <- inner_join(merged_data, data, by = "ID")
    }
  }
  
  return(merged_data)
}




# Example usage:
file_names <- c("./KOLF_IPSC_GBA1_D448H_HET1/KOLF_IPSC_GBA1_D448H_HET1.genes.csv", "./KOLF_IPSC_GBA1_D448H_HET2/KOLF_IPSC_GBA1_D448H_HET2.genes.csv", "./KOLF_IPSC_GBA1_D448H_HET3/KOLF_IPSC_GBA1_D448H_HET3.genes.csv", "./KOLF_IPSC_GBA1_D448H_HOM1/KOLF_IPSC_GBA1_D448H_HOM1.genes.csv", "./KOLF_IPSC_GBA1_D448H_HOM2/KOLF_IPSC_GBA1_D448H_HOM2.genes.csv", "./KOLF_IPSC_GBA1_D448H_HOM3/KOLF_IPSC_GBA1_D448H_HOM3.genes.csv", "./KOLF_IPSC_GBA1_KO1/KOLF_IPSC_GBA1_KO1.genes.csv", "./KOLF_IPSC_GBA1_KO2/KOLF_IPSC_GBA1_KO2.genes.csv", "./KOLF_IPSC_GBA1_KO3/KOLF_IPSC_GBA1_KO3.genes.csv", "./KOLF_IPSC_GBA1_D448V_HOM1/KOLF_IPSC_GBA1_D448V_HOM1.genes.csv", "./KOLF_IPSC_GBA1_D448V_HOM2/KOLF_IPSC_GBA1_D448V_HOM2.genes.csv", "./KOLF_IPSC_GBA1_D448V_HOM3/KOLF_IPSC_GBA1_D448V_HOM3.genes.csv", "./KOLF_IPSC_GBA1_D448V_HET1/KOLF_IPSC_GBA1_D448V_HET1.genes.csv", "./KOLF_IPSC_GBA1_D448V_HET2/KOLF_IPSC_GBA1_D448V_HET2.genes.csv", "./KOLF_IPSC_GBA1_D448V_HET3/KOLF_IPSC_GBA1_D448V_HET3.genes.csv", "./KOLF_IPSC_WT_rep1/KOLF_IPSC_WT_rep1.genes.csv", "./KOLF_IPSC_WT_rep2/KOLF_IPSC_WT_rep2.genes.csv", "./KOLF_IPSC_WT_rep3/KOLF_IPSC_WT_rep3.genes.csv")
merged_data <- merge_files_by_id(file_names)

# Keep only IDS + TPM
merged_data<- merged_data[, c(5,1,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70)]


write.table(merged_data, "all_samples_salmon_genes_ips.txt", row.names=FALSE)
write.csv(merged_data, "all_samples_salmon_genes_ips.csv", row.names=FALSE)
