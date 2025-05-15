# Setwd
setwd("/data/metabolomics_figure/")

# Load packages
library(ComplexHeatmap)
library(dplyr)  # For data manipulation
library(circlize)  # For colorRamp2 function
library(grid)  # For adding the title with grid.text()
library(purrr)

# Read in the three datasets
ips <- read.csv("metabolomics_ips_ceramides.csv", sep=",", header=TRUE, check.names=FALSE)
neuron <- read.csv("metabolomics_neuron_ceramides.csv", sep=",", header=TRUE, check.names=FALSE)
mgl <- read.csv("metabolomics_microglia_ceramides.csv", sep=",", header=TRUE, check.names=FALSE)

# Rename sample columns to include the dataset suffix
colnames(ips)[4:ncol(ips)] <- paste0(colnames(ips)[4:ncol(ips)], "_IPS")
colnames(neuron)[4:ncol(neuron)] <- paste0(colnames(neuron)[4:ncol(neuron)], "_NEURON")
colnames(mgl)[4:ncol(mgl)] <- paste0(colnames(mgl)[4:ncol(mgl)], "_MGL")

# Merge datasets using a full join to keep all metabolites
merged_data <- reduce(list(ips, neuron, mgl), 
                      function(x, y) full_join(x, y, by = c("Metabolite", "Metabolite_Full", "Sub_group")))

# Replace all NA values with 0
merged_data[is.na(merged_data)] <- 0

# Group by Metabolite to remove duplicates by summing values across duplicate rows
merged_data <- merged_data %>%
  group_by(Metabolite, Metabolite_Full, Sub_group) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")  #


# Save the merged dataset to a CSV file
write.csv(merged_data, "joined_metabolomics_ceramide_data.csv", row.names=FALSE)


### Create heatmap

# Extracting the part of the data for heatmap (numeric values and rownames as Metabolite_Full)
heatmap_data <- merged_data[, 4:ncol(merged_data)] # Select only numeric columns
rownames(heatmap_data) <- merged_data$Metabolite_Full

# Calculate averages for the replicates
#IPS
heatmap_data$WT_IPS <- rowMeans(heatmap_data[, c("WT1_IPS", "WT2_IPS", "WT3_IPS")], na.rm = TRUE)
heatmap_data$D409H_HET_IPS <- rowMeans(heatmap_data[, c("D448H_HET1_IPS", "D448H_HET2_IPS", "D448H_HET3_IPS")], na.rm = TRUE)
heatmap_data$D409V_HET_IPS <- rowMeans(heatmap_data[, c("D448V_HET1_IPS", "D448V_HET2_IPS", "D448V_HET3_IPS")], na.rm = TRUE)
heatmap_data$D409H_HOM_IPS <- rowMeans(heatmap_data[, c("D448H_HOM1_IPS", "D448H_HOM2_IPS", "D448H_HOM3_IPS")], na.rm = TRUE)
heatmap_data$D409V_HOM_IPS <- rowMeans(heatmap_data[, c("D448V_HOM1_IPS", "D448V_HOM2_IPS", "D448V_HOM3_IPS")], na.rm = TRUE)
heatmap_data$KO_IPS <- rowMeans(heatmap_data[, c("KO1_IPS", "KO2_IPS", "KO3_IPS")], na.rm = TRUE)
# Neuron
heatmap_data$WT_NEURON <- rowMeans(heatmap_data[, c("WT1_NEURON", "WT2_NEURON", "WT3_NEURON")], na.rm = TRUE)
heatmap_data$D409H_HET_NEURON <- rowMeans(heatmap_data[, c("D448H_HET1_NEURON", "D448H_HET2_NEURON", "D448H_HET3_NEURON")], na.rm = TRUE)
heatmap_data$D409V_HET_NEURON <- rowMeans(heatmap_data[, c("D448V_HET1_NEURON", "D448V_HET2_NEURON", "D448V_HET3_NEURON")], na.rm = TRUE)
heatmap_data$D409H_HOM_NEURON <- rowMeans(heatmap_data[, c("D448H_HOM1_NEURON", "D448H_HOM2_NEURON", "D448H_HOM3_NEURON")], na.rm = TRUE)
heatmap_data$D409V_HOM_NEURON <- rowMeans(heatmap_data[, c("D448V_HOM1_NEURON", "D448V_HOM2_NEURON", "D448V_HOM3_NEURON")], na.rm = TRUE)
heatmap_data$KO_NEURON <- rowMeans(heatmap_data[, c("KO1_NEURON", "KO2_NEURON", "KO3_NEURON")], na.rm = TRUE)
# MGL
heatmap_data$WT_MGL <- rowMeans(heatmap_data[, c("WT1_MGL", "WT2_MGL", "WT3_MGL")], na.rm = TRUE)
heatmap_data$D409H_HET_MGL <- rowMeans(heatmap_data[, c("D448H_HET1_MGL", "D448H_HET2_MGL", "D448H_HET3_MGL")], na.rm = TRUE)
heatmap_data$D409V_HET_MGL <- rowMeans(heatmap_data[, c("D448V_HET1_MGL", "D448V_HET2_MGL", "D448V_HET3_MGL")], na.rm = TRUE)
heatmap_data$D409H_HOM_MGL <- rowMeans(heatmap_data[, c("D448H_HOM1_MGL", "D448H_HOM2_MGL", "D448H_HOM3_MGL")], na.rm = TRUE)
heatmap_data$D409V_HOM_MGL <- rowMeans(heatmap_data[, c("D448V_HOM1_MGL", "D448V_HOM2_MGL", "D448V_HOM3_MGL")], na.rm = TRUE)
heatmap_data$KO_MGL <- rowMeans(heatmap_data[, c("KO1_MGL", "KO2_MGL", "KO3_MGL")], na.rm = TRUE)


# Add back the Metabolite_Full and Sub_group to the heatmap data
heatmap_data$Metabolite_Full <- merged_data$Metabolite_Full
heatmap_data$Sub_group <- merged_data$Sub_group

# First, reorder the data by Sub_group and then alphabetically by Metabolite_Full
heatmap_data <- heatmap_data %>%
  arrange(Sub_group, Metabolite_Full)

# Ensure columns are explicitly ordered after sorting the rows
averaged_heatmap_data <- heatmap_data[, c("WT_IPS", "D409H_HET_IPS", "D409V_HET_IPS", "D409H_HOM_IPS", "D409V_HOM_IPS", "KO_IPS",
                                          "WT_NEURON", "D409H_HET_NEURON", "D409V_HET_NEURON", "D409H_HOM_NEURON", "D409V_HOM_NEURON", "KO_NEURON",
                                          "WT_MGL", "D409H_HET_MGL", "D409V_HET_MGL", "D409H_HOM_MGL", "D409V_HOM_MGL", "KO_MGL")]
rownames(averaged_heatmap_data) <- heatmap_data$Metabolite_Full

# Apply log2 transformation to the data (adding a small constant to avoid log(0) issues)
log_transformed_data <- log2(averaged_heatmap_data + 1)

# Create a custom color ramp (from white to red)
color_fun <- colorRamp2(c(0, 5), c("white", "red"))  # Fixed to range from 0 to 5


# Create a column annotation for the subgroups
subgroup_annotation <- rowAnnotation(
  `Metabolite Pathway` = heatmap_data$Sub_group,  # Renamed the key
  col = list(`Metabolite Pathway` = c("Ceramides" = "#ABDDA4",  # Peach
                                      "Hexosylceramides (HCER)" = "#3288BD",  # Light Coral
                                      "Dihydroceramides" = "#66C2A5",  # Peach Orange
                                      "Lactosylceramides (LCER)" = "#5E4FA2",  # Light Red
                                      "Ceramide PEs" = "#E6F598")),
  show_annotation_name = FALSE  # Keep this to show the new label
)

# Define column split points (index of last column before transition)
ko_ips_index <- which(colnames(log_transformed_data) == "KO_IPS")
ko_neuron_index <- which(colnames(log_transformed_data) == "KO_NEURON")

# Get the total number of columns
num_cols <- ncol(log_transformed_data)

# Create the heatmap
ht <- Heatmap(as.matrix(log_transformed_data),  
              name = "Expression (log2)", 
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 10),
              col = color_fun,
              right_annotation = subgroup_annotation,
              cluster_columns = FALSE,
              cluster_rows = FALSE)

# Save the heatmap as a PNG file
png("heatmap_all_joined_409.png", width = 6000, height = 5000, res = 600)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", padding = unit(c(10, 10, 10, 10), "mm"))


# Add vertical lines at the right edge of KO_IPS and KO_NEURON
decorate_heatmap_body("Expression (log2)", {
  grid.lines(x = unit((ko_ips_index + 0) / num_cols, "npc"),  # Shift line to the right edge
             y = unit(c(0, 1), "npc"), 
             gp = gpar(lwd = 2, col = "black"))
  
  grid.lines(x = unit((ko_neuron_index + 0) / num_cols, "npc"),  # Shift line to the right edge
             y = unit(c(0, 1), "npc"), 
             gp = gpar(lwd = 2, col = "black"))
})


# Add title
grid.text("All Ceramides", 
          x = unit(0.13, "npc"), 
          y = unit(1, "npc") - unit(5, "mm"), 
          gp = gpar(fontsize = 12, fontface = "bold"))

dev.off()