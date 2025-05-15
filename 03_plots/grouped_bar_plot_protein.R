# Grouped bar plot for protein expression --> you can do smae thing for RNA
setwd("/data/tables/")
# Load necessary library
library(dplyr)
library(ggplot2)

# Read in the files
ips_data <- read.csv("analysis_table_ips_protein.csv")
microglia_data <- read.csv("analysis_table_mgl_protein_switch.csv")
neuron_data <- read.csv("analysis_table_neuron_protein.csv")

# SWITCH TO 409 variant ID
ips_data$Group <- gsub("D448", "D409", ips_data$Group)
microglia_data$Group <- gsub("D448", "D409", microglia_data$Group)
neuron_data$Group <- gsub("D448", "D409", neuron_data$Group)

# Add a Cell_Type column to each dataset
ips_data$Cell_Type <- "IPSC"
microglia_data$Cell_Type <- "Microglia"
neuron_data$Cell_Type <- "Neuron"

# Combine all datasets
combined_data <- bind_rows(ips_data, microglia_data, neuron_data)



# Ensure the Group variable is treated as a factor and ordered as desired
combined_data$Group <- factor(combined_data$Group, levels = c("WT", "D409H_HET", "D409V_HET", "D409H_HOM", "D409V_HOM", "KO"))

# Create a new grouping factor for the x-axis
combined_data <- combined_data %>%
  mutate(Grouping = factor(paste(Cell_Type, Group, sep = "_"), 
                           levels = c(
                             paste("IPSC", c("WT", "D409H_HET", "D409V_HET", "D409H_HOM", "D409V_HOM", "KO"), sep = "_"),
                             paste("Neuron", c("WT", "D409H_HET", "D409V_HET", "D409H_HOM", "D409V_HOM", "KO"), sep = "_"),
                             paste("Microglia", c("WT", "D409H_HET", "D409V_HET", "D409H_HOM", "D409V_HOM", "KO"), sep = "_")
                           )))

# Update x-axis labels to display only mutation names
combined_data <- combined_data %>%
  mutate(Label = Group)

# Ensure the Cell_Type variable is a factor and ordered as desired
combined_data$Cell_Type <- factor(combined_data$Cell_Type, levels = c("IPSC", "Neuron", "Microglia"))

# Define custom colors for the cell types
custom_colors <- c("IPSC" = "#E41A1C", "Neuron" = "#377EB8", "Microglia" = "#4DAF4A")  # Red, Blue, Green

# Plot grouped boxplot with aligned individual points and ordered facets
plot <- ggplot(combined_data, aes(x = Group, y = GBA, fill = Cell_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +  # Add boxplot
  geom_point(aes(color = Cell_Type), 
             position = position_dodge(width = 0.75), size = 1.5, alpha = 0.8) +  # Add aligned points
  facet_wrap(~Cell_Type, scales = "free_x", strip.position = "top") +
  labs(
    title = "GCase Protein Expression",
    x = NULL,
    y = "GCase Expression"
  ) +
  scale_fill_manual(values = custom_colors) +  # Custom fill colors for boxplots
  scale_color_manual(values = custom_colors) +  # Custom colors for points
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=18),
    axis.text.y = element_text(size=16),
    axis.title.y = element_text(size=16),
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add borders around facets
    legend.position = "none"  # Remove the legend
  )

plot
ggsave("grouped_boxplot_protein.png", plot = plot, width = 14, height = 10, dpi = 600)
