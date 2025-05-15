# Check RNA marker genes --> can do the same thing for protein

library(tidyr)

# Read in the files
ips_data <- read.csv("marker_gene_table_ips.csv")
microglia_data <- read.csv("marker_gene_table_mgl.csv")
neuron_data <- read.csv("marker_gene_table_neuron.csv")

# Add a Cell_Type column to each dataset
ips_data$Cell_Type <- "IPSC"
microglia_data$Cell_Type <- "Microglia"
neuron_data$Cell_Type <- "Neuron"

# Combine all datasets
combined_data <- bind_rows(ips_data, microglia_data, neuron_data)


# Subset the data to include only the desired genes
genes_of_interest <- c("POU5F1", "MAP2", "TREM2")

long_data <- combined_data %>%
  select(Cell_Type, Group, POU5F1, MAP2, TREM2) %>%
  pivot_longer(
    cols = c(POU5F1, MAP2, TREM2),
    names_to = "Gene",
    values_to = "Expression"
  )

# Ensure Gene is a factor and ordered as desired
long_data$Gene <- factor(long_data$Gene, levels = c("POU5F1", "MAP2", "TREM2"))

# Ensure Cell_Type is ordered as IPS, Neuron, Microglia
long_data$Cell_Type <- factor(long_data$Cell_Type, levels = c("IPSC", "Neuron", "Microglia"))

# Define custom colors for the cell types
custom_colors <- c("IPSC" = "#E41A1C", "Neuron" = "#377EB8", "Microglia" = "#4DAF4A")  # Red, Blue, Green

# Plot grouped violin plots with consistent y-axis and ordered x-axis
plot <- ggplot(long_data, aes(x = Cell_Type, y = Expression, fill = Cell_Type)) +
  geom_violin(alpha = 0.7, scale = "width", trim = FALSE) +  # Add violin plot
  geom_point(
    aes(color = Cell_Type), 
    position = position_dodge(width = 0.2), size = 1.5, alpha = 0.8  # Add jittered points
  ) +
  facet_wrap(~Gene, scales = "fixed", strip.position = "top") +  # Use fixed scales for consistent y-axis
  labs(
    title = "Gene Expression Across Cell Types (Violin Plot)",
    x = NULL,
    y = "Expression (TPM)"
  ) +
  scale_fill_manual(values = custom_colors) +  # Custom fill colors for violin plots
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

# Display the plot
print(plot)
ggsave("marker_violin_rna.png", plot = plot, width = 10, height = 8, dpi = 600)
