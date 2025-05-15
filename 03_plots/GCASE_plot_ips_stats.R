library(ggplot2)
library(ggpubr)

# GCASE data values
group <- c(
  rep("KO", 3),
  rep("WT", 3),
  rep("D409H_HET", 3),
  rep("D409V_HET", 3),
  rep("D409V_HOM", 3),
  rep("D409H_HOM", 3)
)

gcase <- c(
  221.137, 253.7345, 181.20825,             # KO
  9204.766, 8846.687, 8672.02575,           # WT
  5557.3865, 5226.4265, 6838.7605,          # D409H_HET
  4691.223, 5754.28875, 5275.0845,          # D409V_HET
  480.6715, 352.863, 471.531,               # D409V_HOM
  706.091, 730.8755, 727.382                # D409H_HOM
)

df_raw <- data.frame(group, gcase)

# ANOVA + Tukey
anova_res <- aov(gcase ~ group, data = df_raw)
tukey_res <- TukeyHSD(anova_res)
print(tukey_res)


# Extract significant comparisons
sig_results <- as.data.frame(tukey_res$group)
sig_results$comparison <- rownames(sig_results)
sig_results <- subset(sig_results, `p adj` < 0.05)

# Convert comparison string to list of pairs
comparisons <- strsplit(sig_results$comparison, "-")

library(ggpubr)

# Build significance label dataframe
sig_labels <- data.frame(
  group1 = sapply(comparisons, `[`, 1),
  group2 = sapply(comparisons, `[`, 2),
  y.position = max(data$GCase + data$SD) * seq(1.05, 1.2, length.out = length(comparisons)),  # stack labels
  p.adj = sig_results$`p adj`,
  p.label = ifelse(sig_results$`p adj` < 0.001, "***",
                   ifelse(sig_results$`p adj` < 0.01, "**",
                          ifelse(sig_results$`p adj` < 0.05, "*", "ns")))
)

# Add to your existing plot
plot +
  stat_pvalue_manual(sig_labels, label = "p.label", tip.length = 0.01)


# Make 'group' a factor in desired order
df_raw$group <- factor(df_raw$group, levels = c("WT", "D409H_HET", "D409V_HET", "D409H_HOM", "D409V_HOM", "KO"))

# Define comparisons
comparisons <- list(
  c("WT", "KO"),
  c("WT", "D409V_HOM"),
  c("WT", "D409H_HOM"),
  c("WT", "D409V_HET"),
  c("WT", "D409H_HET")
)

# Plot
plot <- ggplot(df_raw, aes(x = group, y = gcase, fill = group)) +
  geom_boxplot(fill = "#E41A1C",width = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0, size = 2, alpha = 0.6) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") +
  scale_y_continuous(breaks = seq(0, 10000, by = 2000)) +
  labs(title = "GCase Activity in IPSC Samples", y = "GCase Activity", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y= element_text(size=14),
    axis.title.y= element_text(size=14),
    plot.title = element_text(hjust = 0.5, face = "bold", size=18),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

plot

# Save as PNG with transparent background
ggsave("GCase_barplot_ips_stats.png", plot = plot, width = 10, height = 6, dpi = 600, bg = "transparent")
