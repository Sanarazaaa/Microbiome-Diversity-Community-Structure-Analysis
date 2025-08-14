# --- Install and load libraries ---
install.packages(c("tidyverse", "ggplot2", "reshape2", "pheatmap", "vegan", "patchwork", "grid", "png"))
library(tidyverse)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(vegan)
library(patchwork)
library(grid)
library(png)

# --- 1. Load data ---
file_path <- file.choose()
microbiome_data <- read_csv(file_path)
microbiome_data[,-1] <- sapply(microbiome_data[,-1], as.numeric)

# --- 2. Prepare data ---
count_matrix <- microbiome_data %>%
  column_to_rownames(var = "Taxonomy") %>%
  as.matrix()

# --- 3. Shannon diversity ---
shannon_values <- diversity(t(count_matrix), index = "shannon")
shannon_df <- data.frame(Sample = names(shannon_values),
                         Shannon_Diversity = shannon_values)

# --- 4. Relative abundance ---
rel_abund <- microbiome_data
rel_abund[,-1] <- apply(microbiome_data[,-1], 2, function(x) x / sum(x, na.rm = TRUE))
sample_order <- shannon_df %>%
  arrange(desc(Shannon_Diversity)) %>%
  pull(Sample)
rel_abund <- rel_abund %>%
  select(Taxonomy, all_of(sample_order))

# Top 20 taxa
taxa_totals <- rowSums(rel_abund[,-1], na.rm = TRUE)
top_taxa <- names(sort(taxa_totals, decreasing = TRUE))[1:20]
rel_abund$Taxonomy <- ifelse(rel_abund$Taxonomy %in% top_taxa, rel_abund$Taxonomy, "Other")
rel_abund <- rel_abund %>%
  group_by(Taxonomy) %>%
  summarise(across(everything(), sum), .groups = "drop")

# --- 5. Bar plot ---
long_data <- rel_abund %>%
  pivot_longer(cols = -Taxonomy, names_to = "Sample", values_to = "Relative_Abundance")
p_bar <- ggplot(long_data, aes(x = Sample, y = Relative_Abundance, fill = Taxonomy)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative Abundance of Top 20 Taxa per Sample (Ordered by Diversity)",
       y = "Relative Abundance", x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 6. Heatmap (as grob for patchwork) ---
abundance_matrix <- rel_abund %>%
  column_to_rownames(var = "Taxonomy") %>%
  as.matrix()
heatmap_file <- tempfile(fileext = ".png")
png(heatmap_file, width = 1200, height = 900, res = 150)
pheatmap(log1p(abundance_matrix),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of Top 20 Taxa Relative Abundance",
         fontsize_row = 8,
         fontsize_col = 8)
dev.off()
p_heatmap <- ggplot() +
  annotation_custom(rasterGrob(readPNG(heatmap_file), interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()

# --- 7. Shannon diversity plot ---
p_shannon <- ggplot(shannon_df, aes(x = reorder(Sample, -Shannon_Diversity), y = Shannon_Diversity)) +
  geom_col(fill = "skyblue") +
  labs(title = "Shannon Diversity per Sample",
       y = "Shannon Index", x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 8. PCoA plot ---
bray_dist <- vegdist(t(count_matrix), method = "bray")
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_df <- data.frame(Sample = rownames(pcoa_res$points),
                      PC1 = pcoa_res$points[,1],
                      PC2 = pcoa_res$points[,2])
p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3, color = "tomato") +
  geom_text(vjust = -1, size = 3) +
  labs(title = "PCoA (Bray-Curtis)",
       x = "PCoA 1", y = "PCoA 2") +
  theme_minimal()

# --- 9. Rank abundance curve ---
taxa_means <- rowMeans(count_matrix)
taxa_sorted <- sort(taxa_means, decreasing = TRUE)
rank_df <- data.frame(Rank = 1:length(taxa_sorted), Abundance = taxa_sorted)
p_rank <- ggplot(rank_df, aes(x = Rank, y = Abundance)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Rank Abundance Curve",
       x = "Taxa Rank (most to least abundant)",
       y = "Mean Abundance (log scale)")

# --- 10. Combine all plots into 2x3 grid ---
combined_plot <- (p_bar | p_heatmap | p_shannon) /
                 (p_pcoa | p_rank | plot_spacer())
combined_plot

# Save combined plot
ggsave("combined_microbiome_plots.png", combined_plot, width = 18, height = 12, dpi = 300)
