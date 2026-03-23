# ==============================================================================
# Final Proteomics Analysis: Bonded (B) vs. Bond Disrupted (BD)
# Project: IS-Bonding Studies
# Features: Supervised Heatmap, Volcano, GO Pathways, and Top Marker with P-Value
# ==============================================================================

# 1. SETUP & LIBRARIES
set.seed(42) 
output_dir <- "IS-BD-proteomics"
if(!dir.exists(output_dir)) dir.create(output_dir)

library(openxlsx)
library(dplyr)
library(ggplot2)
library(vegan)
library(ggpubr)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)

# 2. DATA LOADING & CLEANING
pro_data <- read.xlsx('proteomics.xlsx')
numeric_data <- pro_data[, 3:10] 
log_data <- log2(numeric_data + 1)

gene_symbols <- gsub(".*GN=([^ |]*).*", "\\1", pro_data$Fasta.headers)
gene_symbols[gene_symbols == pro_data$Fasta.headers] <- gsub(";.*", "", pro_data$Protein.IDs[gene_symbols == pro_data$Fasta.headers])
unique_genes <- make.unique(as.character(gene_symbols))
rownames(log_data) <- unique_genes
rownames(numeric_data) <- unique_genes

clean_labels <- c("B_1", "B_2", "B_3", "BD_1", "BD_2", "BD_3", "BD_4", "BD_5")
colnames(log_data) <- clean_labels
colnames(numeric_data) <- clean_labels

groups <- factor(c(rep("Bonded", 3), rep("Bond Disrupted", 5)), levels = c("Bonded", "Bond Disrupted"))
my_comparisons <- list(c("Bonded", "Bond Disrupted"))

# ==============================================================================
# SECTION 1: DIVERSITY SUITE
# ==============================================================================

t_data <- t(numeric_data)
counts_data <- round(t_data)
min_depth <- min(rowSums(counts_data))

diversity_results <- data.frame(
  Sample = clean_labels,
  Group = groups,
  Shannon_Normal = diversity(t_data, index = "shannon"),
  Shannon_Rarefied = diversity(rrarefy(counts_data, min_depth), index = "shannon"),
  Evenness = diversity(t_data, index = "shannon") / log(specnumber(t_data))
)

plot_pro_metric <- function(df, col, title) {
  ggplot(df, aes_string(x = "Group", y = col, fill = "Group")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 3) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
    theme_pubr() + scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
    labs(title = title, y = col) + theme(legend.position = "none")
}

ggsave(file.path(output_dir, "01a_Shannon_Standard.jpeg"), plot_pro_metric(diversity_results, "Shannon_Normal", "Shannon Entropy"), width = 6, height = 5)
ggsave(file.path(output_dir, "01b_Shannon_Rarefied.jpeg"), plot_pro_metric(diversity_results, "Shannon_Rarefied", "Rarefied Shannon Entropy"), width = 6, height = 5)
ggsave(file.path(output_dir, "01c_Evenness.jpeg"), plot_pro_metric(diversity_results, "Evenness", "Proteome Evenness"), width = 6, height = 5)

# ==============================================================================
# SECTION 2: INSTABILITY (TII) & VOLCANO
# ==============================================================================

baseline_mean <- rowMeans(log_data[, 1:3])
results_df <- data.frame(
  Gene = unique_genes,
  log2FC = rowMeans(log_data[, 4:8]) - baseline_mean,
  p_val = apply(log_data, 1, function(x) {
    if(sd(x[4:8]) == 0 && sd(x[1:3]) == 0) return(1)
    t.test(x[4:8], x[1:3])$p.value
  }),
  Instability_Score = apply(log_data[, 4:8], 1, function(x) mean((x - baseline_mean)^2))
)

fc_threshold <- log2(1.5) 
results_df <- results_df %>%
  mutate(Significant = ifelse(p_val < 0.05 & abs(log2FC) > fc_threshold, "Significant", "Non-Significant"))

p_volc <- ggplot(results_df, aes(x = log2FC, y = -log10(p_val), color = Significant)) +
  geom_point(aes(size = Instability_Score), alpha = 0.5) +
  scale_color_manual(values = c("grey", "red3")) + theme_minimal() +
  geom_text_repel(data = filter(results_df, Significant == "Significant"), aes(label = Gene), fontface = "bold") +
  labs(title = "Proteome Volcano Plot", subtitle = "Significant markers labeled")
ggsave(file.path(output_dir, "02_Proteome_Volcano.jpeg"), p_volc, width = 8, height = 7)

write.csv(results_df %>% arrange(desc(Instability_Score)), file.path(output_dir, "02_Proteome_Master_Results.csv"), row.names = FALSE)

# ==============================================================================
# DATA BRIDGE: Create z_plot_data for subsequent sections
# ==============================================================================
top_unstable_genes <- head(results_df %>% arrange(desc(Instability_Score)), 12)$Gene
z_plot_list <- list()
for (gene in top_unstable_genes) {
  b_mu <- mean(as.numeric(log_data[gene, 1:3]))
  b_sd <- sd(as.numeric(log_data[gene, 1:3]))
  if(is.na(b_sd) || b_sd == 0) b_sd <- 0.01 
  z_vals <- (as.numeric(log_data[gene, ]) - b_mu) / b_sd
  z_plot_list[[gene]] <- data.frame(Gene = gene, Group = groups, Z_Score = z_vals)
}
z_plot_data <- do.call(rbind, z_plot_list)

# ==============================================================================
# SECTION 3: SUPERVISED HEATMAP 
# ==============================================================================

top_50 <- results_df %>% arrange(desc(Instability_Score)) %>% head(50) %>% pull(Gene)
heatmap_matrix <- log_data[top_50, ]
heatmap_matrix <- heatmap_matrix[apply(heatmap_matrix, 1, function(x) var(x, na.rm=T) > 0), ]

annotation_df <- data.frame(Group = groups, row.names = clean_labels)
jpeg(file.path(output_dir, "03_Proteome_Heatmap_Supervised.jpeg"), width = 1000, height = 1200, res = 150)
pheatmap(heatmap_matrix, annotation_col = annotation_df, scale = "row", 
         cluster_cols = FALSE, show_colnames = TRUE, fontsize_col = 12,
         main = "Top 50 Unstable Proteins",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()



