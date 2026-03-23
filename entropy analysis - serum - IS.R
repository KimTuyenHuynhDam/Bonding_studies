# ==============================================================================
# Consolidated RNAseq Analysis: Serum-treated A549 Cells
# Final Version: Supervised Grouping & Top Marker with P-Values
# ==============================================================================

# 1. SETUP & LIBRARIES
set.seed(42) # ENSURES P-VALUES AND RAREFACTION ARE IDENTICAL EVERY RUN
output_dir <- "IS-BD-sera-RNAseq"
if(!dir.exists(output_dir)) dir.create(output_dir)

library(openxlsx)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(vegan)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

# 2. DATA LOADING & CLEANING
serum_data <- read.xlsx('20201116Asieh_fpkm_indv.xlsx')
numeric_data <- serum_data[, 2:25] # B_01 through C_06
log_data <- log2(numeric_data + 1)

# Ensure unique gene names
unique_gene_names <- make.unique(as.character(serum_data$X1))
rownames(log_data) <- unique_gene_names
rownames(numeric_data) <- unique_gene_names

# 3. GROUPING & ORDERING
sample_names <- colnames(numeric_data)
groups <- case_when(
  grepl("^B_",  sample_names) ~ "Bonded",
  grepl("^BD_", sample_names) ~ "Bond Disrupted",
  grepl("^V_",  sample_names) ~ "Virgin",
  grepl("^C_",  sample_names) ~ "Control (FBS)"
)
groups <- factor(groups, levels = c("Bonded", "Bond Disrupted", "Virgin", "Control (FBS)"))

# Create clean sample labels
clean_labels <- c(paste0("B_", 1:6), paste0("BD_", 1:6), paste0("V_", 1:6), paste0("C_", 1:6))
colnames(log_data) <- clean_labels
colnames(numeric_data) <- clean_labels

# Comparison pairs for statistics
my_comparisons <- list(
  c("Bonded", "Bond Disrupted"), 
  c("Bonded", "Virgin"), 
  c("Bond Disrupted", "Virgin"),
  c("Bond Disrupted", "Control (FBS)")
)

# ==============================================================================
# SECTION 1: GLOBAL DIVERSITY
# ==============================================================================

t_data <- t(numeric_data)
counts_data <- round(t_data)
min_depth <- min(rowSums(counts_data))

diversity_results <- data.frame(
  Sample = clean_labels,
  Group  = groups,
  Shannon_Standard = diversity(t_data, index = "shannon"),
  Shannon_Rarefied = diversity(rrarefy(counts_data, min_depth), index = "shannon"),
  Simpson          = diversity(t_data, index = "simpson"),
  Richness         = specnumber(t_data),
  Evenness         = diversity(t_data, index = "shannon") / log(specnumber(t_data))
)

write.csv(diversity_results, file.path(output_dir, "01_Sera_Diversity_Data.csv"), row.names = FALSE)

plot_diversity <- function(df, y_col, title_text) {
  ggplot(df, aes_string(x = "Group", y = y_col, fill = "Group")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
    theme_minimal() + scale_fill_brewer(palette = "Set2") +
    labs(title = title_text, y = y_col) + theme(legend.position = "none")
}

ggsave(file.path(output_dir, "01a_Shannon_Normal.jpeg"), plot_diversity(diversity_results, "Shannon_Standard", "Shannon Entropy (Standard)"), width = 7, height = 5)
ggsave(file.path(output_dir, "01b_Shannon_Rarefied.jpeg"), plot_diversity(diversity_results, "Shannon_Rarefied", "Shannon Entropy (Rarefied)"), width = 7, height = 5)
ggsave(file.path(output_dir, "01c_Evenness.jpeg"), plot_diversity(diversity_results, "Evenness", "Pielou's Evenness"), width = 7, height = 5)
ggsave(file.path(output_dir, "01d_Simpson.jpeg"), plot_diversity(diversity_results, "Simpson", "Simpson Diversity"), width = 7, height = 5)
ggsave(file.path(output_dir, "01e_Richness.jpeg"), plot_diversity(diversity_results, "Richness", "Gene Richness"), width = 7, height = 5)

# ==============================================================================
# SECTION 2: INSTABILITY INDEX (TII)
# ==============================================================================

baseline_mean <- rowMeans(log_data[, groups == "Bonded"])
diversity_results$Instability <- apply(log_data, 2, function(x) mean((x - baseline_mean)^2))

p_tii <- ggplot(diversity_results, aes(x = Group, y = Instability, fill = Group)) +
  geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  theme_minimal() + scale_fill_brewer(palette = "Set2") +
  labs(title = "Transcriptome Instability Index (TII)", y = "Instability Score (MSE)")
ggsave(file.path(output_dir, "02_Sera_Instability.jpeg"), p_tii, width = 8, height = 6)

gene_scores <- rowMeans((log_data[, groups == "Bond Disrupted"] - baseline_mean)^2)
marker_list <- data.frame(Gene = rownames(log_data), Instability_Score = gene_scores) %>% 
  arrange(desc(Instability_Score))
write.csv(marker_list, file.path(output_dir, "02_Sera_Instability_Markers.csv"), row.names = FALSE)

# ==============================================================================
# SECTION 3: SUPERVISED HEATMAPS & PCA
# ==============================================================================

pca_res <- prcomp(t(log_data[apply(log_data, 1, var) > 0, ]), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Group = groups)
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 5) + theme_bw() + labs(title = "PCA: Serum Sample Clustering")
ggsave(file.path(output_dir, "03_Sera_PCA.jpeg"), p_pca, width = 8, height = 6)

# Supervised Heatmap
top_50 <- head(marker_list$Gene, 50)
annotation_df <- data.frame(Group = groups, row.names = clean_labels)
jpeg(file.path(output_dir, "04_Sera_Heatmap_Supervised.jpeg"), width = 1000, height = 1200, res = 150)
pheatmap(log_data[top_50, ], annotation_col = annotation_df, scale = "row", 
         cluster_cols = FALSE, show_colnames = TRUE,
         main = "Top 50 Unstable Markers",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

# Correlation Heatmap
jpeg(file.path(output_dir, "05_Sera_Correlation.jpeg"), width = 1000, height = 1000, res = 150)
pheatmap(cor(log_data), display_numbers = TRUE, main = "Sample Correlation")
dev.off()

# ==============================================================================
# SECTION 4: FUNCTIONAL ANALYSIS
# ==============================================================================

go_res <- enrichGO(gene = head(marker_list$Gene, 200), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
ggsave(file.path(output_dir, "06_Sera_Pathways.jpeg"), dotplot(go_res), width = 10, height = 7)

# ==============================================================================
# SECTION 5: PCA LOADINGS & TOP MARKER CHECK (With P-Values)
# ==============================================================================

pc_loadings <- as.data.frame(pca_res$rotation) %>% mutate(Gene = rownames(.)) %>%
  dplyr::select(Gene, PC1, PC2) %>% arrange(desc(abs(PC2))) 
write.csv(pc_loadings, file.path(output_dir, "07_Sera_PCA_Loadings.csv"), row.names = FALSE)

# Extract and plot the single most influential gene with significance brackets
top_gene <- pc_loadings$Gene[1]
gene_plot_data <- data.frame(
  Group = groups, 
  Expression = as.numeric(log_data[top_gene, ])
)

p_top_marker <- ggplot(gene_plot_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2.5) +
  # Adding statistical significance brackets for each pair
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  theme_minimal() + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = paste("Top PCA Driver (PC2):", top_gene), 
       subtitle = "P-values calculated via t-test for each group pair",
       y = "log2(FPKM+1)") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "08_Sera_Top_Marker_with_P.jpeg"), p_top_marker, width = 8, height = 7, dpi = 300)

print("Full scientific analysis complete. All graphs (including Top Marker with P) are saved.")