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
# SECTION 3: SUPERVISED HEATMAP & PCA
# ==============================================================================

top_50 <- results_df %>% arrange(desc(Instability_Score)) %>% head(50) %>% pull(Gene)
heatmap_matrix <- log_data[top_50, ]
heatmap_matrix <- heatmap_matrix[apply(heatmap_matrix, 1, function(x) var(x, na.rm=T) > 0), ]

annotation_df <- data.frame(Group = groups, row.names = clean_labels)
jpeg(file.path(output_dir, "03_Proteome_Heatmap_Supervised.jpeg"), width = 1000, height = 1200, res = 150)
pheatmap(heatmap_matrix, annotation_col = annotation_df, scale = "row", 
         cluster_cols = FALSE, show_colnames = TRUE, fontsize_col = 12,
         main = "Top 50 Unstable Proteins (Grouped)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

pca_res <- prcomp(t(log_data[apply(log_data, 1, var) > 0, ]), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Group = groups)
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 5) + theme_minimal() + labs(title = "PCA: Proteome Separation")
ggsave(file.path(output_dir, "03b_Proteome_PCA.jpeg"), p_pca, width = 8, height = 6)

# ==============================================================================
# SECTION 4: TOP MARKER CHECK
# ==============================================================================

top_marker_name <- results_df$Gene[which.max(results_df$Instability_Score)]
top_marker_data <- data.frame(Group = groups, Expression = as.numeric(log_data[top_marker_name, ]))

p_top_marker <- ggplot(top_marker_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  theme_pubr() + scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  labs(title = paste("Top Unstable Protein:", top_marker_name),
       subtitle = "This protein most strongly defines the Bond Disruption shift",
       y = "log2 Intensity") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "04_Proteome_Top_Marker_with_P.jpeg"), p_top_marker, width = 7, height = 6, dpi = 300)

# ==============================================================================
# SECTION 5: GO PATHWAYS
# ==============================================================================

target_proteins <- results_df %>% dplyr::filter(p_val < 0.1) %>% pull(Gene)
go_res <- enrichGO(gene = target_proteins, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
ggsave(file.path(output_dir, "05_Proteome_GO_Pathways.jpeg"), dotplot(go_res), width = 10, height = 8)

# ==============================================================================
# FIXED Z-SCORE DRIFT PLOT
# ==============================================================================

p_values_f <- z_plot_data %>%
  group_by(Gene) %>%
  summarise(p_val = var.test(Z_Score ~ Group)$p.value) %>%
  mutate(p_label = sprintf("p (F-test) = %.3f", p_val))

z_plot_data_with_p <- z_plot_data %>%
  left_join(p_values_f, by = "Gene") %>%
  mutate(Gene_P = paste0(Gene, "\n", p_label))

p_zscore_final <- ggplot(z_plot_data_with_p, aes(x = Group, y = Z_Score, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, shape = 21, color = "black") +
  facet_wrap(~Gene_P, scales = "free_y") + 
  theme_pubr() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  labs(title = "Proteomic Drift (Z-Scores)",
       subtitle = "P-value (F-test) measures if the BD group is significantly more 'scattered' than Bonded",
       y = "Z-Score (Distance from Bonded Mean)") +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 10))

ggsave(file.path(output_dir, "04e_Zscore_Drift_Final_P.jpeg"), 
       p_zscore_final, width = 12, height = 6, dpi = 300)

# ==============================================================================
# GLOBAL DRIFT (TOP 50)
# ==============================================================================

top_50_genes_drift <- head(results_df %>% arrange(desc(Instability_Score)), 50)$Gene
z_matrix_drift <- apply(log_data[top_50_genes_drift, ], 1, function(x) {
  b_mu <- mean(x[1:3]); b_sd <- sd(x[1:3])
  if(b_sd == 0) b_sd <- 0.01
  return(abs((x - b_mu) / b_sd))
})
drift_df <- data.frame(Group = groups, Drift = rowMeans(z_matrix_drift))

p_global_drift <- ggplot(drift_df, aes(x = Group, y = Drift, fill = Group)) +
  geom_boxplot() + geom_jitter(width = 0.1, size = 3) +
  stat_compare_means(method = "t.test") + 
  labs(title = "Global Proteomic Drift (Top 50 Genes)",
       y = "Mean Absolute Z-Score (Distance from Baseline)")

ggsave(file.path(output_dir, "04f_Zscore_Drift_top50_genes_Final_P.jpeg"), 
       p_global_drift, width = 12, height = 6, dpi = 300)

# ==============================================================================
# AUTOMATED SIGNIFICANCE PLOTTER (ALL HERO MARKERS)
# ==============================================================================

significant_markers <- results_df %>%
  dplyr::filter(Significant == "Significant") %>%
  arrange(desc(Instability_Score))

n_sig <- nrow(significant_markers)
top_n <- min(n_sig, 12) 
selected_genes <- head(significant_markers$Gene, top_n)

plot_list_sig <- list()
for (gene in selected_genes) {
  b_mu <- mean(as.numeric(log_data[gene, 1:3]))
  b_sd <- sd(as.numeric(log_data[gene, 1:3]))
  if(is.na(b_sd) || b_sd == 0) b_sd <- 0.01
  z_vals <- (as.numeric(log_data[gene, ]) - b_mu) / b_sd
  plot_list_sig[[gene]] <- data.frame(Gene = gene, Group = groups, Z_Score = z_vals)
}
multi_sig_data <- do.call(rbind, plot_list_sig)

p_all_sig <- ggplot(multi_sig_data, aes(x = Group, y = Z_Score, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, shape = 21, color = "black", stroke = 0.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = 3) + 
  stat_compare_means(method = "t.test", label = "p.format", size = 3) +
  theme_pubr() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  labs(title = "Proteomic Hero Markers: Significant Instability",
       subtitle = "Genes filtered by p < 0.05 and Fold Change > 1.5. Plotted as Z-score drift.",
       y = "Z-Score (Distance from Bonded Mean)") +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "04j_All_Significant_Markers_Zscore.jpeg"), 
       p_all_sig, width = 10, height = 3.5 * ceiling(top_n/3), dpi = 300)

# ==============================================================================
# PROTEOMICS HERO TRIO
# ==============================================================================
# 1. Filter and immediately "clean" the data frame to avoid Rle/List errors
results_sig_trio <- results_df %>%
  dplyr::filter(p_val < 0.05 & abs(log2FC) > fc_threshold) %>%
  as.data.frame()

# Ensure Gene column is a simple character vector, not a factor or Rle
results_sig_trio$Gene <- as.character(results_sig_trio$Gene)

# 2. Winner 1: Most Unstable
hero_instability <- results_sig_trio %>% 
  dplyr::arrange(desc(Instability_Score)) %>% 
  head(1)
primary_gene <- hero_instability$Gene

# 3. Winner 2: Strongest Magnitude (Check for Duplicates)
# We sort the whole list and find the first one that isn't the primary_gene
mag_candidates <- results_sig_trio %>% 
  dplyr::arrange(desc(abs(log2FC)))

if (mag_candidates$Gene[1] == primary_gene) {
  hero_magnitude <- mag_candidates[2, , drop = FALSE] # Base R slice is safer here
} else {
  hero_magnitude <- mag_candidates[1, , drop = FALSE]
}
secondary_gene <- hero_magnitude$Gene

# 4. Winner 3: Most Confident (Check against both Winner 1 and 2)
conf_candidates <- results_sig_trio %>% 
  dplyr::arrange(p_val)

# Logic: Filter out the first two winners, then take the top of what's left
hero_confidence <- conf_candidates %>%
  dplyr::filter(!(Gene %in% c(primary_gene, secondary_gene))) %>%
  head(1)

# 5. Combine into the list for your plotting loop
heroes_list <- list(
  "Most Unstable" = hero_instability,
  "Strongest Magnitude" = hero_magnitude,
  "Most Confident" = hero_confidence
)


plot_list_trio <- list()
for (i in 1:length(heroes_list)) {
  name <- names(heroes_list)[i]
  gene_info <- heroes_list[[i]]
  if(nrow(gene_info) == 0) next
  gene_name <- as.character(gene_info$Gene)
  real_fc <- 2^(gene_info$log2FC)
  
  b_mu <- mean(as.numeric(log_data[gene_name, 1:3]))
  b_sd <- sd(as.numeric(log_data[gene_name, 1:3]))
  if(is.na(b_sd) || b_sd == 0) b_sd <- 0.01
  z_vals <- (as.numeric(log_data[gene_name, ]) - b_mu) / b_sd
  
  df_p <- data.frame(Group = factor(groups, levels = c("Bonded", "Bond Disrupted")), Z = z_vals)
  stat_text <- paste0("p-val = ", format.pval(gene_info$p_val, digits=3), 
                      "\nReal FC = ", round(real_fc, 2), "x",
                      "\nInstab. = ", round(gene_info$Instability_Score, 1))
  
  p <- ggplot(df_p, aes(x = Group, y = Z, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 3, shape = 21, color = "black", stroke = 1) +
    theme_pubr() +
    scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
    labs(title = paste0(gsub("_", " ", name), ": ", gene_name),
         subtitle = stat_text, y = "Z-Score (Distance from Baseline)") +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9, family = "mono"))
  plot_list_trio[[name]] <- p
}

p_final_trio_realFC <- wrap_plots(plot_list_trio, ncol = 3) + 
  plot_annotation(title = "Top Significant Proteomic Markers (Real Fold Change)")

ggsave(file.path(output_dir, "04k_Proteome_Hero_Trio_RealFC.jpeg"), 
       p_final_trio_realFC, width = 14, height = 5.5, dpi = 300)