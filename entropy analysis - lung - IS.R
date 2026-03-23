# ==============================================================================
# Refined RNAseq Analysis: Lung Tissue (B, BD, V)
# Project: IS-Bonding Studies
# Final Version: COMPLETE SUITE (Diversity, TII, Volcano,  Heatmaps, GO)
# ==============================================================================

# 1. SETUP & LIBRARIES
set.seed(42) # ENSURES P-VALUES AND RAREFACTION ARE IDENTICAL EVERY RUN
output_dir <- "IS-BD-lung-RNAseq"
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
library(org.Mm.eg.db) # Mouse database for lung tissue

# 2. DATA LOADING & CLEANING
lung_data <- read.xlsx('20221221_fpkm_indv_IS lungs.xlsx')
numeric_data <- lung_data[, 2:18] # B1-B6, BD1-BD6, V2-V6
log_data <- log2(numeric_data + 1)

# Handle Duplicate Gene Symbols for Rownames
rownames(log_data) <- make.unique(as.character(lung_data$X1))
rownames(numeric_data) <- rownames(log_data)

# 3. GROUPING & REFINED LABELS
sample_names <- colnames(numeric_data)
groups <- case_when(
  grepl("^B[0-9]",  sample_names) ~ "Bonded",
  grepl("^BD",      sample_names) ~ "Bond Disrupted",
  grepl("^V",       sample_names) ~ "Virgin"
)
groups <- factor(groups, levels = c("Bonded", "Bond Disrupted", "Virgin"))

# Clean Names: B_1-6, BD_1-6, V_1-5
clean_labels <- c(paste0("B_", 1:6), paste0("BD_", 1:6), paste0("V_", 1:5))
colnames(log_data) <- clean_labels
colnames(numeric_data) <- clean_labels

# Statistics comparisons
my_comparisons <- list(
  c("Bonded", "Bond Disrupted"), 
  c("Bonded", "Virgin"), 
  c("Bond Disrupted", "Virgin")
)

# ==============================================================================
# SECTION 1: COMPLETE DIVERSITY SUITE
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

write.csv(diversity_results, file.path(output_dir, "01_Lung_Full_Diversity_Data.csv"), row.names = FALSE)

plot_diversity <- function(df, y_col, title_text) {
  ggplot(df, aes_string(x = "Group", y = y_col, fill = "Group")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") +
    theme_minimal() + scale_fill_brewer(palette = "Set1") +
    labs(title = title_text, y = y_col) + theme(legend.position = "none")
}

ggsave(file.path(output_dir, "01a_Shannon_Normal.jpeg"),   plot_diversity(diversity_results, "Shannon_Standard", "Shannon Entropy (Standard)"), width = 7, height = 5)
ggsave(file.path(output_dir, "01b_Shannon_Rarefied.jpeg"), plot_diversity(diversity_results, "Shannon_Rarefied", "Shannon Entropy (Rarefied)"), width = 7, height = 5)
ggsave(file.path(output_dir, "01c_Evenness.jpeg"),         plot_diversity(diversity_results, "Evenness", "Pielou's Evenness"), width = 7, height = 5)
ggsave(file.path(output_dir, "01d_Simpson.jpeg"),          plot_diversity(diversity_results, "Simpson", "Simpson Diversity"), width = 7, height = 5)
ggsave(file.path(output_dir, "01e_Richness.jpeg"),         plot_diversity(diversity_results, "Richness", "Gene Richness"), width = 7, height = 5)

# ==============================================================================
# SECTION 2: INSTABILITY & VOLCANO PLOT
# ==============================================================================

# A. Instability Index (TII)
baseline_mean <- rowMeans(log_data[, groups == "Bonded"])
diversity_results$Instability <- apply(log_data, 2, function(x) mean((x - baseline_mean)^2))

p_tii <- ggplot(diversity_results, aes(x = Group, y = Instability, fill = Group)) +
  geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, size = 3) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") +
  theme_minimal() + scale_fill_brewer(palette = "Set1") +
  labs(title = "Transcriptome Instability Index (TII)", y = "Instability Score (MSE)")
ggsave(file.path(output_dir, "02a_Lung_Instability.jpeg"), p_tii, width = 8, height = 6)

# B. Volcano Plot (Bond Disrupted vs Bonded)
results_volc <- data.frame(
  Gene = rownames(log_data),
  log2FC = rowMeans(log_data[, groups == "Bond Disrupted"]) - baseline_mean,
  p_val = apply(log_data, 1, function(x) t.test(x[groups=="Bond Disrupted"], x[groups=="Bonded"])$p.value)
) %>% mutate(Significant = ifelse(p_val < 0.05 & abs(log2FC) > 0.5, "Yes", "No"))

p_volc <- ggplot(results_volc, aes(x = log2FC, y = -log10(p_val), color = Significant)) +
  geom_point(alpha = 0.4) + theme_minimal() + scale_color_manual(values = c("grey", "red3")) +
  geom_text_repel(data = head(arrange(results_volc, p_val), 10), aes(label = Gene), color = "black") +
  labs(title = "Volcano Plot: Lung (BD vs B)", x = "Log2 Fold Change", y = "-Log10 P-value")
ggsave(file.path(output_dir, "02b_Lung_Volcano.jpeg"), p_volc, width = 8, height = 6)

# ==============================================================================
# SECTION 3: SUPERVISED HEATMAPS & CORRELATION
# ==============================================================================

# A. Sample Correlation
jpeg(file.path(output_dir, "03a_Lung_Correlation.jpeg"), width = 1000, height = 1000, res = 150)
pheatmap(cor(log_data), display_numbers = TRUE, main = "Lung Sample Correlation")
dev.off()

# B. Supervised Expression Heatmap
gene_scores <- rowMeans((log_data[, groups == "Bond Disrupted"] - baseline_mean)^2)
marker_list <- data.frame(Gene = rownames(log_data), Instability_Score = gene_scores) %>% arrange(desc(Instability_Score))
write.csv(marker_list, file.path(output_dir, "03b_Lung_Instability_Markers.csv"), row.names = FALSE)

top_50 <- head(marker_list$Gene, 50)
annotation_df <- data.frame(Group = groups, row.names = clean_labels)
jpeg(file.path(output_dir, "03c_Lung_Heatmap_Supervised.jpeg"), width = 1000, height = 1200, res = 150)
pheatmap(log_data[top_50, ], annotation_col = annotation_df, scale = "row", 
         cluster_cols = FALSE, show_colnames = TRUE, main = "Top 50 Unstable Lung Genes",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()



# ==============================================================================
# SECTION 4: TOP MARKER CHECK (With Significance Brackets)
# ==============================================================================

# 1. Identify the absolute most unstable gene
top_gene <- marker_list$Gene[1]

# 2. Prepare the data for plotting
# We use diversity_results because it already contains Sample, Group, and clean labels
top_gene_data <- diversity_results %>%
  mutate(Expression = as.numeric(log_data[top_gene, ]))

# 3. Create the Boxplot with p-values for each pair
p_top_marker <- ggplot(top_gene_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) +
  # Add the statistical comparisons
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test", 
                     label = "p.format",
                     tip_length = 0.02) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = paste("Top Unstable Marker in Lung:", top_gene),
       subtitle = "P-values calculated via Wilcoxon test for each group pair",
       y = "log2(FPKM + 1)",
       x = "Experimental Group") +
  theme(legend.position = "none")

# 4. Save the high-resolution version
ggsave(file.path(output_dir, "04c_Lung_Top_Marker_Check_with_P.jpeg"), 
       p_top_marker, width = 8, height = 7, dpi = 300)

# Display result
print(p_top_marker)



# ==============================================================================
# SECTION 5: FUNCTIONAL PATHWAYS (GO)
# ==============================================================================

go_res <- enrichGO(gene = head(marker_list$Gene, 200), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
ggsave(file.path(output_dir, "05_Lung_Pathways.jpeg"), dotplot(go_res), width = 10, height = 7)

