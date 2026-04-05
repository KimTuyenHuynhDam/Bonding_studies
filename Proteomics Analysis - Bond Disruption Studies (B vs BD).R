# ==============================================================================
# Final Proteomics Analysis: Bonded (B) vs. Bond Disrupted (BD)
# Project: IS-Bonding Studies
# Features: Supervised Heatmap, Volcano, GO Pathways, Diversity & Variance
# ==============================================================================

# 1. SETUP & LIBRARIES
set.seed(100) 
output_dir <- "IS-BD-proteomics"
if(!dir.exists(output_dir)) dir.create(output_dir)

library(openxlsx)
library(dplyr)
library(ggplot2)
library(vegan)
library(ggpubr)
library(ggrepel)
library(ggsignif) # Required for the variance statistical brackets


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
# SECTION 2: FLUCTUATION & VARIANCE ANALYSIS (STABILITY VS NOISE)
# ==============================================================================

# 1. Summarize the CV data for the proteomics groups
variance_summary <- diversity_results %>%
  group_by(Group) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Normal),
    SD_Shannon = sd(Shannon_Normal),
    Variance = var(Shannon_Normal),
    CV_Percent = (SD_Shannon / Mean_Shannon) * 100
  ) %>%
  mutate(Group = factor(Group, levels = c("Bonded", "Bond Disrupted")))

# Export CV results to Excel
write.xlsx(variance_summary, file.path(output_dir, "02_Proteomics_Variance_Summary.xlsx"))
print("Saved proteomics variance summary to Excel.")

# 2. Helper function: Calculate F-test (variance) and force exactly 4-decimal formatting
get_var_pval <- function(g1, g2) {
  v1 <- diversity_results$Shannon_Normal[diversity_results$Group == g1]
  v2 <- diversity_results$Shannon_Normal[diversity_results$Group == g2]
  res <- var.test(v1, v2)
  return(sprintf("p = %.4f", res$p.value))
}

# 3. Calculate exact p-value for the single comparison (Bonded vs Bond Disrupted)
p_B_BD <- get_var_pval("Bonded", "Bond Disrupted")

# 4. Dynamically calculate y-axis height for the bracket
max_cv <- max(variance_summary$CV_Percent)
y_pos1 <- max_cv * 1.20
y_max_limit <- max_cv * 1.40

# 5. Build the final CV plot
plot_cv_stats <- ggplot(variance_summary, aes(x = Group, y = CV_Percent, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.9, color = "black", width = 0.5) +
  geom_text(aes(label = round(CV_Percent, 3)), vjust = -0.5, size = 4, fontface = "bold") +
  theme_pubr() + 
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) + # Matches the exact hex colors from your boxplots
  
  # Expand the Y-axis to make room for the statistical bracket
  coord_cartesian(ylim = c(0, y_max_limit)) + 
  
  labs(title = "Proteome Fluctuation (Variance) by Group",
       subtitle = "F-test comparing variance of Shannon Entropy",
       y = "Coefficient of Variation (%)", 
       x = "Group") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 11, face = "bold")) +
  
  # 6. Add the precisely formatted bracket (Only 1 bracket needed for 2 groups)
  geom_signif(
    comparisons = list(c("Bonded", "Bond Disrupted")),
    y_position = y_pos1, 
    annotations = p_B_BD,
    tip_length = 0.02, 
    textsize = 4,
    vjust = -0.2
  )

# 7. Save the CV Bar Chart
ggsave(file.path(output_dir, "02a_Fluctuation_CV_Exact_Pvalues.jpeg"), plot_cv_stats, width = 6, height = 5, dpi = 300)

print("Statistically annotated Proteomics CV graph successfully saved!")