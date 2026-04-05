# ==============================================================================
# Consolidated RNAseq Analysis: Serum-treated A549 Cells

# ==============================================================================

# 1. SETUP & LIBRARIES
set.seed(100) # ENSURES P-VALUES AND RAREFACTION ARE IDENTICAL EVERY RUN
output_dir <- "IS-BD-sera-RNAseq"
if(!dir.exists(output_dir)) dir.create(output_dir)

library(openxlsx)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(vegan)
library(ggsignif)
library(car)

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

write.csv(diversity_results, file.path(output_dir, "Sera_Diversity_Data.csv"), row.names = FALSE)

plot_diversity <- function(df, y_col, title_text) {
  ggplot(df, aes_string(x = "Group", y = y_col, fill = "Group")) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
    theme_minimal() + scale_fill_brewer(palette = "Set2") +
    labs(title = title_text, y = y_col) + theme(legend.position = "none")
}

ggsave(file.path(output_dir, "Shannon_Normal.jpeg"), plot_diversity(diversity_results, "Shannon_Standard", "Shannon Entropy (Standard)"), width = 7, height = 5)
ggsave(file.path(output_dir, "Shannon_Rarefied.jpeg"), plot_diversity(diversity_results, "Shannon_Rarefied", "Shannon Entropy (Rarefied)"), width = 7, height = 5)
ggsave(file.path(output_dir, "Evenness.jpeg"), plot_diversity(diversity_results, "Evenness", "Pielou's Evenness"), width = 7, height = 5)
ggsave(file.path(output_dir, "Simpson.jpeg"), plot_diversity(diversity_results, "Simpson", "Simpson Diversity"), width = 7, height = 5)
ggsave(file.path(output_dir, "Richness.jpeg"), plot_diversity(diversity_results, "Richness", "Gene Richness"), width = 7, height = 5)

# ==============================================================================
# SECTION 2: VARIANCE AND FLUCTUATION ANALYSIS (STABILITY VS NOISE)
# ==============================================================================


# 1. Descriptive Statistics: Calculate Standard Deviation (SD) and CV
# This gives you the unbiased numbers to report in your tables
variance_summary <- diversity_results %>%
  group_by(Group) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Standard),
    SD_Shannon = sd(Shannon_Standard),
    Variance = var(Shannon_Standard),
    CV_Percent = (SD_Shannon / Mean_Shannon) * 100
  ) %>%
  arrange(CV_Percent) # Orders from neatest (lowest CV) to most fluctuating

print("--- Descriptive Fluctuation (Coefficient of Variation) ---")
print(variance_summary)

# 2. Inferential Statistics: Pairwise F-Tests for Equality of Variance
# We compare the variance of the "neat" Bonded group against the others.

bonded_data <- diversity_results$Shannon_Standard[diversity_results$Group == "Bonded"]
bd_data <- diversity_results$Shannon_Standard[diversity_results$Group == "Bond Disrupted"]
virgin_data <- diversity_results$Shannon_Standard[diversity_results$Group == "Virgin"]

print("--- F-Test: Bonded vs Bond Disrupted Variance ---")
ftest_bd <- var.test(bonded_data, bd_data)

print(ftest_bd)

print("--- F-Test: Bonded vs Virgin Variance ---")
ftest_v <- var.test(bonded_data, virgin_data)
print(ftest_v)

# 3. Global Robust Variance Test (Levene's Test)
# This tests if the fluctuation is significantly different across ALL groups
print("--- Levene's Test for Homogeneity of Variance (Global) ---")
levene_res <- leveneTest(Shannon_Standard ~ Group, data = diversity_results)
print(levene_res)

# ==============================================================================
# SECTION 3: CV BAR CHART WITH EXACT 4-DIGIT P-VALUES
# ==============================================================================


# 1. Define the group order to keep the graph logical
group_order <- c("Bonded", "Bond Disrupted", "Virgin", "Control (FBS)")

# 2. Re-summarize the CV data to ensure correct factor ordering
variance_summary <- diversity_results %>%
  group_by(Group) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Standard),
    SD_Shannon = sd(Shannon_Standard),
    CV_Percent = (SD_Shannon / Mean_Shannon) * 100
  ) %>%
  mutate(Group = factor(Group, levels = group_order))

write.xlsx(variance_summary, file.path(output_dir, "Sera_Variance_Summary.xlsx"))

# 3. Helper function: Calculate F-test and force 4-decimal formatting
get_var_pval <- function(g1, g2) {
  v1 <- diversity_results$Shannon_Standard[diversity_results$Group == g1]
  v2 <- diversity_results$Shannon_Standard[diversity_results$Group == g2]
  res <- var.test(v1, v2)
  
  # Format the p-value to exactly 4 decimal places
  return(sprintf("p = %.4f", res$p.value))
}

# 4. Calculate exact p-values for all 6 combinations
p_B_BD <- get_var_pval("Bonded", "Bond Disrupted")
p_B_V  <- get_var_pval("Bonded", "Virgin")
p_B_C  <- get_var_pval("Bonded", "Control (FBS)")
p_BD_V <- get_var_pval("Bond Disrupted", "Virgin")
p_BD_C <- get_var_pval("Bond Disrupted", "Control (FBS)")
p_V_C  <- get_var_pval("Virgin", "Control (FBS)")

# 5. Build the final plot
plot_cv_stats <- ggplot(variance_summary, aes(x = Group, y = CV_Percent, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.9, color = "black", width = 0.6) +
  geom_text(aes(label = round(CV_Percent, 3)), vjust = -0.5, size = 4, fontface = "bold") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  
  # Expand the Y-axis to make room for 6 layers of brackets
  coord_cartesian(ylim = c(0, 0.85)) + 
  
  labs(title = "Transcriptomic Fluctuation (Variance) by Group",
       subtitle = "Pairwise F-tests comparing variance of Shannon Entropy",
       y = "Coefficient of Variation (%)", 
       x = "Group") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 11, face = "bold")) +
  
  # 6. Add the precisely formatted brackets
  geom_signif(
    comparisons = list(
      c("Bonded", "Bond Disrupted"),     
      c("Virgin", "Control (FBS)"),      
      c("Bond Disrupted", "Virgin"),     
      c("Bonded", "Virgin"),             
      c("Bond Disrupted", "Control (FBS)"), 
      c("Bonded", "Control (FBS)")       
    ),
    y_position = c(0.50, 0.50, 0.57, 0.64, 0.71, 0.78), 
    annotations = c(p_B_BD, p_V_C, p_BD_V, p_B_V, p_BD_C, p_B_C),
    tip_length = 0.02, 
    textsize = 3.5,
    vjust = -0.2
  )

# 7. Save the graph
ggsave(file.path(output_dir, "Fluctuation_CV_Exact_Pvalues.jpeg"), plot_cv_stats, width = 8, height = 6, dpi = 300)


