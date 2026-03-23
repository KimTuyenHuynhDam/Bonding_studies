# ==============================================================================
# FINAL MANUSCRIPT FIGURES: MOLECULAR FIDELITY & STOCHASTICITY
# ==============================================================================

library(openxlsx)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# 1. SETUP & DATA LOADING
# ------------------------------------------------------------------------------
output_folder <- "Integrative_omics"
if(!dir.exists(output_folder)) dir.create(output_folder)

# Load Omics Data
rna_raw <- read.xlsx('20221221_fpkm_indv_IS lungs.xlsx')
rna_log <- log2(rna_raw[, 2:13] + 1)
rownames(rna_log) <- toupper(make.unique(as.character(rna_raw$X1)))

pro_raw <- read.xlsx('proteomics.xlsx')
pro_log <- log2(pro_raw[, 3:10] + 1)
pro_sym <- gsub(".*GN=([^ |]*).*", "\\1", pro_raw$Fasta.headers)
pro_sym[pro_sym == pro_raw$Fasta.headers] <- gsub(";.*", "", pro_raw$Protein.IDs[pro_sym == pro_raw$Fasta.headers])
rownames(pro_log) <- toupper(make.unique(as.character(pro_sym)))

common <- intersect(rownames(rna_log), rownames(pro_log))

# Calculate Distances
df_fidelity <- data.frame(
  Gene = common,
  RNA_B = rowMeans(rna_log[common, 1:6]),
  RNA_BD = rowMeans(rna_log[common, 7:12]),
  Prot_B = rowMeans(pro_log[common, 1:3]),
  Prot_BD = rowMeans(pro_log[common, 4:8])
) %>%
  mutate(
    Gap_Bonded = sqrt(RNA_B^2 + Prot_B^2),
    Gap_Disrupted = sqrt(RNA_BD^2 + Prot_BD^2)
  )

# 1.5 DYNAMIC FILTERING (High-Fidelity Threshold)
# ------------------------------------------------------------------------------
# We remove the bottom 10% of features to focus on robust signals
rna_cutoff  <- quantile(df_fidelity$RNA_B, 0.10)
prot_cutoff <- quantile(df_fidelity$Prot_B, 0.10)

df_final <- df_fidelity %>%
  filter(RNA_B > rna_cutoff & RNA_BD > rna_cutoff & 
           Prot_B > prot_cutoff & Prot_BD > prot_cutoff)

# Print final gene count
cat("Initial common genes:", length(common), "\n")
cat("Final high-fidelity genes after 10% cutoff:", nrow(df_final), "\n")

# 2. SUBFIGURE GENERATION
# ------------------------------------------------------------------------------

# --- SUBFIGURE A: Density Plot ---
ks_test <- ks.test(df_final$Gap_Bonded, df_final$Gap_Disrupted)

fig_density <- ggplot(df_final) +
  geom_density(aes(x = Gap_Bonded, fill = "Bonded"), alpha = 0.5) +
  geom_density(aes(x = Gap_Disrupted, fill = "Bond Disrupted"), alpha = 0.5) +
  scale_fill_manual(values = c("Bonded" = "#E41A1C", "Bond Disrupted" = "#377EB8"), name = "Group") +
  theme_bw() +
  labs(title = paste0("Stabilization of Molecular Distance by Bonding\n(KS-test p = ", round(ks_test$p.value, 3), ")"),
       x = "Regulatory Gap (Distance)", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")

print(fig_density)
ggsave(file.path(output_folder, "Fidelity_Density_Shift.jpg"), fig_density, width = 7, height = 5, dpi = 300)


# --- SUBFIGURE B: Boxplot ---
df_long <- df_final %>%
  dplyr::select(Gene, Gap_Bonded, Gap_Disrupted) %>%
  pivot_longer(cols = starts_with("Gap"), names_to = "Group", values_to = "Gap") %>%
  mutate(Group = factor(ifelse(Group == "Gap_Bonded", "Bonded", "Bond Disrupted"), 
                        levels = c("Bonded", "Bond Disrupted")))

t_test_res <- t.test(df_final$Gap_Disrupted, df_final$Gap_Bonded, paired = TRUE)

fig_boxplot <- ggplot(df_long, aes(x = Group, y = Gap, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  stat_compare_means(method = "t.test", paired = TRUE, label = "p.format", label.x = 1.4) +
  scale_fill_manual(values = c("Bonded" = "#E41A1C", "Bond Disrupted" = "#377EB8")) +
  theme_bw() +
  labs(title = "Systemic Regulatory Drift Post-Disruption",
       x = "", y = "Molecular Distance (mRNA-Protein)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

print(fig_boxplot)
ggsave(file.path(output_folder, "Fidelity_Boxplot_Comparison.jpg"), fig_boxplot, width = 6, height = 7, dpi = 300)


# --- SUBFIGURE C: Bonded Scatter ---
cor_B <- cor.test(df_final$RNA_B, df_final$Prot_B)

fig_scatter_B <- ggplot(df_final, aes(x = RNA_B, y = Prot_B)) +
  geom_point(alpha = 0.5, color = "#E41A1C") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  theme_bw() +
  labs(title = paste0("Bonded Fidelity (Baseline Support)\n(r = ", round(cor_B$estimate, 3), ", p = ", format.pval(cor_B$p.value), ")"),
       x = "mRNA Abundance (Log2)", y = "Protein Abundance (Log2)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))

print(fig_scatter_B)
ggsave(file.path(output_folder, "Fidelity_Scatter_Bonded.jpg"), fig_scatter_B, width = 6, height = 6, dpi = 300)


# --- SUBFIGURE D: Disrupted Scatter ---
cor_BD <- cor.test(df_final$RNA_BD, df_final$Prot_BD)

fig_scatter_BD <- ggplot(df_final, aes(x = RNA_BD, y = Prot_BD)) +
  geom_point(alpha = 0.5, color = "#377EB8") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  theme_bw() +
  labs(title = paste0("Bond Disrupted Fidelity (Regulatory Drift)\n(r = ", round(cor_BD$estimate, 3), ", p = ", format.pval(cor_BD$p.value), ")"),
       x = "mRNA Abundance (Log2)", y = "Protein Abundance (Log2)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))

print(fig_scatter_BD)
ggsave(file.path(output_folder, "Fidelity_Scatter_Disrupted.jpg"), fig_scatter_BD, width = 6, height = 6, dpi = 300)