# Bonding_studies: Social experiences modulate gene canalization

This repository contains the raw data, R analytical scripts, and high-resolution figures for the study: **"Social experiences modulate gene canalization in monogamous deer mice"**.

---

## 📊 Project Overview
The research investigates how pair bonding acts as a homeostatic anchor in monogamous deer mice (*Peromyscus californicus*). We demonstrate that bond disruption (widowhood) increases molecular noise (stochasticity) and promotes decanalization across the transcriptome and proteome.

---

## 📁 Repository Structure

### 🧪 Raw Data Files
These files provide the foundation for all calculations and figures:
* **`20201116Asieh_fpkm_indv.xlsx`**: RNA-seq FPKM matrix for human A549 cells cultured in deer mouse sera.
* **`20221221_fpkm_indv_IS lungs.xlsx`**: RNA-seq FPKM matrix for deer mouse lung tissue.
* **`proteomics.xlsx`**: Mass spectrometry intensity data for the blood serum proteome.

### 💻 Analysis Scripts (R)
Primary R scripts for calculating diversity and instability indices:
* **`entropy analysis - serum - IS.R`**: Processes A549 cell RNA-seq data to calculate Shannon’s entropy and TII.
* **`entropy analysis - lung - IS.R`**: Performs diversity and instability analysis on lung transcriptomes.
* **`Proteomics Analysis - Bond Disruption Studies.R`**: Analyzes serum proteome drift and instability.
* **`Integrative_omics.R`**: Multi-omic integration script to calculate the "regulatory gap" (Euclidean distance) between mRNA and protein levels.

### 🖼️ Output Directories
Processed results and publication-quality figures are organized into specialized folders:
* **`IS-BD-sera-RNAseq/`**: Results for cell culture transcriptomics.
* **`IS-BD-lung-RNAseq/`**: Results for deer mouse lung transcriptomics.
* **`IS-BD-proteomics/`**: Results for blood serum proteomics.
* **`Integrative_omics/`**: Integrated multi-omic figures and distribution analysis.

---

## 🛠️ Installation & Requirements
To replicate the analysis, clone the repository and ensure you have the following R environment:
* **R version**: 4.5.2
* **Required Libraries**: `vegan`, `ggplot2`, `dplyr`, `pheatmap`, `clusterProfiler`, `org.Mm.eg.db`, `org.Hs.eg.db`.

---  

## 🧬 Key Methodology
The study utilizes several unique indices to quantify biological decanalization:

### 1. Transcriptome Instability Index (TII)
The TII measures the cumulative deviation of an experimental sample from the "Bonded" control centroid. It is calculated as the **Mean Squared Error (MSE)**:

$$TII_i = \frac{1}{G} \sum_{g=1}^{G} (x_{i,g} - \mu_{g,B})^2$$

* **$G$**: Total number of genes/proteins.
* **$x_{i,g}$**: Log-transformed expression for sample $i$ and gene $g$.
* **$\mu_{g,B}$**: Baseline mean established from the Bonded group.

### 2. Rarefied Shannon Entropy
To account for varying sequencing depths or library sizes, counts were rarefied to the minimum sample depth using **`rrarefy`**. This ensures unbiased comparisons of stochasticity across experimental groups.

### 3. Regulatory Gap
We quantified the discordance between transcriptional and translational outputs—the **"regulatory gap"**—by calculating the Euclidean distance ($d$) for each gene in a two-dimensional multi-omic state space:

$$d = \sqrt{(\text{mRNA}_{log2})^2 + (\text{Protein}_{log2})^2}$$

An increase in this distance signifies an escape from homeostatic boundaries, providing a geometric representation of **decanalization**.

---

## 🔗 Data Availability
* **NCBI GEO (A549 Transcriptome)**: [GSE167827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167827) 
* **NCBI GEO (Lung Transcriptome)**: [GSE229537](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229537) 
* **Proteomics and Raw Matrices**: Available in the root directory of this repository.

---
## 📬 Contact

**Dr. Kim-Tuyen Huynh-Dam** Department of Drug Discovery and Biomedical Sciences, University of South Carolina  
**Email**: [kimtuyenhuynhdam@gmail.com](mailto:kimtuyenhuynhdam@gmail.com) | [huynhdam@email.sc.edu](mailto:huynhdam@email.sc.edu)

**Dr. Hippokratis Kiaris** Department of Drug Discovery and Biomedical Sciences, University of South Carolina  
**Email**: [kiarish@cop.sc.edu](mailto:kiarish@cop.sc.edu) 
