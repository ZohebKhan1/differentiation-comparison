# Principal Component Analysis-based Differentiation Progression Comparison Workflow

This repository is an example workflow for a PCA-based method for assessing the progression of differentiation between two sets of samples (bulk RNA-seq). 

The basic overview is as follows:

- You have a **baseline differentiation time course** (at least 3 timepoints).
- You have one or more **experimental groups** (genotype, treatment, etc.).
- You use the baseline time course to define a **reference trajectory** in PCA space.
- You project all samples onto that trajectory to obtain a **maturation score** that reflects how far along differentiation each sample is.

For the following workflow, I used publicly available data that I processed from the  reference below. This dataset has 3 timepoints: Day 6, Day 10, and Day 17. The authors transcriptionally profiled two sets of samples: D21 (which is our control/baseline for this workflow) and T21, which is the perturbation/experimental group that we will be projecting onto the control principal components space.

Martinez JL, Piciw JG, Crockett M, Sorci IA, Makwana N, Sirois CL, Giffin-Rao Y, Bhattacharyya A. Transcriptional consequences of trisomy 21 on neural induction. Front Cell Neurosci. 2024 Jan 30;18:1341141. doi: 10.3389/fncel.2024.1341141. PMID: 38357436; PMCID: PMC10865501.

### NOTE

- D21 refers to control that we use to define baseline differentiation trajectory
- T21 refers to experimental group that we project **onto** the control differentiation trajectory

## Method Overview

1.  **Find temporally changing genes in the control time course**  
   Use DESeq2’s likelihood ratio test (LRT) on the D21 samples only (Day 6 → Day 10 → Day 17). We then select the top 1,000 genes with the most significant timepoint effect. These genes define the main transcriptional changes during baseline differentiation (which will we will compare experimental group to)
   - LRT reference: [DESeq2 LRT tutorial](https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html) 

2. **Perform PCA on dynamic genes**  Run PCA on DESeq2-calculated variance-stabilized transformed count of the top 1000 LRT genes on just D21 samples. 

3. **Define a differentiation trajectory in PC space**  
   Using only the baseline control differentiation, we:
   - Compute the centroid of Day 6 (earliest timepoint)samples.
   - Compute the centroid of Day 17 samples (latest timepoint).
   - Define a **reference vector in PC space** from the Day 6 centroid to the Day 17 centroid . This is essentially our 'differentiation trajectory' that we will project experimental samples onto and compare with.

4. **Compute a PC-based maturation score each every sample**  
   For each sample (control + experimental):
   - Take its coordinates on the PC1, PC2, and PC3 (to capture largest % variation).
   - Project that sample's position onto the control Day 6 → Day 17 vector we defined earlier .  
   - Normalize this projection in such a way that:
     - ~0 ≈ Day 6 centroid control,
     - ~1 ≈ Day 17 centroid control.  
     - So a maturation score of 1 indicates the differentiation at the end of the control timecourse. We will thus be able to compare in the experimental group how far long the T21 samples are relative to thise control trajectory.

## Files

- `tutorial_workflow.Rmd` – fullly annotated workflow for the maturation score calculation
- `dat/` –  location for input data
  - `dat/metadata/WC24_metadata_clean.rds`
  - `dat/counts/raw/WC24_filt_raw_counts.rds`
  - `dat/counts/vst/WC24_vst_counts.rds`
- `plots/` – directory where figures will be saved
- `results/` – directory where the maturation score table will be written

### 1. Load data & dependencies

Loading metadata and counts, making sure factors are correct (Day 6 -> 10 -> 17)

```r
# Load libs
library(DESeq2)
library(ggplot2)
library(dplyr)

# Load data
metadata <- readRDS("dat/metadata/WC24_metadata_clean.rds")
raw_counts <- readRDS("dat/counts/raw/WC24_filt_raw_counts.rds")
vst_counts <- readRDS("dat/counts/vst/WC24_vst_counts.rds")
```

### 2. Likelihood Ratio Test (LRT) from DESeq2

Using DESeq2 LRT to find top 1000 genes changing in control differentiation (temporal genes)

```r

# Subset D21
dds <- DESeqDataSetFromMatrix(countData = raw_ctrl, colData = meta_ctrl, design = ~ timepoint)
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

# Top 1000
res <- results(dds)

# we define top 1000 by padj because LRT does not output log2FoldChange values since we're just running it on one control set of samples

top_genes <- rownames(res[order(res$padj), ])[1:1000]

```

### 3. PCA & Vector

Using top 1000 DESEq2 genes to define the PC space and also the vector for 'control differentiation trajectory' from Day 6 centroid to Day 17 centroid of the control samples.

```r
# Centroids
origin <- centroids %>% filter(timepoint == 6)
endpoint <- centroids %>% filter(timepoint == 17)

# Vector
ref_vec <- endpoint_vec - origin_vec
```

### 4. Plots

#### Trajectory

Here is a visualization of the control samples in PC space for PC1/PC2/PC3. The arrow represets the differentiation 'path' or maturation 'path' that control samples follow as they differentiate. This will be our baseline to compare to the experimental group..

<img src="plots/pca_trajectory_plot.png" width="500">

#### Scores

Normalized scores (0-1) per sample based on their projection onto the control differentiation trajectory

Based on these results, it looks like the experimental group (T21) differentiates faster initially compared to the reference differentiation, but at the latest timepoint, the control samples differntiate faster than the T21 samples.

Specifically, since we defined the control differentiation trajectory as the baseline, we can state that:
- T21 samples are 13.3% further ahead in differentiation compared to D21 Day 6
- T21 samples are 10.9% further ahead in differentiation compared to D21 Day 10
- T21 samples are 9.9% behind the differentiation compared to D21 Day 17

**Barplot:**
<br>
<img src="plots/maturation_scores_barplot.png" width="400">

**Trajectory:**
<br>
<img src="plots/maturation_scores_lineplot.png" width="400">

## Dependencies for this workflow
- Performed in R
- DESeq2 (LRT)
- ggplot2
- dplyr
- tibble
- tidyr
- scatterplot3d
