# maturation score workflow

repo for calculating a pca-based maturation score. basically projecting samples onto a trajectory defined by a control timecourse.

## overview

1.  **find dynamic genes**: run lrt on control group (d21) to find stuff changing over time.
2.  **define trajectory**: pca on those genes. vector is start (day 6) to end (day 17).
3.  **score**: project everyone onto that vector. 0 is early, 1 is late.

## steps

check `tutorial_workflow.Rmd` for the code.

### 1. load & prep

loading metadata and counts. making sure factors are right (day 6 -> 10 -> 17).

```r
# load libs
library(DESeq2)
library(ggplot2)
library(dplyr)

# load data
metadata <- readRDS("dat/metadata/WC24_metadata_clean.rds")
raw_counts <- readRDS("dat/counts/raw/WC24_filt_raw_counts.rds")
vst_counts <- readRDS("dat/counts/vst/WC24_vst_counts.rds")
```

### 2. lrt

finding top 1000 genes changing in d21 controls.

```r
# subset d21
dds <- DESeqDataSetFromMatrix(countData = raw_ctrl, colData = meta_ctrl, design = ~ timepoint)
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

# top 1000
res <- results(dds)
top_genes <- rownames(res[order(res$padj), ])[1:1000]
```

### 3. pca & vector

pca on those genes. defining the vector from day 6 centroid to day 17 centroid.

```r
# centroids
origin <- centroids %>% filter(timepoint == 6)
endpoint <- centroids %>% filter(timepoint == 17)

# vector
ref_vec <- endpoint_vec - origin_vec
```

### 4. plots

#### trajectory

samples in pc space. arrow is the maturation path.

![pca trajectory](plots/pca_trajectory_plot.png)

#### scores

normalized scores (0-1).

**barplot:**
![barplot](plots/maturation_scores_barplot.png)

**trajectory:**
![lineplot](plots/maturation_scores_lineplot.png)

## usage

1.  clone this.
2.  data goes in `dat/`.
3.  run `tutorial_workflow.Rmd`.

## deps
- R
- DESeq2
- ggplot2
- dplyr
- tibble
- tidyr
- scatterplot3d
