# LOY (Loss of Y Chromosome) Mapping Pipeline

This repository contains scripts and tools to quantify and analyze Loss of Y (LOY) in genomic datasets, primarily from WGS, WES, and RNA-seq data. It is designed for CCLE datasets and other human genomic data sources.

---

## Table of Contents
- [Overview](#overview)  
- [Scripts](#scripts)  
- [Usage](#usage)  
- [Data](#data)  
- [Contributing](#contributing)  
- [License](#license)  

---

## Overview
Loss of Y (LOY) is a common somatic chromosomal aberration with implications in aging and cancer. This pipeline provides tools to:

- Extract reads mapped to chromosome Y.  
- Quantify LOY across samples and tissues.  
- Integrate LOY with RNA-seq for differential expression analysis.  
- Generate visualizations for LOY patterns.

---

## Scripts

### `ExtractYpositions.sh`
Extracts all reads mapped to chromosome Y from a directory of BAM files and records their positions.

**Usage:**
```bash
bash ExtractYpositions.sh /path/to/bam/directory /path/to/output/file.txt
```
### `AnalyzeYPositions.R`
Analyzes chromosome Y read positions to identify sparse regions and low-density genes.  
Generates histograms, cumulative distributions, and visualizations of read density.  

**Dependencies:**  
- R packages: `dplyr`, `tidyr`, `ggplot2`

**Key outputs:**
- `ChromosomeDensityHist.png` – Histogram of read density across chromosome Y  
- `low_density_2000genes.txt` – Genes overlapping low-density regions (below threshold)  
- `ChromosomeCDF.png` – Cumulative distribution function of read densities  
- `orderedSamples.txt` – Sample IDs ordered by number of sparse regions  
- `ChromosomePosition.png` – Scatter plot of reads per sample, ordered by sparsity  
- `ChromosomeBinSummary.png` – Distribution of reads across chromosome Y bins  

**Usage:**
```bash
Rscript AnalyzeYPositions.R
```
## Data

This repository contains data files required for the LOY mapping pipeline.

### Required files

- **`gene_chrY_positions.txt`** – Contains the start and end positions of genes on chromosome Y. 
- **BAM files** – High-coverage WGS/WES BAM files.  
- **RNA-seq counts** – Processed RNA-seq count files for integration with LOY analyses.

> Full CCLE or other controlled-access datasets are not included due to data usage restrictions.


