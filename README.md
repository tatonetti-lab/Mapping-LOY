# LOY (Loss of Y Chromosome) Mapping Pipeline

This repository contains scripts and tools to quantify and analyze Loss of Y (LOY) in genomic datasets, primarily from WGS, WES, and RNA-seq data. It is designed for CCLE datasets and other human genomic data sources.

---

## Table of Contents
- [Overview](#overview)  
- [Repository Structure](#repository-structure)  
- [Installation](#installation)  
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

### `extract_Y_positions.sh`
Extracts all reads mapped to chromosome Y from a directory of BAM files and records their positions.

**Usage:**
```bash
bash scripts/extract_Y_positions.sh /path/to/bam/directory /path/to/output/file.txt

