# scImputeBench

A comprehensive benchmarking pipeline for single-cell RNA-seq imputation methods.

## Overview

This repository provides a reproducible pipeline for benchmarking **11 scRNA-seq imputation methods** across **3 datasets** using **13 evaluation metrics**. The pipeline covers the full workflow from data generation and preprocessing through imputation to multi-dimensional evaluation.

## Benchmarked Methods

| Method | Type |
|--------|------|
| RAW (no imputation) | Baseline |
| ALRA | Low-rank approximation |
| DrImpute | Clustering-based |
| MAGIC | Diffusion-based |
| SAVER | Model-based |
| scVI | Deep generative model |
| DCA | Autoencoder |
| scIGANs | GAN-based |
| scMultiGAN | GAN-based |
| scSTD | Diffusion model |
| scIDPMs | Diffusion model |

## Datasets

- **Simulated data** — Generated using the Splatter R package with known ground truth
- **Chu** — Human embryonic stem cell differentiation (Chu et al., 2016)
- **Zheng** — PBMC 68k dataset (Zheng et al., 2017)

## Evaluation Metrics (13 metrics across 3 categories)

### Accuracy
| Metric | Description |
|--------|-------------|
| Cell-level PCC | Pearson correlation between imputed and true expression per cell |
| Pseudobulk AUPRC | AUPRC of pseudobulk-level signal recovery |
| Cell-level AUPRC | AUPRC of cell-level signal recovery |
| logRMSE | Log-transformed root mean squared error (lower is better) |

### Structural Preservation
| Metric | Description |
|--------|-------------|
| ARI | Adjusted Rand Index for clustering agreement |
| NMI | Normalized Mutual Information |
| Purity | Cluster purity score |
| SCC-FC | Spearman correlation of fold changes |

### Biological Consistency
| Metric | Description |
|--------|-------------|
| Specificity | Marker gene specificity |
| Coexpression | Gene co-expression preservation |
| DEG F1 | F1 score for differentially expressed gene detection |
| DEG SCC | Spearman correlation of DEG statistics |
| DEG AUPRC | AUPRC for DEG detection |

## Repository Structure

```
scBenchmark/
├── README.md
├── environment.yml              # Conda environment for the evaluation pipeline
├── notebooks/
│   ├── 01_simulation_data_generation.ipynb
│   ├── 02_zheng_preprocessing.ipynb
│   ├── 03_imputation.ipynb
│   ├── 04_eval_simulation.ipynb
│   ├── 05_eval_chu.ipynb
│   └── 06_eval_zheng.ipynb
├── scripts/
│   └── summary_visualization.py
├── figures/
└── .gitignore
```

## Installation

### Python environment

```bash
conda env create -f environment.yml
conda activate scbenchmark
```

### R packages (required for simulation and some imputation methods)

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("splatter", "DESeq2", "monocle3"))

# CRAN / GitHub packages
install.packages("SAVER")
install.packages("DrImpute")
devtools::install_github("KlugerLab/ALRA")
```

## Usage

Run the notebooks in order:

1. **`01_simulation_data_generation.ipynb`** — Generate simulated scRNA-seq data using Splatter
2. **`02_zheng_preprocessing.ipynb`** — Download and preprocess the Zheng PBMC dataset
3. **`03_imputation.ipynb`** — Run all 11 imputation methods on each dataset
4. **`04_eval_simulation.ipynb`** — Evaluate imputation results on simulated data
5. **`05_eval_chu.ipynb`** — Evaluate imputation results on the Chu dataset
6. **`06_eval_zheng.ipynb`** — Evaluate imputation results on the Zheng dataset
7. **`scripts/summary_visualization.py`** — Generate the summary benchmark visualization

> **Note:** Each notebook defines a `BASE_DIR` variable at the top. Update this path to point to your local data directory before running.

## Results

The summary visualization ranks all 11 methods across 13 metrics. Key findings:

- **DrImpute** achieves the best overall accuracy (highest cell PCC, lowest logRMSE) and strong structural preservation
- **scIGANs** excels at structural preservation (highest ARI, NMI, Purity) and DEG F1
- **MAGIC** shows strong biological consistency (coexpression, specificity)
- **SAVER** provides balanced performance across all categories
- GAN-based and diffusion-based methods show complementary strengths

## Citation

If you use this benchmarking pipeline in your research, please cite:

```
@misc{scbenchmark2025,
  title={scBenchmark: A Comprehensive Benchmarking Pipeline for Single-Cell RNA-seq Imputation Methods},
  year={2025}
}
```

## License

This project is released under the MIT License.
