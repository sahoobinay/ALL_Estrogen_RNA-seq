# ERα RNA-seq Analysis Pipeline

RNA-seq analysis pipeline for the Nalm6 estrogen receptor (ERα) project.  
Covers the full workflow from Salmon gene counts to publication-ready figures.

**Conditions:** Vehicle (Veh) · 17β-Estradiol (E2) · Tamoxifen (Tam) · Tamoxifen + E2 (TamE2)  
**Replicates:** n = 2 per condition (r1, r2)  
**Cell line:** Nalm6 (B-cell precursor)

---

## Pipeline overview

| Section | Analysis |
|---|---|
| 1–4 | Data loading, edgeR QL model (`~ rep + condition`) |
| 5 | ENSEMBL → gene symbol mapping (versioned IDs handled) |
| 5B | Ribosomal (RPL\*, RPS\*) and mitochondrial (MT-\*) gene removal |
| 6 | Batch correction for visualization (`removeBatchEffect`) |
| 7 | QC figures: library size, PCA (before/after correction), correlation heatmap |
| 8 | Differential expression: E2 vs Veh, Tam vs Veh, TamE2 vs Veh |
| 9 | Volcano plots + effect size distribution |
| 10 | Heatmap of top 100 DE genes |
| 11 | GSEA — MSigDB Hallmark gene sets (fgsea) |
| 12 | UpSet plot — DE gene overlap across contrasts |
| 13–14 | TF activity (DoRothEA + VIPER) + network centrality |
| 15–17 | Stemness scoring + marker expression vs Vehicle |
| 18 | Cell state axis (Proliferation vs Persistence) |
| 19 | Integrated program correlation (Stemness × Persistence) |

---

## Repository structure

```
rnaseq_pipeline/
├── pipeline.R              # Main analysis script
├── install_packages.R      # R package installer
├── session_info.txt        # R session info for reproducibility
├── data/
│   └── salmon.merged.gene_counts.tsv   # Input (not included — see below)
├── results/
│   ├── figures/            # All output figures (generated on run)
│   └── tables/             # All output tables (generated on run)
└── README.md
```

---

## Input data

Place your Salmon gene count matrix at:

```
data/salmon.merged.gene_counts.tsv
```

The file should have genes as rows (versioned ENSEMBL IDs, e.g. `ENSG00000000419.14`)  
and samples as columns matching the names in the metadata section of `pipeline.R`.

> Raw data are not included in this repository. Contact the corresponding author for access.

---

## Installation

### 1. Install R packages

```r
source("install_packages.R")
```

### 2. Run the pipeline

```r
source("pipeline.R")
```

Or from the terminal:

```bash
Rscript pipeline.R
```

---

## Key design decisions

**Batch correction strategy:**  
Replicate (r1 / r2) corresponds to library preparation batch. The design matrix includes `rep` as a covariate (`~ rep + condition`) so that differential expression testing statistically accounts for batch. A separate batch-corrected expression matrix (`logcpm_bc`) is computed with `limma::removeBatchEffect` for visualization only — it is never used for DE testing.

**Ribosomal/mitochondrial filtering:**  
Ribosomal protein genes (RPL\*, RPS\*) and mitochondrial genes (MT-\*) are removed from the visualization matrix to prevent them from dominating PCA and heatmap variance.

**Gene symbol deduplication for GSEA:**  
Where multiple ENSEMBL IDs map to the same gene symbol, the entry with the highest absolute log2FC is retained to ensure a non-redundant rank vector.

**Namespace conflicts:**  
`igraph`, `viper`, and `IRanges` all export functions that conflict with `dplyr` and `base`. Canonical assignments are made at the top of the script (e.g. `select <- dplyr::select`, `intersect <- base::intersect`).

---

## Output

All figures are saved to `results/figures/` and all tables to `results/tables/`.

| Figure | Description |
|---|---|
| Fig1A–D | QC: library size, PCA (×2), sample correlation |
| Fig2A–B | Volcano plot + effect size histogram (E2 vs Veh) |
| Fig3 | Heatmap: top 100 DE genes (E2 vs Veh) |
| Fig4 | GSEA Hallmark barplot (E2 vs Veh) |
| Fig5 | UpSet: DE overlap across contrasts |
| Fig6–7 | TF activity heatmaps (VIPER) |
| Fig8 | Network centrality — top hub TFs |
| Fig9 | Composite stemness score per condition |
| Fig10 | Stemness marker heatmap |
| Fig11 | Stemness marker Δlog-CPM dotplot vs Vehicle |
| Fig12 | Stemness gene set activity barplot vs Vehicle |
| Fig13 | Per-gene stemness marker boxplots |
| Fig14 | Cell state axis (Proliferation vs Persistence) |
| Fig15–16 | Program correlation heatmap + scatter |

---

## Citation

If you use this pipeline, please cite the relevant tools:

- **edgeR**: Robinson MD et al. *Bioinformatics* 2010
- **limma/voom**: Ritchie ME et al. *Nucleic Acids Res* 2015
- **fgsea**: Korotkevich G et al. *bioRxiv* 2021
- **DoRothEA**: Garcia-Alonso L et al. *Genome Res* 2019
- **VIPER**: Alvarez MJ et al. *Nat Genet* 2016
- **MSigDB Hallmark**: Liberzon A et al. *Cell Syst* 2015

---

## Contact

Bikash Sahoo — UVM Larner College of Medicine
