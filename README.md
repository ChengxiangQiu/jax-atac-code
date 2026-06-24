# Evolutionary transfer learning enables organism-wide inference of mammalian enhancer landscapes

> **Note:** This repository is under active construction and will be updated with additional code and documentation as the manuscripts are finalized.

In this study, we profiled whole-embryo chromatin accessibility across 4 million nuclei from 36 precisely staged mouse embryos at 6-hour temporal resolution, resolving 36 cell classes and 140 cell types (Qiu, Daza, Welsh, *bioRxiv*, 2026). Together with our previously published scRNA-seq atlas (Qiu, Martin, Welsh, *Nature*, 2024), this constitutes the most comprehensive multi-omic resource to date for decoding mammalian embryogenesis. We further developed STEAM-v1, a deep learning framework that augments training with synteny-anchored orthologs from up to 241 mammalian genomes, enabling genome-wide inference of enhancer landscapes across major cell classes in individual species. This repository contains the custom Python and R scripts used for data analysis and model training.

> **Scope:** This repository contains **analysis scripts** specific to this study, not a general-purpose software package. The scripts are provided for transparency and reproducibility. Raw sequencing data are processed by a separate pipeline (see below). All published software dependencies are listed under [System requirements](#1-system-requirements).

## Papers

- *Evolutionary transfer learning enables organism-wide inference of mammalian enhancer landscapes* — [bioRxiv (2026)](https://www.biorxiv.org/content/10.64898/2026.04.07.717039v2.abstract)

## Related repositories

- **Raw-data processing pipeline:** [shendurelab/sciatac_pipeline](https://github.com/shendurelab/sciatac_pipeline) — converts raw BCL / FASTQ to fragment files and peak-by-cell matrices.

## Repository layout

```
jax-atac-code/
├── bin/                              # Computational article creation
├── demo/                             # Two demo scripts: one for dimension reduction, one for model prediction.
├── help_code/                        # Utility scripts shared across the analysis pipeline
├── reproducing_figures/              # Scripts to reproduce each main figure
├── section_1_data_processing/        # Per-sample QC, doublet removal, peak calling
├── section_2_basic_analysis/         # Clustering, annotation, RNA–ATAC integration
├── section_3_naive_model/            # Evolution-naive model training
├── section_4_filtering_windows/      # Selection / filtering of training regions
├── section_5_aware_model/            # Evolution-aware model training
├── section_6_augmented_model/        # Full STEAM-v1 (241-genome ortholog-augmented model)
└── section_7_compare_human_mouse/    # Cross-species enhancer-landscape comparisons
```


---

## 1. System requirements

### Software dependencies and operating systems

**Operating systems tested on:**
- Linux (Ubuntu 22.04 / CentOS 7) — primary

**Programming languages:**
- Python (tested with 3.12.1, 3.12.5, and 3.12.13)
- R (version 4.3.2)

**Command-line tools (raw-data processing):**

| Tool        | Version | Source |
|-------------|---------|--------|
| bcl2fastq   | 2.20    | https://support.illumina.com |
| Trimmomatic | 0.36    | https://github.com/usadellab/trimmomatic |
| Bowtie2     | 2.5.3   | https://github.com/BenLangmead/bowtie2 |
| samtools    | 1.19    | https://github.com/samtools/samtools |
| sambamba    | 0.6.8   | https://github.com/biod/sambamba |
| bedtools    | 2.31.1  | https://github.com/arq5x/bedtools2 |
| AMULET      | 1.1     | https://github.com/UcarLab/AMULET |

**Python packages:**

| Package    | Version | Source |
|------------|---------|--------|
| Scanpy     | 1.10.0  | https://github.com/scverse/scanpy |
| scGLUE     | 0.3.2   | https://github.com/gao-lab/GLUE |
| CREsted    | 1.4.0   | https://github.com/aertslab/crested |
| TF-MINDI   | 1.2.0   | https://github.com/aertslab/TF-MINDI |

**R packages:**

| Package | Version        | Source |
|---------|----------------|--------|
| Seurat  | 5.2.1          | https://github.com/satijalab/seurat |
| Signac  | 1.9.0 & 1.14.0 | https://github.com/stuart-lab/signac |
| BPCells | 0.1.0          | https://github.com/bnprks/BPCells |

### Versions the software has been tested on

All scripts were developed and tested with the exact versions listed above. Newer versions of most packages should be compatible, but Seurat / Signac in particular have breaking API changes across major releases — we recommend pinning to the listed versions when reproducing our analyses.

### Required non-standard hardware

- **GPU recommended** for STEAM-v1 / CREsted / TF-MINDI deep-learning steps. Models were trained across NVIDIA A100, L40, and H200 GPUs. Any modern NVIDIA GPU with ≥ 40 GB VRAM should be sufficient.
- For full-dataset analyses we recommend ≥ 100 GB RAM. Subset and demo analyses run on a standard 16 GB laptop.
- No other non-standard hardware is required.

---

## 2. Installation guide

### Instructions

We recommend isolated environments for Python and R.

**Clone the repository:**
```bash
git clone https://github.com/ChengxiangQiu/jax-atac-code.git
cd jax-atac-code
```

**Python environment (conda):**
```bash
conda create -n jax-atac python=3.12
conda activate jax-atac
pip install scanpy==1.10.0
pip install scglue==0.3.2
pip install crested==1.4.0
pip install tf-mindi==1.2.0
# AMULET: follow instructions at https://github.com/UcarLab/AMULET
```

**R environment:**
```r
install.packages("Seurat")                       # 5.2.1
install.packages("Signac")                       # 1.9.0 or 1.14.0
remotes::install_github("bnprks/BPCells/r")      # 0.1.0
```

**Command-line tools** (bcl2fastq, Trimmomatic, Bowtie2, samtools, sambamba, bedtools) can be installed via your system package manager or `conda install -c bioconda ...`.

### Typical install time on a "normal" desktop computer

Approximately **30–60 minutes** on a standard desktop with a good internet connection, dominated by compiling deep-learning dependencies (TensorFlow / PyTorch).

---

## 3. Demo

A small demo dataset and walkthrough are provided in `demo/`.

### Instructions to run on data

```bash
# Example: performing dimension reduction on a subset of the kidney trajectory
Rscript demo/Demo_1_dimension_reduction.R

# Example: model prediction using STEAM-v1
python demo/Demo_2_model_prediction.py
```

### Expected output

**Demo 1 (dimension reduction):**
- `pca_coor.csv`, `umap_coor.csv`, `clustering_result.csv`

**Demo 2 (model prediction):**
- `demo_predictions.csv`

### Expected run time for demo on a "normal" desktop computer

Approximately 5–15 minutes on a 16 GB laptop. Demo 1 is CPU-only. Demo 2 (model prediction) runs on CPU but is ~5–10× faster on a GPU.

---

## 4. Instructions for use

### How to run the software on your data

The analyses in this study were carried out in the following order. To apply the same workflow to your own sci-ATAC-seq3 data:

1. **Raw-data processing** — Run the pipeline at [shendurelab/sciatac_pipeline](https://github.com/shendurelab/sciatac_pipeline) on your raw BCL / FASTQ files to obtain fragment files and a peak-by-cell matrix.
2. **Section 1 (`section_1_data_processing/`)** — Per-sample QC, doublet removal, peak calling, and matrix construction.
3. **Section 2 (`section_2_basic_analysis/`)** — Clustering, cell-type annotation, and integration with paired scRNA-seq data via scGLUE.
4. **Sections 3–6 (`section_3_naive_model/` → `section_6_augmented_model/`)** — Progressive training of STEAM-v1 model variants (Evolution-naive → Evolution-aware → STEAM-v1).
5. **Section 7 (`section_7_compare_human_mouse/`)** — Cross-species comparison of inferred enhancer landscapes.
6. **Figure reproduction (`reproducing_figures/`)** — Regenerate manuscript figure panels.
7. **Demo (`demo/`)** — Two demo scripts: one for dimension reduction, one for model prediction.

Each script's header documents its expected inputs and outputs. See the **Methods** section of the manuscript for full parameter justifications and biological rationale.

### (Optional) Reproduction instructions

To reproduce the figures and quantitative results from the manuscript:

1. Download the processed data from **GEO: GSE325776**.
2. Update file paths to point to your local copy.
3. Run the scripts in `reproducing_figures/` in numerical order.

---

## Data

Large processed data files (peak-by-cell matrices, model checkpoints, genome-wide predictions) are hosted at **GEO: GSE325776** and https://chengxiangqiu.github.io/jax-atac-website/. Individual scripts document which files they consume.

## License

Released under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

## Citation

If you use this code, please cite:

> Qiu, Daza, Welsh, et al. Evolutionary transfer learning enables organism-wide inference of mammalian enhancer landscapes. *bioRxiv* (2026). [link](https://www.biorxiv.org/content/10.64898/2026.04.07.717039v2.abstract)

## Contact

For questions, please open a GitHub issue or contact Chengxiang Qiu (Chengxiang.Qiu [at] dartmouth.edu).
