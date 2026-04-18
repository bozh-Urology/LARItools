# LARItools <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/bozh2/LARItools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bozh2/LARItools/actions/workflows/R-CMD-check.yaml)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0005--6298--1953-green.svg)](https://orcid.org/0009-0005-6298-1953)
<!-- badges: end -->

## Overview

**LARItools** is an R package that provides three complementary scoring
pipelines for transcriptomic data, centred on the lactylation-associated
stromal-immune program. The core metric, the **Lactylation-Associated Risk
Index (LARI)**, is a continuous prognostic score derived from a pre-trained
Random Survival Forest (RSF) classifier. The scoring model implemented in
this package operates on mRNA expression profiles for broad compatibility
with publicly available datasets. The index has been validated as an
independent prognostic indicator across multiple external cohorts spanning
several cancer types.

| Function | Description | Underlying Method |
|---|---|---|
| `calc_LARI()` | Lactylation-Associated Risk Index | Pre-trained Random Survival Forest (RSF) |
| `calc_LacCore_score()` | LacCore signature score | ssGSEA via GSVA ≥ 1.48 new API |
| `run_TIDE()` | Tumour Immune Dysfunction & Exclusion | Python `tidepy` via `reticulate` |

---

## System Requirements

| Requirement | Version | Notes |
|---|---|---|
| R | ≥ 4.1.0 | Required |
| Bioconductor | ≥ 3.14 | For GSVA |
| Python | ≥ 3.7 | Only for `run_TIDE()` |
| tidepy | latest | Only for `run_TIDE()`; install via `pip` |

---

## Installation

### Step 1 — Install Bioconductor dependency (GSVA ≥ 1.48)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA")   # requires GSVA >= 1.48 (new API)
```

### Step 2 — Install CRAN dependencies

```r
install.packages(c("randomForestSRC", "reticulate", "dplyr", "tibble"))
```

### Step 3 — Install LARItools from GitHub

```r
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("bozh-Urology/LARItools")
```

### Step 4 — Python dependency (required only for `run_TIDE()`)

```bash
pip install tidepy
```

> **Note:** `run_TIDE()` requires Python ≥ 3.7 with `tidepy` installed.
> The function auto-detects your Python environment across Windows, Linux,
> and macOS (including Conda / Miniconda / Miniforge). You may also specify
> the interpreter explicitly via `python_path = "/path/to/python"`.

---

## Usage

### `calc_LARI()` — LARI Risk Score

Predicts per-sample LARI scores using a bundled pre-trained RSF model.
The model is built on a gene signature identified in the
accompanying manuscript (Bo *et al.*, in preparation). Missing model
features are imputed with zero and a warning is issued.

```r
library(LARItools)

# expr_matrix: numeric matrix, rownames = gene symbols, columns = samples
# Recommended scale: log2(TPM + 1) or log2(FPKM + 1)
result <- calc_LARI(expr_matrix)
head(result)
#         sample LARI_score
# 1  SAMPLE_001     0.7823
# 2  SAMPLE_002     0.3145
```

> The complete feature gene list and model derivation are described in
> the accompanying manuscript. The bundled model file
> (`inst/extdata/mRNA_model.rds`) contains all information required
> for score computation. See `?calc_LARI` after installation for details.

**Key arguments:**

| Argument | Default | Description |
|---|---|---|
| `expr_matrix` | — | Numeric matrix, genes × samples |
| `model_path` | bundled RDS | Path to pre-trained RSF model |
| `scale` | `TRUE` | Z-score normalisation per feature (within-dataset) |
| `verbose` | `TRUE` | Print progress messages |

---

### `calc_LacCore_score()` — LacCore Signature Score

Computes per-sample LacCore Signature scores using the GSVA ≥ 1.48 new API.
The `scale()` step performs z-score standardisation across samples for each
gene set, consistent with published immune infiltration analyses.

```r
result <- calc_LacCore_score(expr_matrix)
head(result)
#         sample LacCore_score
# 1  SAMPLE_001        1.2341
# 2  SAMPLE_002       -0.8823
```

**Scoring pipeline:**

```
ssgseaParam(exprData, geneSets)
  → gsva()
  → t() → scale() → t()
  → NA → 0
```

> The LacCore gene set and full lactylation-related gene panel (LRG_PANEL)
> are defined in the accompanying manuscript. See `?LACCORE_GENES` and
> `?LRG_PANEL` after installation for the complete gene lists.

**Key arguments:**

| Argument | Default | Description |
|---|---|---|
| `expr_matrix` | — | Numeric matrix, genes × samples |
| `genes` | `LACCORE_GENES` | Gene symbols for the signature |
| `alpha` | `0.25` | ssGSEA weighting exponent |
| `verbose` | `TRUE` | Print progress messages |

---

### `run_TIDE()` — TIDE Immune Dysfunction & Exclusion Score

Calls the Python `tidepy` package via `reticulate` to compute TIDE scores.
Expression format is auto-detected and appropriate pre-processing is applied.

```r
result <- run_TIDE(expr_matrix, cancer = "LUAD")
head(result[, c("Sample", "TIDE", "Dysfunction", "Exclusion", "Responder")])
```

**Cancer type mapping:**

| `cancer` argument | TIDE parameter used |
|---|---|
| `"SKCM"`, `"UVM"` | `"Melanoma"` |
| `"LUAD"`, `"LUSC"` | `"NSCLC"` |
| All others | `"Other"` |

**Input format auto-detection logic:**

| Detected format | Trigger condition | Pre-processing applied |
|---|---|---|
| `raw_counts` | > 1% of values exceed 100 | Delegated to TIDE internal normalisation |
| `log2` | Range [0, 20], mean ≥ 4 | Row-mean centring applied in R |
| `log2_centered` | ≥ 10% negative values | Passed through as-is |

**Key arguments:**

| Argument | Default | Description |
|---|---|---|
| `expr_matrix` | — | Numeric matrix, genes × samples |
| `cancer` | `"Other"` | TCGA code or TIDE label |
| `input_type` | `"auto"` | `"auto"` / `"raw_counts"` / `"log2"` / `"log2_centered"` |
| `pretreat` | `FALSE` | Set `TRUE` for pre-immunotherapy cohorts |
| `vthres` | `0.0` | TIDE response threshold |
| `python_path` | `NULL` | Path to Python executable with tidepy |
| `verbose` | `TRUE` | Print progress messages |

---

## Input Data Requirements

| Requirement | Detail |
|---|---|
| Object class | Numeric `matrix` or `data.frame` |
| Row names | Gene symbols (HGNC), required |
| Column names | Sample identifiers, required |
| Recommended scale | `log2(TPM + 1)` or `log2(FPKM + 1)` |
| Platform | RNA-seq or microarray |
| Missing values | `NA` values are replaced with `0` internally |

---

## Dependencies

### R packages

| Package | Source | Role | Min. version |
|---|---|---|---|
| `GSVA` | Bioconductor | ssGSEA scoring (`ssgseaParam` + `gsva`) | **≥ 1.48** |
| `randomForestSRC` | CRAN | RSF model prediction | ≥ 3.2.0 |
| `reticulate` | CRAN | Python bridge for TIDE | any |
| `dplyr` | CRAN | Data manipulation | ≥ 1.0.0 |
| `tibble` | CRAN | Tidy data frames | ≥ 3.0.0 |
| `stats` | base R | `scale()`, `sd()`, `predict()`, `median()` | — |

### Python (optional — required only for `run_TIDE()`)

| Package | Install command | Role |
|---|---|---|
| `tidepy` | `pip install tidepy` | TIDE immune dysfunction algorithm |
| `pandas` | auto-installed with tidepy | R ↔ Python data exchange |

---

## Citation

If you use **LARItools** in your research, please cite the primary paper:

> Bo Z. *et al.* (2025) Lactylation-Associated Risk Index (LARI) for
> Pan-Cancer Prognosis. *(journal / DOI to be updated upon acceptance)*

Please also cite the following method papers as appropriate:

**ssGSEA:**
> Barbie, D.A. *et al.* (2009) Systematic RNA interference reveals that
> oncogenic KRAS-driven cancers require TBK1.
> *Nature*, **462**, 108–112.
> <https://doi.org/10.1038/nature08460>

**GSVA:**
> Hänzelmann, S., Castelo, R. & Guinney, J. (2013) GSVA: gene set
> variation analysis for microarray and RNA-seq data.
> *BMC Bioinformatics*, **14**, 7.
> <https://doi.org/10.1186/1471-2105-14-7>

**TIDE:**
> Jiang, P. *et al.* (2018) Signatures of T cell dysfunction and exclusion
> predict cancer immunotherapy response.
> *Nature Medicine*, **24**, 1550–1558.
> <https://doi.org/10.1038/s41591-018-0136-1>
>
> tidepy Python package: <http://tide.dfci.harvard.edu>

**Random Survival Forest:**
> Ishwaran, H. *et al.* (2008) Random survival forests.
> *The Annals of Applied Statistics*, **2**(3), 841–860.
> <https://doi.org/10.1214/08-AOAS169>

---

## Availability

This repository is currently **private** pending acceptance of the
accompanying manuscript. The package will be made publicly available
upon publication. Reviewers requiring access should contact the
corresponding author.

---

## License

This package is distributed under the **GNU General Public License,
version 2 or version 3** (GPL-2 | GPL-3).

This licensing choice is required by the inclusion of GPL-licensed
dependencies (`GSVA` [GPL-2] and `randomForestSRC` [GPL ≥ 3]).
Users may choose to comply with either GPL-2 or GPL-3.

- GPL-2: <https://www.gnu.org/licenses/old-licenses/gpl-2.0.html>
- GPL-3: <https://www.gnu.org/licenses/gpl-3.0.html>

© 2025 Zhihao Bo

---

## Contact

**Zhihao Bo**
Tianjin Medical University
📧 bozhihao@tmu.edu.cn
🔗 ORCID: [0009-0005-6298-1953](https://orcid.org/0009-0005-6298-1953)

For bug reports and feature requests, please open an issue at:
<https://github.com/bozh2/LARItools/issues>
