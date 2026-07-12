# LARItools

LARItools is an R package for transcriptomic analysis of the
lactylation-associated stromal-immune program. It provides functions for LARI
scoring and subtype classification, LacCore signature scoring, and TIDE
analysis.

## Functions

| Function | Purpose | Method |
|---|---|---|
| `calc_LARI()` | Calculate a continuous LARI score | Pre-trained random forest classifier |
| `classify_LARI()` | Assign CS1 or CS2 using the locked discovery-cohort threshold | Fixed cutoff stored in the model bundle |
| `calc_LacCore_score()` | Calculate a LacCore signature score | ssGSEA with GSVA |
| `run_TIDE()` | Calculate TIDE immune dysfunction and exclusion scores | Python `tidepy` through `reticulate` |

## Installation

Install the R dependencies first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("GSVA")
install.packages(c("randomForestSRC", "reticulate", "remotes"))
```

Install LARItools from GitHub:

```r
remotes::install_github("bozh2/LARItools")
```

Alternatively, install a local source package:

```r
install.packages(
  "LARItools_1.2.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

TIDE analysis additionally requires Python and `tidepy`:

```bash
python -m pip install tidepy
```

## Input Format

Expression matrices must have gene symbols as row names and sample identifiers
as column names. For RNA-seq data, use `log2(TPM + 1)` or an equivalent
log-scale expression measure.

```text
              Sample_1  Sample_2  Sample_3
COL5A1           4.21      5.03      3.84
GLT8D2           2.76      3.12      2.45
COL3A1           6.18      7.01      5.92
...
```

LARI normalization uses gene means calculated from the submitted cohort. A
single sample is therefore not a valid input for cohort-level LARI analysis.
Samples that will be compared should be processed together using the same
expression unit and transformation.

## LARI Score

`calc_LARI()` uses the following 10 genes:

```text
COL5A1, GLT8D2, COL3A1, COL1A2, COL1A1,
PDGFRB, OLFML1, FAP, TIMP2, EMILIN1
```

For each gene, expression is centered using the mean of the submitted cohort
and divided by the fixed standard deviation estimated in the discovery
training set:

```text
z = (expression - cohort gene mean) / training-set gene SD
```

The resulting matrix is passed to the bundled random forest model. The class
`"1"` prediction probability is returned as `LARI_score`.

```r
library(LARItools)

expr_log2 <- log2(tpm_matrix + 1)
lari_scores <- calc_LARI(expr_log2)
head(lari_scores)
```

The required model, gene list, training-set scaling parameters, and locked
threshold are stored in:

```text
inst/extdata/lari_model_bundle_v1.2.rds
```

If one or more model genes are absent, `calc_LARI()` issues a warning and
imputes those features with zero. Results from incomplete gene panels should
be interpreted cautiously.

## Locked CS1/CS2 Classification

The default classification threshold is:

```text
0.3180523827398828
```

This value is the median out-of-bag predicted probability for class `"1"`
across 6,798 samples in the TCGA pan-cancer discovery cohort. It is stored in
the model bundle and applied without recalibration to validation cohorts.

```r
lari_result <- classify_LARI(lari_scores)
table(lari_result$LARI_subtype)
```

Classification follows these rules:

```text
LARI_score <= 0.3180523827398828  -> CS1
LARI_score >  0.3180523827398828  -> CS2
```

CS1 represents lower stromal/interstitial activity. CS2 represents higher
stromal/interstitial activity and an immune-excluded phenotype.

For independent validation, use `classify_LARI(lari_scores)` without supplying
a cohort-derived cutoff. The optional `cutoff` argument is provided for
explicit sensitivity analyses and should be reported whenever it is used.

## LacCore Signature Score

`calc_LacCore_score()` calculates an ssGSEA score for the six-gene LacCore
signature and standardizes the score across samples.

```r
lacore_result <- calc_LacCore_score(expr_log2)
head(lacore_result)
```

The default signature is:

```text
CREBBP, EP300, LDHA, PKM, SLC16A1, SLC16A3
```

The exported objects `LACCORE_GENES` and `LRG_PANEL` provide the LacCore genes
and the broader lactylation-related gene panel.

## TIDE Analysis

`run_TIDE()` calls the Python `tidepy` package through `reticulate`.

```r
tide_result <- run_TIDE(
  expr_matrix = tpm_matrix,
  cancer = "LUAD",
  input_type = "raw_counts"
)
```

Supported TIDE cancer mappings are:

| Input | TIDE parameter |
|---|---|
| `SKCM`, `UVM` | `Melanoma` |
| `LUAD`, `LUSC` | `NSCLC` |
| Other values | `Other` |

`input_type` can be `"auto"`, `"raw_counts"`, `"log2"`, or
`"log2_centered"`. A Python interpreter containing `tidepy` can be selected
with `python_path`.

## Model Version

Version 1.2.0 introduced:

- cohort mean centering with fixed discovery-set standard deviations;
- the locked discovery-cohort OOB threshold;
- `classify_LARI()` for reproducible CS1/CS2 assignment;
- a single model bundle containing the model and preprocessing metadata.

## Citation

If you use LARItools, cite the accompanying LARI study and the relevant method
papers:

- Hanzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for
  microarray and RNA-seq data. *BMC Bioinformatics*. 2013;14:7.
  <https://doi.org/10.1186/1471-2105-14-7>
- Barbie DA, et al. Systematic RNA interference reveals that oncogenic
  KRAS-driven cancers require TBK1. *Nature*. 2009;462:108-112.
  <https://doi.org/10.1038/nature08460>
- Jiang P, et al. Signatures of T cell dysfunction and exclusion predict
  cancer immunotherapy response. *Nature Medicine*. 2018;24:1550-1558.
  <https://doi.org/10.1038/s41591-018-0136-1>

The citation for the primary LARI study will be updated when the final
bibliographic record is available.

## License

LARItools is distributed under GPL-2 or GPL-3.

## Contact

Zhihao Bo  
Tianjin Medical University  
Email: bozhihao@tmu.edu.cn  
ORCID: [0009-0005-6298-1953](https://orcid.org/0009-0005-6298-1953)

Issues: <https://github.com/bozh2/LARItools/issues>
