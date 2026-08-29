# LARItools

LARItools is an R package for cohort-level transcriptomic characterization of a
**lactylation-program-associated stromal-immune tumor state**.

Its primary model is a locked 10-gene mRNA random-forest classifier trained to
distinguish the CS1 and CS2 multi-omics consensus subtypes identified in the
TCGA discovery cohort. `calc_LARI()` returns the predicted probability of the
CS2-associated class, whereas `classify_LARI()` applies a discovery-cohort
threshold to generate predicted **CS1-like** or **CS2-like** labels.

The LARI score is a transcriptomic classifier output. It should not be
interpreted as a direct biochemical measurement of histone lactylation, a
standalone survival probability, or a substitute for direct multi-omics
subtyping.

## Functions

| Function | Purpose | Method |
|---|---|---|
| `calc_LARI()` | Estimate a continuous score for the CS2-associated class | Locked 10-gene mRNA random-forest classifier |
| `classify_LARI()` | Predict a CS1-like or CS2-like class | Fixed discovery-cohort cutoff stored in the model bundle |
| `calc_LacCore_score()` | Calculate a cohort-standardized LacCore signature score | ssGSEA implemented with GSVA |
| `run_TIDE()` | Calculate exploratory TIDE dysfunction and exclusion metrics | Python `tidepy` through `reticulate` |

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
remotes::install_github("bozh-Urology/LARItools")
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

## Input format

Expression matrices must contain gene symbols as row names and sample
identifiers as column names.

For RNA-seq data, use `log2(TPM + 1)` or another comparable log-scale
expression measure. For microarray data, use an appropriately normalized
log-scale expression matrix.

```text
              Sample_1  Sample_2  Sample_3
COL5A1           4.21      5.03      3.84
GLT8D2           2.76      3.12      2.45
COL3A1           6.18      7.01      5.92
...
```

Samples intended for comparison should be processed together using the same
expression unit, normalization procedure, and transformation.

## LARI score

`calc_LARI()` uses the following 10 mRNA features:

```text
COL5A1, GLT8D2, COL3A1, COL1A2, COL1A1,
PDGFRB, OLFML1, FAP, TIMP2, EMILIN1
```

For each model gene, expression is centered using the mean of the submitted
cohort and divided by the fixed standard deviation estimated in the discovery
training set:

```text
z = (expression - submitted-cohort gene mean) / discovery-training gene SD
```

The scaled expression matrix is then passed to the bundled random-forest
classifier. The predicted probability of class `"1"` is returned as
`LARI_score`.

```r
library(LARItools)

expr_log2 <- log2(tpm_matrix + 1)

lari_scores <- calc_LARI(expr_log2)
head(lari_scores)
```

The model, required gene list, discovery-training scaling parameters, and
locked classification threshold are stored in:

```text
inst/extdata/lari_model_bundle_v1.2.rds
```

### Important cohort-level limitation

Because centering is performed using the mean of the submitted cohort, the LARI
score is **cohort-relative**. The same sample can receive a different score if
it is analyzed together with a substantially different sample set.

Therefore:

- single-sample input is not valid for LARI analysis;
- unrelated studies should not be pooled solely for normalization;
- samples being compared should be processed as one analytically coherent
  cohort;
- the locked cutoff removes outcome-based threshold recalibration, but it does
  not make the preprocessing fully single-sample independent.

### Missing model genes

If one or more model genes are absent, `calc_LARI()` issues a warning and
imputes the corresponding standardized feature values with zero. Results from
incomplete gene panels should be interpreted cautiously and the missing genes
should be reported.

## Locked CS1-like/CS2-like classification

The default locked threshold is stored at full precision in the model bundle
and is approximately:

```text
0.3181
```

The underlying full-precision value is:

```text
0.3180523827398828
```

For exact classification, LARI scores should be retained at full numerical
precision internally and rounded only for display.

This threshold was defined as the median out-of-bag predicted probability for
class `"1"` across 6,798 samples in the TCGA discovery training cohort. It is
applied to external cohorts without deriving a new cutoff from their outcomes
or subtype distribution.

```r
lari_result <- classify_LARI(lari_scores)
table(lari_result$LARI_subtype)
```

Classification follows these rules:

```text
LARI_score <= 0.3180523827398828  -> predicted CS1-like
LARI_score >  0.3180523827398828  -> predicted CS2-like
```

The labels produced in external cohorts are classifier-predicted analogues of
the discovery multi-omics consensus subtypes. They are not de novo MOVICS
consensus assignments.

The predicted CS1-like class generally represents lower stromal/interstitial
activity, whereas the predicted CS2-like class represents higher
stromal/interstitial activity and an immune-excluded phenotype.

For independent validation, use:

```r
classify_LARI(lari_scores)
```

without supplying a cohort-derived cutoff. The optional `cutoff` argument is
intended only for explicitly labeled sensitivity analyses and should always be
reported when used.

## LacCore signature score

`calc_LacCore_score()` calculates an ssGSEA score for the six-gene LacCore
signature and then standardizes the score across the submitted samples.

```r
lacore_result <- calc_LacCore_score(expr_log2)
head(lacore_result)
```

The default LacCore signature is:

```text
CREBBP, EP300, LDHA, PKM, SLC16A1, SLC16A3
```

The exported objects `LACCORE_GENES` and `LRG_PANEL` provide the LacCore genes
and the broader lactylation-related gene panel.

### Important interpretation note

Because the ssGSEA result is standardized across the submitted samples, the
reported LacCore score is also cohort-relative. A single-sample input is not
meaningful for the standardized score.

The LacCore score is a transcriptomic surrogate of a lactylation-related
program. It is not a direct measurement of histone lactylation abundance.

## TIDE analysis

`run_TIDE()` calls the Python `tidepy` package through `reticulate`.

Use an `input_type` that matches the actual expression matrix. Do not label TPM
values as raw read counts.

Example for log2-transformed expression:

```r
tide_result <- run_TIDE(
  expr_matrix = expr_log2,
  cancer = "LUAD",
  input_type = "log2"
)
```

Example for already centered log2 expression:

```r
tide_result <- run_TIDE(
  expr_matrix = expr_log2_centered,
  cancer = "LUAD",
  input_type = "log2_centered"
)
```

Supported `input_type` values are:

```text
"auto", "raw_counts", "log2", "log2_centered"
```

A Python interpreter containing `tidepy` can be selected with `python_path`.

Supported cancer mappings are:

| Input | TIDE parameter |
|---|---|
| `SKCM`, `UVM` | `Melanoma` |
| `LUAD`, `LUSC` | `NSCLC` |
| Other values | `Other` |

### Important TIDE limitation

TIDE was originally developed and validated primarily in melanoma and NSCLC
cohorts. Results for other cancer types or treatment settings should be
interpreted as exploratory computational estimates rather than validated
clinical predictions.

## Model version

Version 1.2.0 introduced:

- cohort mean centering with fixed discovery-training standard deviations;
- a locked discovery-cohort out-of-bag threshold;
- `classify_LARI()` for reproducible external-cohort classification;
- a single model bundle containing the classifier and preprocessing metadata.

For reproducible analyses, report:

- the LARItools version;
- the model-bundle version;
- the expression unit and transformation;
- the number of submitted samples;
- any missing model genes;
- whether the default locked threshold or a custom cutoff was used.

## Interpretation summary

LARItools should be interpreted as follows:

- `LARI_score` is the predicted probability of the CS2-associated class from a
  locked 10-gene mRNA classifier;
- predicted CS1-like/CS2-like labels are mRNA-based approximations of the
  discovery multi-omics consensus subtypes;
- LARI and LacCore scores are cohort-relative under the current preprocessing
  implementation;
- neither score directly measures histone lactylation;
- TIDE results outside melanoma and NSCLC should be treated as exploratory.

## Citation

The primary study associated with this package is:

> Bo Z, Zhang S, Zheng Z, Ma J, Feng Y, Li J, Liu Y, Fu Z, Xing H, Yu S, Guo S, Si X, Wang J, Wang R, Wang R, Zhang X, Yue D, Wang Y
> **Pan-cancer multi-omics machine learning defines a lactylation-associated immune-excluded tumor state with proteomic and experimental corroboration**  
> omput Biol Chem. 2026 Jul 31;125:109294. doi: 10.1016/j.compbiolchem.2026.109294. PMID: 42556016.

Please also cite the relevant method papers:

- Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for
  microarray and RNA-seq data. *BMC Bioinformatics*. 2013;14:7.  
  https://doi.org/10.1186/1471-2105-14-7

- Barbie DA, et al. Systematic RNA interference reveals that oncogenic
  KRAS-driven cancers require TBK1. *Nature*. 2009;462:108-112.  
  https://doi.org/10.1038/nature08460

- Jiang P, et al. Signatures of T cell dysfunction and exclusion predict
  cancer immunotherapy response. *Nature Medicine*. 2018;24:1550-1558.  
  https://doi.org/10.1038/s41591-018-0136-1

## License

LARItools is distributed under GPL-2 or GPL-3.

## Contact

Zhihao Bo  
Tianjin Medical University  
Email: bozhihao@tmu.edu.cn  
ORCID: [0009-0005-6298-1953](https://orcid.org/0009-0005-6298-1953)

Issues: <https://github.com/bozh-Urology/LARItools/issues>
