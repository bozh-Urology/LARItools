# ============================================================
#  File: R/calc_LacCore_score.R
#  LacCore Signature Score
#  Pipeline: ssgseaParam -> gsva -> t -> scale -> t -> NA->0
#  Requires: GSVA >= 1.48 (new API)
# ============================================================

#' The Canonical Six-Gene LacCore Signature
#'
#' @format A character vector of length 6.
#' @export
LACCORE_GENES <- c("CREBBP", "EP300", "LDHA", "PKM", "SLC16A1", "SLC16A3")


#' The Full Panel of Lactylation-Related Genes (LRGs)
#'
#' @format A named list with five elements.
#' @export
LRG_PANEL <- list(
  glycolytic_enzymes   = c("HK2", "PFKP", "PKM", "ALDOA", "GAPDH",
                           "ENO1", "LDHA", "PGAM1"),
  lactate_transporters = c("SLC16A1", "SLC16A3"),
  pyruvate_oxidation   = c("LDHB", "PDHA1", "PDHB", "DLAT", "DLD"),
  lac_writers          = c("EP300", "CREBBP"),
  lac_erasers          = c("HDAC1", "HDAC2", "HDAC3",
                           "SIRT1", "SIRT2", "SIRT3")
)


#' Compute LacCore Signature Scores via ssGSEA (GSVA new API)
#'
#' @description
#' Calculates per-sample LacCore Signature scores using the official
#' \pkg{GSVA} package (>= 1.48) new API. The pipeline strictly follows:
#'
#' \preformatted{
#'   param     <- ssgseaParam(exprData = mat, geneSets = gene_list)
#'   imm.score <- gsva(param) |> t() |> scale() |> t()
#'   imm.score[is.na(imm.score)] <- 0
#' }
#'
#' The \code{scale()} step performs z-score standardisation across
#' samples for each gene set (centre and scale), which is the standard
#' post-processing used in published immune infiltration analyses.
#'
#' @param expr_matrix A numeric matrix with gene symbols as rownames
#'   and samples as columns. Log2-transformed values (e.g. log2(TPM+1))
#'   are recommended.
#' @param genes Character vector of LacCore gene symbols.
#'   Defaults to \code{LACCORE_GENES}.
#' @param alpha Numeric. ssGSEA weighting exponent. Default \code{0.25}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{data.frame} with one row per sample and two columns:
#'   \describe{
#'     \item{sample}{Character. Sample identifier.}
#'     \item{LacCore_score}{Numeric. Z-score-standardised ssGSEA score.
#'       \code{NA} values are replaced with \code{0}.}
#'   }
#'
#' @references
#' Barbie, D.A. et al. (2009) Systematic RNA interference reveals that
#' oncogenic KRAS-driven cancers require TBK1.
#' \emph{Nature}, 462, 108-112. \doi{10.1038/nature08460}
#'
#' Hanzelmann, S., Castelo, R. & Guinney, J. (2013) GSVA: gene set
#' variation analysis for microarray and RNA-seq data.
#' \emph{BMC Bioinformatics}, 14, 7. \doi{10.1186/1471-2105-14-7}
#'
#' @examples
#' \dontrun{
#' result <- calc_LacCore_score(expr_matrix)
#' head(result)
#' #         sample LacCore_score
#' # 1  SAMPLE_001        1.2341
#' # 2  SAMPLE_002       -0.8823
#' }
#'
#' @importFrom GSVA gsva ssgseaParam
#' @export
calc_LacCore_score <- function(
    expr_matrix,
    genes   = LACCORE_GENES,
    alpha   = 0.25,
    verbose = TRUE
) {
  # -- 0. Check dependency --------------------------------------------------
  if (!requireNamespace("GSVA", quietly = TRUE))
    stop("Please install GSVA: BiocManager::install('GSVA')", call. = FALSE)
  
  # -- 1. Input validation --------------------------------------------------
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix))
    stop("`expr_matrix` must be a matrix or data.frame.", call. = FALSE)
  if (is.null(rownames(expr_matrix)))
    stop("`expr_matrix` must have rownames (gene symbols).", call. = FALSE)
  
  mat       <- as.matrix(expr_matrix)
  mode(mat) <- "numeric"
  
  # -- 2. Gene availability check -------------------------------------------
  missing_g <- setdiff(genes, rownames(mat))
  if (length(missing_g) > 0L)
    warning(length(missing_g), " gene(s) not found and excluded: ",
            paste(missing_g, collapse = ", "), call. = FALSE)
  
  genes_use <- intersect(genes, rownames(mat))
  if (length(genes_use) == 0L)
    stop("No LacCore genes found in expr_matrix.", call. = FALSE)
  
  if (verbose)
    message("[LacCore] Using ", length(genes_use), "/", length(genes),
            " genes: ", paste(genes_use, collapse = ", "))
  
  # -- 3. Build gene set (named list; accepted directly by GSVA new API) ----
  gene_list <- list(LacCore = genes_use)
  
  # -- 4. Core three-line pipeline ------------------------------------------
  #
  #   param     <- ssgseaParam(exprData = mat, geneSets = gene_list)
  #   imm.score <- gsva(param) |> t() |> scale() |> t()
  #   imm.score[is.na(imm.score)] <- 0
  #
  if (verbose)
    message("[LacCore] Calling GSVA::ssgseaParam() + gsva() ...")
  
  param <- GSVA::ssgseaParam(
    exprData = mat,
    geneSets = gene_list,
    alpha    = alpha
  )
  
  imm.score <- GSVA::gsva(param)    # gene-sets x samples matrix
  imm.score <- t(imm.score)         # -> samples x gene-sets
  imm.score <- scale(imm.score)     # z-score per gene-set (base::scale)
  imm.score <- t(imm.score)         # -> gene-sets x samples
  imm.score[is.na(imm.score)] <- 0  # replace NA with 0
  
  # -- 5. Build output data.frame with LacCore_score column only ------------
  result <- data.frame(
    sample        = colnames(mat),
    LacCore_score = as.numeric(imm.score["LacCore", ]),
    row.names     = NULL,
    stringsAsFactors = FALSE
  )
  
  if (verbose)
    message("[LacCore] Done. Score range: [",
            round(min(result$LacCore_score, na.rm = TRUE), 4), ", ",
            round(max(result$LacCore_score, na.rm = TRUE), 4), "]")
  
  return(result)
}