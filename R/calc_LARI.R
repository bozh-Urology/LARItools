#' Calculate LARI Risk Score
#'
#' @description
#' Computes per-sample LARI (Lung Adenocarcinoma Risk Index) risk scores
#' from a gene expression matrix using a pre-trained random forest model.
#' Missing model features are imputed with zero and a warning is issued.
#'
#' @param expr_matrix A matrix or data.frame of gene expression values,
#'   with rownames as gene symbols and columns as samples.
#' @param model_path Path to the trained LARI model RDS file. Defaults to
#'   the bundled model in the \pkg{LARItools} package.
#' @param scale Logical. Whether to apply z-score normalisation
#'   (within-dataset, per feature). Default \code{TRUE}.
#' @param verbose Logical. Whether to print progress messages.
#'   Default \code{TRUE}.
#'
#' @return A \code{data.frame} with two columns: \code{sample} and
#'   \code{LARI_score} (rounded to 4 decimal places).
#'
#' @importFrom stats sd predict median
#' @export
calc_LARI <- function(expr_matrix,
                      model_path = system.file("extdata", "mRNA_model.rds",
                                               package = "LARItools"),
                      scale   = TRUE,
                      verbose = TRUE) {
  
  # ── Feature gene set ──────────────────────────────────
  LARI_FEATURES <- c("COL5A1", "GLT8D2", "COL3A1", "COL1A2", "COL1A1",
                     "PDGFRB", "OLFML1", "FAP",    "TIMP2",  "EMILIN1")
  
  # ── Input validation ──────────────────────────────────
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix))
    stop("`expr_matrix` must be a matrix or data.frame.")
  if (is.null(rownames(expr_matrix)))
    stop("`expr_matrix` must have rownames (gene symbols).")
  if (ncol(expr_matrix) == 0)
    stop("`expr_matrix` has no samples (columns).")
  if (!file.exists(model_path))
    stop("Model file not found: ", model_path)
  
  # ── Load model ────────────────────────────────────────
  if (verbose) message("[LARItools] Loading LARI model...")
  model <- readRDS(model_path)
  
  # ── Handle missing features ───────────────────────────
  missing_genes <- setdiff(LARI_FEATURES, rownames(expr_matrix))
  if (length(missing_genes) > 0) {
    warning(
      "The following feature gene(s) were not found in `expr_matrix` ",
      "and have been imputed with 0: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
    zero_mat <- matrix(
      0,
      nrow     = length(missing_genes),
      ncol     = ncol(expr_matrix),
      dimnames = list(missing_genes, colnames(expr_matrix))
    )
    expr_matrix <- rbind(expr_matrix, zero_mat)
  }
  
  # ── Subset and transpose to sample x feature ─────────
  input_df <- as.data.frame(t(expr_matrix[LARI_FEATURES, , drop = FALSE]))
  
  # ── Optional z-score normalisation ───────────────────
  if (scale) {
    if (verbose)
      message("[LARItools] Applying z-score normalisation (within-dataset, per feature)...")
    feat_mean             <- colMeans(input_df)
    feat_sd               <- apply(input_df, 2, sd)
    feat_sd[feat_sd == 0] <- 1   # avoid division by zero for constant features
    input_df <- as.data.frame(
      sweep(sweep(input_df, 2, feat_mean, "-"), 2, feat_sd, "/")
    )
  }
  
  # ── Predict LARI scores ───────────────────────────────
  if (verbose)
    message("[LARItools] Predicting LARI scores for ",
            nrow(input_df), " sample(s)...")
  
  pred_obj   <- predict(model, newdata = input_df)
  pred_prob  <- pred_obj$predicted
  risk_score <- pred_prob[, "1"]
  
  # ── Resolve sample identifiers ────────────────────────
  sample_ids <- colnames(expr_matrix)
  if (is.null(sample_ids) || all(nchar(trimws(sample_ids)) == 0)) {
    sample_ids <- paste0("Sample_", seq_len(nrow(input_df)))
  }
  
  # ── Assemble result data.frame ────────────────────────
  result <- data.frame(
    sample     = sample_ids,          # ✅ sample_id -> sample
    LARI_score = round(as.numeric(risk_score), 4),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  
  if (verbose)
    message("[LARItools] Done. Returning ", nrow(result), " sample(s).")
  
  return(result)
}