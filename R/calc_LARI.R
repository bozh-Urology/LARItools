#' Calculate LARI Score
#'
#' Computes the LARI (Lactylation-Associated Reactive Interstitial) score
#' for each sample using a pre-trained random forest model. Gene expression
#' values are centered using cohort-level gene means (to correct for additive
#' platform effects) and scaled using fixed standard deviations from the
#' training set (to ensure consistent gene weighting across cohorts).
#'
#' @param expr A matrix or data.frame of gene expression values.
#'   Rows must be genes (with gene symbols as row names),
#'   columns must be samples. Values should be log2(TPM+1) or equivalent.
#' @param model_path Path to the model bundle RDS file. Defaults to the
#'   bundled model in the package.
#' @param verbose Logical. If TRUE, prints diagnostic information.
#'   Default is FALSE.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{sample}{Sample identifiers}
#'     \item{LARI_score}{Continuous LARI score (0 to 1). Higher values
#'       indicate greater stromal activity.}
#'   }
#'
#' @details
#' Normalization strategy:
#' Each sample's expression values are first mean-centered using the
#' cohort-level gene means, which removes additive technical batch effects
#' without utilizing any outcome-related information. The centered values
#' are then scaled using fixed standard deviations derived from the TCGA
#' pan-cancer training set (n=6,798), which preserves the relative
#' importance of each gene as learned during training. This approach is
#' analogous to the median-centering step employed in PAM50 subtype
#' classification.
#'
#' @examples
#' \dontrun{
#' result <- calc_LARI(expr_matrix)
#' head(result)
#' }
#'
#' @importFrom randomForestSRC predict.rfsrc
#' @export
calc_LARI <- function(
    expr,
    model_path = system.file(
      "extdata",
      "lari_model_bundle_v1.2.rds",
      package = "LARItools"
    ),
    verbose = FALSE
) {
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("'expr' must be a matrix or data.frame with genes as rows.")
  }
  if (is.null(rownames(expr))) {
    stop("'expr' must have gene symbols as row names.")
  }
  if (ncol(expr) == 0) {
    stop("'expr' has no samples (columns).")
  }
  if (!file.exists(model_path)) {
    stop(
      "Model bundle not found: ", model_path,
      "\nPlease reinstall LARItools or provide a valid model_path."
    )
  }

  bundle <- readRDS(model_path)

  if (inherits(bundle, "rfsrc")) {
    warning(
      "Legacy model RDS detected. Using default gene list and ",
      "cohort-internal scaling. For locked predictions, please ",
      "use the updated model bundle."
    )
    model <- bundle
    genes <- c(
      "COL5A1", "GLT8D2", "COL3A1", "COL1A2", "COL1A1",
      "PDGFRB", "OLFML1", "FAP", "TIMP2", "EMILIN1"
    )
    train_scale <- NULL
  } else {
    model <- bundle$model
    genes <- bundle$genes
    train_scale <- bundle$train_scale
  }

  missing_genes <- setdiff(genes, rownames(expr))
  if (length(missing_genes) > 0) {
    warning(
      "The following feature genes are missing and will be ",
      "imputed with 0: ",
      paste(missing_genes, collapse = ", ")
    )
    zero_mat <- matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(expr),
      dimnames = list(missing_genes, colnames(expr))
    )
    expr <- rbind(expr, zero_mat)
  }

  sub_mat <- as.matrix(expr[genes, , drop = FALSE])
  input_df <- as.data.frame(t(sub_mat))

  cohort_mean <- colMeans(input_df, na.rm = TRUE)
  input_scaled <- sweep(as.matrix(input_df), 2, cohort_mean, "-")

  if (!is.null(train_scale) && length(train_scale) == length(genes)) {
    input_scaled <- sweep(input_scaled, 2, train_scale, "/")
    if (verbose) {
      cat(
        "[calc_LARI] Normalization: cohort mean centering + ",
        "fixed training-set SD scaling\n"
      )
      cat(
        "[calc_LARI] Mean offset from training center: ",
        round(mean(cohort_mean - bundle$train_center, na.rm = TRUE), 4), "\n"
      )
    }
  } else {
    warning(
      "train_scale not found in bundle. Falling back to ",
      "cohort-internal SD scaling. Results may vary across cohorts."
    )
    cohort_sd <- apply(input_df, 2, stats::sd, na.rm = TRUE)
    cohort_sd[cohort_sd == 0] <- 1
    input_scaled <- sweep(input_scaled, 2, cohort_sd, "/")
  }

  input_scaled_df <- as.data.frame(input_scaled)

  pred_obj <- randomForestSRC::predict.rfsrc(
    model,
    newdata = input_scaled_df
  )

  pred_mat <- pred_obj$predicted
  if (!"1" %in% colnames(pred_mat)) {
    stop(
      "Model output does not contain class '1'. ",
      "Check that the model was trained with binary labels 0/1."
    )
  }
  lari_scores <- as.numeric(pred_mat[, "1"])

  result <- data.frame(
    sample = rownames(input_scaled_df),
    LARI_score = round(lari_scores, 4),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("[calc_LARI] Samples processed:", nrow(result), "\n")
    cat(
      "[calc_LARI] Score range: [",
      round(min(result$LARI_score), 4), ",",
      round(max(result$LARI_score), 4), "]\n"
    )
  }

  return(result)
}
