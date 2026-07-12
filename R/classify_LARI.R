#' Classify Samples into LARI Subtypes
#'
#' Classifies samples into CS1 (low stromal activity) or CS2 (high stromal
#' activity) subtypes based on LARI scores computed by \code{calc_LARI()}.
#' Uses a fixed threshold derived from the discovery cohort to ensure
#' consistent, reproducible classification across all cohorts.
#'
#' @param lari_result A data.frame returned by \code{calc_LARI()}, containing
#'   at minimum a column named \code{LARI_score}.
#' @param cutoff Numeric. Classification threshold. If NULL (default), uses
#'   the locked threshold stored in the model bundle
#'   (0.3181, derived as the median OOB-predicted probability in the
#'   6,798-patient TCGA discovery cohort).
#' @param model_path Path to the model bundle RDS file. Only used when
#'   \code{cutoff = NULL}.
#' @param score_col Character. Name of the score column in \code{lari_result}.
#'   Default is \code{"LARI_score"}.
#'
#' @return The input data.frame with two additional columns:
#'   \describe{
#'     \item{LARI_cutoff}{The threshold value used for classification}
#'     \item{LARI_subtype}{Factor with levels CS1 and CS2.
#'       CS1: score <= cutoff (low stromal activity, immune-inflamed).
#'       CS2: score > cutoff (high stromal activity, immune-excluded).}
#'   }
#'
#' @details
#' The locked threshold (0.3180523827) was defined as the median
#' OOB (out-of-bag) predicted probability for class 1 across the 6,798
#' TCGA pan-cancer training samples. This threshold is applied without
#' recalibration to all validation cohorts, ensuring that the classifier
#' is fully independent of validation cohort outcomes.
#'
#' Biological interpretation:
#' CS1 (LARI_score <= 0.3181): Low stromal/interstitial activity.
#'   Associated with immune-inflamed phenotype and better response to
#'   immune checkpoint blockade.
#' CS2 (LARI_score > 0.3181): High stromal/interstitial activity.
#'   Associated with immune-excluded phenotype, TGF-beta signaling,
#'   and reduced response to immune checkpoint blockade.
#'
#' @references
#' Mariathasan et al. (2018) TGFbeta attenuates tumour response to
#' PD-L1 blockade by contributing to exclusion of T cells. Nature.
#'
#' @examples
#' \dontrun{
#' scores <- calc_LARI(expr_matrix)
#' result <- classify_LARI(scores)
#' table(result$LARI_subtype)
#' }
#'
#' @export
classify_LARI <- function(
    lari_result,
    cutoff = NULL,
    model_path = system.file(
      "extdata",
      "lari_model_bundle_v1.2.rds",
      package = "LARItools"
    ),
    score_col = "LARI_score"
) {
  if (!is.data.frame(lari_result)) {
    stop("'lari_result' must be a data.frame (output of calc_LARI()).")
  }
  if (!score_col %in% colnames(lari_result)) {
    stop(
      "Column '", score_col, "' not found in lari_result. ",
      "Did you run calc_LARI() first?"
    )
  }

  scores <- lari_result[[score_col]]
  if (!is.numeric(scores)) {
    stop("Score column '", score_col, "' must be numeric.")
  }

  if (is.null(cutoff)) {
    if (!file.exists(model_path)) {
      stop(
        "Model bundle not found: ", model_path,
        "\nPlease provide a numeric cutoff value manually."
      )
    }
    bundle <- readRDS(model_path)
    if (is.null(bundle$cutoff)) {
      stop(
        "No locked cutoff found in model bundle. ",
        "Please provide a numeric cutoff value manually."
      )
    }
    cutoff <- bundle$cutoff
    message(
      "[classify_LARI] Using locked cutoff: ", round(cutoff, 4),
      " (source: ", bundle$cutoff_method, ")"
    )
  }

  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("'cutoff' must be a single non-missing numeric value.")
  }

  subtypes <- factor(
    ifelse(scores <= cutoff, "CS1", "CS2"),
    levels = c("CS1", "CS2")
  )

  lari_result$LARI_cutoff <- cutoff
  lari_result$LARI_subtype <- subtypes
  lari_result
}
