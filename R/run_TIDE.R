# ============================================================
#  Helper: Collect all candidate Python executable paths
# ============================================================
.collect_python_candidates <- function() {
  candidates <- character(0)

  if (.Platform$OS.type == "windows") {

    local_app <- Sys.getenv("LOCALAPPDATA")
    if (nzchar(local_app)) {
      py_root <- file.path(local_app, "Python")
      if (dir.exists(py_root)) {
        subdirs <- list.dirs(py_root, recursive = FALSE, full.names = TRUE)
        for (d in subdirs) {
          exe <- file.path(d, "python.exe")
          if (file.exists(exe)) candidates <- c(candidates, exe)
        }
      }
      py_root2 <- file.path(local_app, "Programs", "Python")
      if (dir.exists(py_root2)) {
        subdirs <- list.dirs(py_root2, recursive = FALSE, full.names = TRUE)
        for (d in subdirs) {
          exe <- file.path(d, "python.exe")
          if (file.exists(exe)) candidates <- c(candidates, exe)
        }
      }
    }

    for (drv in c("C:", "D:")) {
      for (d in Sys.glob(file.path(drv, "Python*"))) {
        exe <- file.path(d, "python.exe")
        if (file.exists(exe)) candidates <- c(candidates, exe)
      }
    }

    conda_roots <- c(
      file.path(Sys.getenv("USERPROFILE"), "miniconda3"),
      file.path(Sys.getenv("USERPROFILE"), "anaconda3"),
      file.path(Sys.getenv("USERPROFILE"), "miniforge3"),
      file.path(Sys.getenv("LOCALAPPDATA"), "miniconda3"),
      file.path(Sys.getenv("LOCALAPPDATA"), "anaconda3"),
      "C:/ProgramData/miniconda3",
      "C:/ProgramData/anaconda3"
    )
    for (root in conda_roots) {
      exe <- file.path(root, "python.exe")
      if (file.exists(exe)) candidates <- c(candidates, exe)
      for (d in list.dirs(file.path(root, "envs"), recursive = FALSE, full.names = TRUE)) {
        exe <- file.path(d, "python.exe")
        if (file.exists(exe)) candidates <- c(candidates, exe)
      }
    }

    path_py  <- Sys.which("python")
    path_py3 <- Sys.which("python3")
    if (nzchar(path_py))  candidates <- c(candidates, path_py)
    if (nzchar(path_py3)) candidates <- c(candidates, path_py3)

  } else {
    # Linux / macOS
    for (cmd in c("python3", "python")) {
      p <- Sys.which(cmd)
      if (nzchar(p)) candidates <- c(candidates, p)
    }
    fixed <- c(
      "/usr/bin/python3",
      "/usr/local/bin/python3",
      "/opt/homebrew/bin/python3",   # macOS Homebrew (Apple Silicon / Intel)
      "/usr/bin/python",
      "/usr/local/bin/python"
    )
    candidates <- c(candidates, fixed[file.exists(fixed)])

    home <- Sys.getenv("HOME")
    for (root in c(
      file.path(home, c("miniconda3", "anaconda3", "miniforge3")),
      "/opt/miniconda3", "/opt/anaconda3", "/opt/miniforge3"
    )) {
      exe <- file.path(root, "bin", "python")
      if (file.exists(exe)) candidates <- c(candidates, exe)
      for (d in list.dirs(file.path(root, "envs"), recursive = FALSE, full.names = TRUE)) {
        exe <- file.path(d, "bin", "python")
        if (file.exists(exe)) candidates <- c(candidates, exe)
      }
    }
  }

  unique(candidates[file.exists(candidates)])
}


# ============================================================
#  Helper: Check whether a given Python has tidepy installed
# ============================================================
.python_has_tidepy <- function(python_exe) {
  tryCatch({
    res    <- system2(python_exe,
                      args   = c("-c", shQuote("import tidepy")),
                      stdout = TRUE, stderr = TRUE)
    status <- attr(res, "status")
    is.null(status) || status == 0
  }, error = function(e) FALSE)
}


# ============================================================
#  Helper: Locate a Python installation that has tidepy
# ============================================================
.find_tidepy_python <- function(python_path = NULL, verbose = TRUE) {
  if (!is.null(python_path)) {
    if (!file.exists(python_path))
      stop("Specified python_path does not exist:\n  ", python_path)
    if (!.python_has_tidepy(python_path))
      stop("Specified Python has no 'tidepy' installed:\n  ", python_path,
           "\n  Run: ", python_path, " -m pip install tidepy")
    return(python_path)
  }

  if (reticulate::py_available(initialize = FALSE)) {
    cur <- reticulate::py_config()$python
    if (.python_has_tidepy(cur)) {
      if (verbose) message("[TIDE] Using reticulate's active Python: ", cur)
      return(cur)
    }
  }

  if (verbose) message("[TIDE] Searching for a Python environment with tidepy...")
  candidates <- .collect_python_candidates()
  if (verbose && length(candidates) > 0)
    message("[TIDE] Found ", length(candidates),
            " Python installation(s); checking for tidepy...")

  for (py in candidates) {
    if (.python_has_tidepy(py)) {
      if (verbose) message("[TIDE] tidepy found in: ", py)
      return(py)
    }
  }

  stop(
    "Could not find a Python environment with 'tidepy' installed.\n\n",
    .cli_box("Option 1: Install tidepy in your current Python",
             "  pip install tidepy"),
    .cli_box("Option 2: Specify the Python path explicitly",
             "  run_TIDE(..., python_path = '/path/to/python')"),
    .cli_box("Option 3: Configure reticulate before calling run_TIDE",
             "  reticulate::use_condaenv('your_env')")
  )
}


# ============================================================
#  Helper: Format a CLI hint box for error messages
# ============================================================
.cli_box <- function(title, ...) {
  lines <- paste0("  ", c(...), collapse = "")
  paste0("\n\u2500\u2500 ", title,
         " ", strrep("\u2500", max(0, 44 - nchar(title))),
         "\n", lines, "\n")
}


# ============================================================
#  Helper: Auto-detect expression matrix format
#  Returns: list(type, reason, stats)
# ============================================================
.detect_expr_format <- function(mat, verbose = TRUE) {

  vals        <- as.numeric(mat)
  vals_finite <- vals[is.finite(vals)]

  v_min    <- min(vals_finite)
  v_max    <- max(vals_finite)
  v_mean   <- mean(vals_finite)
  v_median <- median(vals_finite)
  pct_neg  <- mean(vals_finite < 0)
  pct_gt20 <- mean(vals_finite > 20)
  pct_gt100 <- mean(vals_finite > 100)

  if (verbose) {
    message("[TIDE] \u2500\u2500 Data format detection \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
    message(sprintf("[TIDE]   Range  : [%.2f, %.2f]", v_min, v_max))
    message(sprintf("[TIDE]   Mean   : %.2f  |  Median: %.2f", v_mean, v_median))
    message(sprintf("[TIDE]   %% < 0  : %.1f%%", pct_neg * 100))
    message(sprintf("[TIDE]   %% > 20 : %.1f%%", pct_gt20 * 100))
    message(sprintf("[TIDE]   %% > 100: %.1f%%", pct_gt100 * 100))
  }

  if (pct_gt100 > 0.01) {
    detected <- "raw_counts"
    reason   <- sprintf("%.1f%% of values exceed 100 -- likely raw counts or TPM",
                        pct_gt100 * 100)
  } else if (pct_neg >= 0.10) {
    detected <- "log2_centered"
    reason   <- sprintf("%.1f%% negative values -- likely log2 + mean-centred",
                        pct_neg * 100)
  } else if (v_max <= 20 && v_min >= 0 && v_mean >= 4) {
    detected <- "log2"
    reason   <- sprintf("Range [%.1f, %.1f], mean = %.1f -- typical log2 microarray",
                        v_min, v_max, v_mean)
  } else if (v_max <= 20 && v_min >= 0) {
    detected <- "log2"
    reason   <- sprintf("Range [%.1f, %.1f] -- treated as log2", v_min, v_max)
  } else {
    detected <- "raw_counts"
    reason   <- "Format could not be determined; defaulting to raw_counts"
  }

  if (verbose) {
    message(sprintf("[TIDE]   Detected: %-15s (%s)", detected, reason))
    message("[TIDE] \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  }

  list(type   = detected,
       reason = reason,
       stats  = list(min = v_min, max = v_max,
                     mean = v_mean, median = v_median,
                     pct_neg = pct_neg, pct_gt20 = pct_gt20))
}


# ============================================================
#  Main function: run_TIDE
# ============================================================

#' Run TIDE Immune Dysfunction and Exclusion Analysis
#'
#' @description
#' Calls the Python \pkg{tidepy} package via \pkg{reticulate} to compute TIDE
#' (Tumour Immune Dysfunction and Exclusion) scores for each sample.
#' The function automatically detects the expression format (raw counts,
#' log2, or log2 mean-centred) and applies appropriate pre-processing.
#' Supports Windows, Linux, and macOS.
#'
#' @param expr_matrix A numeric matrix or data.frame with genes as rows
#'   (gene symbols as rownames) and samples as columns.
#' @param cancer Character string specifying the cancer type, e.g. \code{"SKCM"}.
#'   Defaults to \code{"Other"}.
#' @param cancer_map Named character vector mapping TCGA cancer codes to TIDE
#'   cancer labels (\code{"Melanoma"}, \code{"NSCLC"}, or \code{"Other"}).
#' @param input_type Expression scale passed to the TIDE algorithm:
#'   \code{"auto"} (default), \code{"raw_counts"}, \code{"log2"},
#'   or \code{"log2_centered"}.
#' @param pretreat Logical; set to \code{TRUE} if samples are from a
#'   pre-immunotherapy cohort. Defaults to \code{FALSE}.
#' @param vthres Numeric TIDE response threshold. Defaults to \code{0.0}.
#' @param python_path Optional path to a Python executable with tidepy
#'   installed. When \code{NULL} (default), the path is resolved automatically.
#' @param verbose Logical; whether to print progress messages.
#'   Defaults to \code{TRUE}.
#'
#' @return A \code{data.frame} with one row per sample containing all
#'   TIDE output columns plus \code{Sample}, \code{Cancer}, and
#'   \code{TIDE_Cancer_Param}.
#'
#' @importFrom stats median
#' @import reticulate

#' @export
run_TIDE <- function(expr_matrix,
                     cancer      = "Other",
                     cancer_map  = c(SKCM = "Melanoma",
                                     UVM  = "Melanoma",
                                     LUAD = "NSCLC",
                                     LUSC = "NSCLC"),
                     input_type  = c("auto", "raw_counts", "log2", "log2_centered"),
                     pretreat    = FALSE,
                     vthres      = 0.0,
                     python_path = NULL,
                     verbose     = TRUE) {

  input_type <- match.arg(input_type)

  # -- 1. Input validation ------------------------------
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix))
    stop("`expr_matrix` must be a matrix or data.frame.")
  if (is.null(rownames(expr_matrix)))
    stop("`expr_matrix` must have rownames (gene symbols or Entrez IDs).")
  if (ncol(expr_matrix) == 0)
    stop("`expr_matrix` has no samples (columns).")

  # -- 2. Check reticulate availability -----------------
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required.\n",
         "  Install with: install.packages('reticulate')")

  # -- 3. Locate a Python installation with tidepy ------
  python_path <- .find_tidepy_python(python_path, verbose)

  # -- 4. Initialise Python session ---------------------
  reticulate::use_python(python_path, required = TRUE)
  if (verbose) message("[TIDE] Python: ", python_path)

  # -- 5. Import tidepy ----------------------------------
  if (verbose) message("[TIDE] Loading tidepy...")
  reticulate::py_run_string("from tidepy.pred import TIDE as _tide_func")
  TIDE_func <- reticulate::py$`_tide_func`
  pd        <- reticulate::import("pandas")
  if (verbose) message("[TIDE] tidepy loaded successfully.")

  # -- 6. Resolve TIDE cancer parameter -----------------
  tide_cancer <- if (cancer %in% names(cancer_map)) {
    cancer_map[[cancer]]
  } else if (cancer %in% c("Melanoma", "NSCLC", "Other")) {
    cancer
  } else {
    if (verbose) message("[TIDE] Unknown cancer type '", cancer,
                         "'; defaulting to 'Other'.")
    "Other"
  }
  if (verbose) message("[TIDE] Cancer parameter -> \"", tide_cancer, "\"")

  # -- 7. Coerce to numeric matrix -----------------------
  mat       <- as.matrix(expr_matrix)
  mode(mat) <- "numeric"

  # -- 8. Auto-detect or apply user-specified format -----
  if (input_type == "auto") {
    detected   <- .detect_expr_format(mat, verbose)
    input_type <- detected$type
    if (verbose)
      message("[TIDE] Auto-detected input_type = \"", input_type, "\"")
  } else {
    if (verbose)
      message("[TIDE] User-specified input_type = \"", input_type, "\"")
  }

  # -- 9. Pre-processing based on detected format --------
  if (input_type == "raw_counts") {
    force_norm  <- TRUE
    ignore_norm <- FALSE
    if (verbose) message("[TIDE] Normalisation: delegating to TIDE internal",
                         " (log2 + mean-centring)")
  } else if (input_type == "log2") {
    if (verbose) message("[TIDE] Normalisation: subtracting row means (mean-centring)...")
    mat         <- mat - rowMeans(mat, na.rm = TRUE)
    force_norm  <- FALSE
    ignore_norm <- TRUE
    if (verbose) message("[TIDE] Normalisation complete.")
  } else {
    # log2_centered: pass through as-is
    force_norm  <- FALSE
    ignore_norm <- TRUE
    if (verbose) message("[TIDE] Normalisation: data already pre-processed; passing through.")
  }

  # -- 10. Convert R matrix to pandas DataFrame ----------
  if (verbose) message("[TIDE] Converting to pandas DataFrame (",
                       nrow(mat), " genes x ", ncol(mat), " samples)...")
  df_py <- pd$DataFrame(
    mat,
    index   = rownames(mat),
    columns = colnames(mat)
  )

  # -- 11. Call TIDE scoring function --------------------
  if (verbose) message("[TIDE] Running TIDE scoring...")
  t0 <- Sys.time()

  result_py <- TIDE_func(
    df_py,
    cancer          = tide_cancer,
    pretreat        = pretreat,
    vthres          = vthres,
    ignore_norm     = ignore_norm,
    force_normalize = force_norm
  )

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (verbose) message("[TIDE] Scoring completed in ", elapsed, " sec.")

  # -- 12. Convert pandas DataFrame back to R data.frame -
  result_r <- as.data.frame(reticulate::py_to_r(result_py))
  result_r <- cbind(
    Sample            = rownames(result_r),
    result_r,
    Cancer            = cancer,
    TIDE_Cancer_Param = tide_cancer,
    stringsAsFactors  = FALSE
  )
  rownames(result_r) <- NULL

  if (verbose) message("[TIDE] Done. Returning ", nrow(result_r), " sample(s).")
  return(result_r)
}
