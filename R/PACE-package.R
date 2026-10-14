#' PACE: Proximity-Associated Changes in Expression
#'
#' PACE quantifies how a cell's gene expression changes with proximity to specific
#' neighbouring cell types in imaging-based spatial transcriptomics. It fits
#' hierarchical negative binomial mixed models with partial pooling across cell
#' types, corrects for transcript contamination between adjacent cells with a
#' per-cell ambient term, and decomposes expression variance into cell-type
#' identity, spatial cell state, contamination, and residual components.
#'
#' The user-facing entry point is [paceFit()], which takes a
#' [SpatialExperiment::SpatialExperiment] and returns a [PACEFit] object. The
#' results are read out with [neighbourSlopes()], [varianceDecomposition()], and
#' [topDrivers()].
#'
#' @keywords internal
#' @name PACE-package
#' @aliases PACE
#'
#' @useDynLib PACE, .registration = TRUE
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom Rcpp evalCpp
#' @import methods
#' @import dplyr
#' @import ggplot2
#' @import tibble
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom mashr mash mash_set_data cov_canonical cov_ed cov_pca get_significant_results estimate_null_correlation_simple
#' @importFrom ashr get_lfsr get_pm get_psd
#' @importFrom stats var sd median quantile pnorm qnorm optimize coef predict
#'   as.formula model.matrix setNames rnorm
#' @importFrom methods new validObject is
"_PACKAGE"

## Column names used inside dplyr/tidyr pipelines. Declaring them keeps
## `R CMD check` from reading each one as an undefined global.
utils::globalVariables(c(
  "Block", "Cell type %", "driver_score", "driver_score_spec4", "R2_RxS", "R2_S", "Residual %",
  "Responder spatial state %", "Spatial state %", "Spillover %", "Total",
  "V_RxS", "V_S", "V_resid", "V_total", "b_clean", "block",
  "celltype_offset_sq", "delta", "drop_patient", "estimate", "estimate_shrunk",
  "focal", "focal_mean", "gene", "group", "is_contaminated", "level", "lfsr",
  "mu_bar", "n_focal", "neighbour", "pct", "pct_total", "per_gene",
  "resp_pct_drop", "scaled_estimate", "sd_shrunk", "spec", "spec_w_sum",
  "std.error", "term", "total_4", "u_raw"
))
