#' Shrunken neighbour proximity slopes
#'
#' Returns the shrunken per-(gene, focal, neighbour) proximity coefficients with
#' their local false sign rate (lfsr). For condition cohorts the responder
#' interaction terms are included.
#'
#' @param object A [PACEFit].
#' @param ... Unused.
#' @return A data frame of shrunken slopes.
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' head(neighbourSlopes(fit))
#' @rdname neighbourSlopes
#' @export
setMethod("neighbourSlopes", "PACEFit", function(object, ...) {
  object@neighbourSlopes
})

#' Per-gene variance decomposition
#'
#' Returns the per-gene, per-focal variance decomposition. By default the
#' observed single-frame decomposition (cell-type identity, spatial cell state,
#' contamination, residual; and a responder spatial block for condition cohorts)
#' is returned; set `which = "blocks"` for the underlying link-scale block tables.
#'
#' @param object A [PACEFit].
#' @param which Either `"perGene"` (default, the observed single-frame table) or
#'   `"blocks"` (the raw link-scale decomposition).
#' @param ... Unused.
#' @return A data frame (for `"perGene"`) or a list (for `"blocks"`).
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' head(varianceDecomposition(fit))
#' @rdname varianceDecomposition
#' @export
setMethod("varianceDecomposition", "PACEFit",
  function(object, which = c("perGene", "blocks"), ...) {
    which <- match.arg(which)
    if (which == "perGene") object@varianceDecomposition$perGene
    else object@varianceDecomposition$blocks
  })

#' Per-pair driver tables
#'
#' Returns the driver-score tables ranking the genes that mediate each
#' focal-neighbour spatial relationship.
#'
#' Each pair's genes are called at lfsr < 0.05. The lfsr is the posterior
#' probability that a call has the wrong sign, so its sum over a pair's calls,
#' `expected_false_sign`, is the expected number of those calls with the wrong
#' sign, and `false_sign_rate` is that count as a fraction of the calls
#' (Stephens 2017). Read a pair's call count against it: many ordered pairs are
#' tested, and a pair with three calls and an expected 1.5 wrong is weak
#' evidence. Both are only as well calibrated as the lfsr itself.
#'
#' On a fit with a condition, the table carries two different slopes and they
#' are not a raw-and-shrunken pair of the same quantity. `b_clean` and `lfsr`
#' are the condition-interaction term, the one the shrinkage was applied to, and
#' the driver score is built from `b_clean`. `u_raw` is the *baseline*
#' slope, the unshrunken BLUP for the neighbour on its own, and it is there
#' because the variance split needs it: `V_S` comes from `u_raw` and `V_RxS`
#' from `b_clean`. On a fit without a condition there is no interaction to
#' separate, so `u_raw` is absent and `b_clean` is the only slope.
#'
#' @param object A [PACEFit].
#' @param ... Unused.
#' @return A named list with one entry per focal-neighbour pair, each a list
#'   with `scores` (the driver table), `status`, `expected_false_sign` and
#'   `false_sign_rate`.
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' names(topDrivers(fit))
#' @rdname topDrivers
#' @export
setMethod("topDrivers", "PACEFit", function(object, ...) {
  object@topDrivers
})

#' @param object A [PACEFit] object.
#' @return The `PACEFit` class itself is a container and is not returned by a
#'   function documented here; construct one with [paceFit()]. The `show` method
#'   is called for its side effect of printing a one-screen summary of the fit --
#'   cell types, kernel bandwidths, contamination and dispersion settings, the
#'   condition term if there is one, and which pipeline stages have run -- and
#'   returns `NULL` invisibly.
#' @importFrom methods show
#' @rdname PACEFit-class
#' @export
setMethod("show", "PACEFit", function(object) {
  p <- object@params
  cat("class: PACEFit\n")
  cat("cell types (", length(object@cellTypes), "): ",
      paste(object@cellTypes, collapse = ", "), "\n", sep = "")
  cat("kernels: h_bio = ", p$h_bio, " um, h_tech = ", p$h_tech,
      " um | contamination: ", p$contamination,
      "; dispersion: ", p$dispersion, "\n", sep = "")
  if (!is.null(p$condition_col))
    cat("condition: ", p$condition_col, " (", p$resp_term, ")\n", sep = "")

  ## pipeline stage status
  done_shrink <- nrow(object@neighbourSlopes) > 0L
  done_decomp <- length(object@varianceDecomposition) > 0L
  done_driv   <- length(object@topDrivers) > 0L
  cat("pipeline: model",
      if (done_shrink) " -> shrink" else "",
      if (done_decomp) " -> decompose" else "",
      if (done_driv)   " -> drivers" else "", "\n", sep = "")
  if (done_shrink) {
    fsr <- expected_false_sign(object@neighbourSlopes$lfsr)
    cat("  neighbour slopes: ", nrow(object@neighbourSlopes),
        " rows (", fsr$n_calls, " at lfsr < 0.05, ",
        sprintf("%.1f", fsr$expected_false_sign),
        " expected to have the wrong sign)\n", sep = "")
  }
  invisible(NULL)
})

#' Per-cell contamination
#'
#' Returns the per-cell contamination loading `rho_i` and the fraction of a
#' cell's expected counts attributed to its local ambient field. Under the
#' per-cell contamination model the expected count is
#' `mu_ig = mu_bio_ig + rho_i * a_ig`, where `a_ig` is the cross-cell-type
#' ambient field at cell `i` and `rho_i` is a single empirical-Bayes shrunken
#' loading shared across genes. The contamination fraction summarises that
#' second term over genes,
#' `contamFraction_i = sum_g mu_spill_ig / sum_g mu_ig`, and is the quantity
#' reported as `contam_frac` in the solver's fitting trace.
#'
#' A high contamination fraction marks a cell whose profile is substantially
#' explained by its neighbours rather than by its own cell type, which is the
#' expected signature of a segmentation or transcript-assignment error. It is
#' not a complete cell-quality score: because `rho_i` carries a fixed
#' gene-direction (the local ambient), it will under-report contamination whose
#' direction departs from the neighbourhood average, and a cell that is wholly a
#' segmentation artefact may fit some incorrect cell type with a low `rho_i`.
#' Pair it with ordinary per-cell quality control.
#'
#' Read `contamFraction` rather than `rho` on its own. The ambient field is
#' cross-cell-type, so a cell with no differently-typed neighbour inside the
#' technical kernel has `a_ig = 0` for every gene and therefore no contamination
#' at all. Its `rho_i` is then unidentified and shrinks to the empirical Bayes
#' prior mean, giving every such cell the same apparently middling loading. On
#' the shipped example fit this is 4,589 of 7,898 cells, 84% of the tumour cells
#' in a tumour-dominated crop. `contamFraction` reports them as zero, which is
#' correct; `rho` alone does not.
#'
#' @param object A [PACEFit] fitted with `contamination = "percell_hc"`.
#' @param ... Unused.
#' @return A data frame with one row per cell, in fitting order: `cell` (the
#'   working-frame cell identifier, or the row index if the frame is unnamed),
#'   `celltype` (the label the cell was fitted under), `rho` (the contamination
#'   loading) and `contamFraction`.
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' cc <- cellContamination(fit)
#' head(cc)
#' # Contamination is highest where cells sit against a different type.
#' tapply(cc$contamFraction, cc$celltype, median)
#' @rdname cellContamination
#' @export
setMethod("cellContamination", "PACEFit", function(object, ...) {
  rho <- object@fit$percell_bleed_rho
  if (is.null(rho))
    stop("no per-cell contamination in this fit; refit with ",
         "contamination = \"percell_hc\".", call. = FALSE)

  ## Both solvers return `contam_frac` directly. Fits made before it was
  ## retained carry the full matrices instead, so recompute from those.
  frac <- object@fit$contam_frac
  if (is.null(frac)) {
    mu_spill <- object@fit$mu_spill
    if (is.null(mu_spill))
      stop("this fit predates `contam_frac` and does not retain `mu_spill`; ",
           "refit to obtain contamination fractions.", call. = FALSE)
    ## `mu` is stored by the dense solver; fall back to the same floor it applies.
    mu <- object@fit$mu
    if (is.null(mu)) mu <- pmax(object@fit$mu_bio + mu_spill, 1e-6)
    frac <- rowSums(mu_spill) / pmax(rowSums(mu), 1e-9)
  }

  df  <- object@context$df
  ids <- if (!is.null(df) && !is.null(rownames(df))) rownames(df) else NULL
  ct  <- if (!is.null(df) && !is.null(df$celltype)) as.character(df$celltype) else NA_character_

  data.frame(
    cell           = if (is.null(ids)) seq_along(rho) else ids,
    celltype       = ct,
    rho            = as.numeric(rho),
    contamFraction = as.numeric(frac),
    stringsAsFactors = FALSE,
    row.names      = NULL
  )
})
