#' PACEFit: a fitted PACE model
#'
#' An S4 container for the result of a PACE fit. It holds the raw penalised
#' quasi-likelihood solver state together with the three interpretable output
#' layers: shrunken per-(gene, focal, neighbour) proximity slopes, the per-gene
#' variance decomposition, and the per-pair driver tables. Construct one with
#' [paceFit()]; read it out with [neighbourSlopes()], [varianceDecomposition()],
#' and [topDrivers()].
#'
#' @slot fit The raw streaming PQL solver state (fixed and random effects,
#'   per-cell contamination loadings, gene-wise overdispersions, variance
#'   components).
#' @slot neighbourSlopes A data frame of shrunken proximity coefficients, one row
#'   per gene, focal cell type and neighbour cell type, plus the condition term
#'   for condition cohorts, with the raw estimate, shrunken estimate, and local
#'   false sign rate (lfsr).
#' @slot varianceDecomposition A list of per-gene, per-focal variance
#'   decomposition tables (cell-type identity, spatial cell state, contamination,
#'   residual; and, for condition cohorts, a responder spatial block).
#' @slot topDrivers A list of per-pair driver-score tables ranking the genes that
#'   mediate each focal-neighbour spatial relationship.
#' @slot cellTypes The cell types modelled, in fitting order.
#' @slot params The fitting parameters actually used (bandwidths, contamination
#'   and dispersion settings, condition column, and so on).
#' @slot context The working frame and fixed-effect design retained from the fit
#'   so that the downstream stages ([paceShrink()], [paceDecompose()],
#'   [paceDrivers()]) can run without refitting.
#'
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' fit
#' neighbourSlopes(fit)
#' @name PACEFit
#' @aliases PACEFit-class
#' @rdname PACEFit-class
#' @exportClass PACEFit
setClass(
  "PACEFit",
  slots = c(
    fit                   = "list",
    neighbourSlopes       = "data.frame",
    varianceDecomposition = "list",
    topDrivers            = "list",
    cellTypes             = "character",
    params                = "list",
    context               = "list"
  ),
  prototype = list(
    neighbourSlopes       = data.frame(),
    varianceDecomposition = list(),
    topDrivers            = list(),
    context               = list()
  )
)

## ----------------------------------------------------------------------------
## Validity
## ----------------------------------------------------------------------------
## Without this, replacing a slot wholesale passes silently and breaks every
## readout downstream: `fit@fit <- list()` is the obvious case. The checks are
## structural rather than numerical -- they ask whether the object can still be
## read out, not whether the fit is any good.

## The solver state every downstream layer reads. A fit missing any of these
## cannot be shrunk, decomposed or scored.
.pace_required_fit <- c("B", "U", "se_B", "se_U", "alpha", "re_meta")

## The columns the shrinkage layer is defined by.
.pace_required_slopes <- c("gene", "focal", "neighbour", "term",
                           "estimate", "estimate_shrunk", "lfsr")

setValidity("PACEFit", function(object) {
  problems <- character(0)

  missing_fit <- setdiff(.pace_required_fit, names(object@fit))
  if (length(missing_fit)) {
    problems <- c(problems, paste0(
      "`fit` is missing the solver state the readouts need: ",
      paste(missing_fit, collapse = ", ")))
  }

  if (!length(object@cellTypes) || anyNA(object@cellTypes)) {
    problems <- c(problems, "`cellTypes` must name at least one cell type and cannot contain NA")
  } else if (anyDuplicated(object@cellTypes)) {
    problems <- c(problems, "`cellTypes` must be unique")
  }

  ## The reporting layers are empty until their stage has run, which is valid.
  ## Once populated they have to carry what the readouts promise.
  if (nrow(object@neighbourSlopes)) {
    missing_cols <- setdiff(.pace_required_slopes, names(object@neighbourSlopes))
    if (length(missing_cols)) {
      problems <- c(problems, paste0(
        "`neighbourSlopes` is missing columns: ", paste(missing_cols, collapse = ", ")))
    }
  }

  if (length(object@varianceDecomposition) &&
      !"perGene" %in% names(object@varianceDecomposition)) {
    problems <- c(problems, "`varianceDecomposition` must contain `perGene`")
  }

  if (length(object@topDrivers)) {
    has_scores <- vapply(object@topDrivers,
                         function(x) is.list(x) && "scores" %in% names(x), logical(1))
    if (!all(has_scores)) {
      problems <- c(problems, paste0(
        "every `topDrivers` entry needs a `scores` table; these do not: ",
        paste(names(object@topDrivers)[!has_scores], collapse = ", ")))
    }
  }

  if (length(problems)) problems else TRUE
})
