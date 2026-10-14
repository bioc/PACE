## ----------------------------------------------------------------------------
## Stashing readouts back into the SpatialExperiment.
## ----------------------------------------------------------------------------
## A PACEFit is deliberately a separate object: a fit on a large cohort runs to
## hundreds of megabytes, it usually outlives the session that made it, and the
## vignette shows stripping and saving one on its own. Folding all of that into
## the SpatialExperiment would tie the two together for good.
##
## What travels well is the readouts, which are small and are what downstream
## work actually reads. stashPACE() puts those in `metadata()`, following the
## convention QC, clustering and deconvolution tools use, and leaves the solver
## state in the fit.

#' Stash PACE readouts in a SpatialExperiment
#'
#' Copies a fit's interpretable output into `metadata(spe)`, so that results
#' travel with the data they came from and downstream code can find them
#' without being handed a second object.
#'
#' The solver state is deliberately not copied. A fit is large, it is often
#' saved and reloaded independently of the data, and nothing downstream reads
#' it; see [fitState()] if you need it. What is stashed is the three reporting
#' layers plus enough provenance to tell which fit produced them: the cell
#' types, the fitting parameters, and the package version.
#'
#' @param spe A [SpatialExperiment::SpatialExperiment] -- the one the fit was
#'   built on. The cell count is checked against the fit.
#' @param fit A [PACEFit].
#' @param name The `metadata()` element to write. Defaults to `"PACE"`.
#' @param overwrite Whether to replace an existing element of that name.
#'
#' @return `spe`, with `metadata(spe)[[name]]` populated.
#'
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' spe <- stashPACE(spe, fit)
#' names(S4Vectors::metadata(spe)$PACE)
#' head(S4Vectors::metadata(spe)$PACE$neighbourSlopes, 3)
#' @seealso [neighbourSlopes()], [varianceDecomposition()], [topDrivers()]
#' @export
stashPACE <- function(spe, fit, name = "PACE", overwrite = FALSE) {
  if (!methods::is(fit, "PACEFit")) {
    stop("`fit` must be a PACEFit; got ", class(fit)[1], ".", call. = FALSE)
  }
  n_fit <- nrow(fit@context$df)
  if (!is.null(n_fit) && ncol(spe) != n_fit) {
    stop("`spe` has ", ncol(spe), " cells but the fit has ", n_fit,
         ". Stash into the SpatialExperiment the fit was built on.", call. = FALSE)
  }
  existing <- S4Vectors::metadata(spe)
  if (!overwrite && name %in% names(existing)) {
    stop("`metadata(spe)$", name, "` already exists. ",
         "Pass overwrite = TRUE to replace it, or choose another `name`.", call. = FALSE)
  }

  existing[[name]] <- list(
    neighbourSlopes       = fit@neighbourSlopes,
    varianceDecomposition = fit@varianceDecomposition,
    topDrivers            = fit@topDrivers,
    cellTypes             = cellTypes(fit),
    params                = paceParams(fit),
    packageVersion        = as.character(utils::packageVersion("PACE")),
    stashedAt             = Sys.time()
  )
  S4Vectors::metadata(spe) <- existing
  spe
}

#' Drop the per-cell matrices from a fit
#'
#' A fit stores the fitted means and the contamination offsets as dense cells by
#' genes matrices. On a real cohort those dominate its size, and nothing
#' downstream needs them kept: [paceDecompose()] and [paceDrivers()] read
#' per-cell-type summaries the fit also holds, and rebuild the matrices from the
#' counts when they are absent.
#'
#' Use this before saving a fit. Every readout still works afterwards, and
#' re-decomposing a compacted fit reproduces the tables it already carries.
#'
#' @param object A [PACEFit].
#' @param which The matrices to drop. The default is all of them.
#'
#' @return `object`, without the named matrices.
#'
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' small <- compactFit(fit)
#' c(before = object.size(fit), after = object.size(small)) |> format(units = "MB")
#' @seealso [stashPACE()] for keeping the readouts with the data.
#' @export
compactFit <- function(object,
                       which = c("mu", "mu_bio", "mu_spill",
                                 "technical_offset_mat", "bleed_offset_mat")) {
  if (!methods::is(object, "PACEFit")) {
    stop("`object` must be a PACEFit; got ", class(object)[1], ".", call. = FALSE)
  }
  which <- match.arg(which, several.ok = TRUE)
  state <- object@fit
  for (nm in which) state[[nm]] <- NULL
  object@fit <- state
  validObject(object)   # the solver state the readouts need must survive
  object
}
