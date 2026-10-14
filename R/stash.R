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
