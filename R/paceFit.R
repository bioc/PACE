## ---------------------------------------------------------------------------
## Shared SPE -> (counts, cell frame) extraction used by the fit entry points.
## ---------------------------------------------------------------------------
.pace_prepare <- function(object, celltype_col, image_col, condition_col,
                          assay_name, resp_term) {
  if (!assay_name %in% SummarizedExperiment::assayNames(object))
    stop("assay '", assay_name, "' not found in the object.", call. = FALSE)
  cd <- as.data.frame(SummarizedExperiment::colData(object))
  if (!celltype_col %in% colnames(cd))
    stop("celltype_col '", celltype_col, "' not found in colData.", call. = FALSE)

  ## counts are genes x cells in an SPE; PACE wants cells x genes, and keeps
  ## them sparse all the way into the solver.
  Y <- .pace_as_dgc(Matrix::t(SummarizedExperiment::assay(object, assay_name)))
  ## Feature names index every reported table, so duplicates would make two
  ## different genes share one row of the statistics and one row of the output.
  if (anyDuplicated(colnames(Y)))
    stop("the assay has duplicated feature names (",
         paste(utils::head(unique(colnames(Y)[duplicated(colnames(Y))]), 3), collapse = ", "),
         "); make them unique before fitting.", call. = FALSE)
  ## a zero library size makes the log offset and log1p CP10k undefined
  empty_cells <- sum(Matrix::rowSums(Y) == 0)
  if (empty_cells > 0L)
    stop(empty_cells, " cell(s) have zero total counts; remove them before fitting.",
         call. = FALSE)
  ## exactly two finite coordinate columns, refused (not dropped) otherwise
  coords <- .pace_coordinate_matrix(SpatialExperiment::spatialCoords(object))
  df <- cd
  df[["x"]] <- coords[, 1]
  df[["y"]] <- coords[, 2]

  ## single-section objects get a constant image grouping.
  if (is.null(image_col)) {
    df[[".image"]] <- factor("all")
    image_col <- ".image"
  }
  ## derive the responder term prefix for condition cohorts, and set the 0/1
  ## condition indicator the decomposition and driver scoring require.
  if (!is.null(condition_col)) {
    ## model.matrix names the term paste0(condition_col, level), so both the
    ## default term and the level it refers to must come from condition_col.
    ## Hardcoding "Responder" here left .resp_dummy all zeros for any other
    ## column name, silently zeroing the whole responder block.
    lev        <- levels(factor(df[[condition_col]]))
    cond_terms <- paste0(condition_col, lev)
    if (is.null(resp_term)) resp_term <- cond_terms[length(cond_terms)]
    hit <- match(resp_term, cond_terms)
    if (is.na(hit)) {
      stop(sprintf(
        "resp_term '%s' does not name a level of condition_col '%s'. Valid terms: %s.",
        resp_term, condition_col, paste(cond_terms, collapse = ", ")),
        call. = FALSE)
    }
    resp_case <- lev[hit]
    df$.resp_dummy <- as.integer(as.character(df[[condition_col]]) == resp_case)
  }
  list(Y = Y, df = df, image_col = image_col, resp_term = resp_term)
}

## Observed single-frame (log1p CP10k) decomposition from the fitted block
## proportions; the manuscript's headline per-gene frame. Internal.
.pace_single_frame <- function(Y, df, celltype_col, dec, has_condition) {
  ## Condition cohorts keep the responder spatial block separate (5-block);
  ## single-condition cohorts have an identically zero responder block, so they
  ## use the 4-block table in which it is merged back into the spatial block.
  block <- if (has_condition) dec$gene_focal_5block else dec$gene_focal_4block
  if (is.null(block)) return(NULL)
  nCount   <- as.numeric(Matrix::rowSums(Y))
  celltype <- df[[celltype_col]]
  ## The fallback exists for fits whose blocks cannot be mapped onto the observed
  ## frame at all; it must say so rather than quietly returning a different table.
  tryCatch(single_frame_decomp_obs(Y, celltype, nCount, block),
           error = function(e) {
             warning("the observed single-frame decomposition failed (",
                     conditionMessage(e), "); reporting the fit's own block table ",
                     "instead, which has different columns.", call. = FALSE)
             block
           })
}

## ===========================================================================
## Stage 1: fit the model
## ===========================================================================

#' Fit the PACE model
#'
#' `paceModel()` builds the biological and technical neighbourhood kernels from
#' the spatial coordinates and cell-type labels and fits the hierarchical
#' negative binomial mixed model by streaming penalised quasi-likelihood (with
#' the per-cell contamination correction). It returns a [PACEFit] carrying only
#' the fitted model; the reporting layers are added by [paceShrink()],
#' [paceDecompose()], and [paceDrivers()] (or all at once by [paceFit()]).
#'
#' @param object A [SpatialExperiment::SpatialExperiment] with a counts assay and
#'   two-dimensional spatial coordinates.
#' @param celltype_col colData column with the discrete cell-type annotation.
#' @param image_col Optional colData column grouping cells into images/samples.
#'   `NULL` (default) treats the object as one section.
#' @param condition_col Optional colData column with a binary condition; enables
#'   the responder spatial block for multi-sample cohorts.
#' @param assay_name Counts assay name (default `"counts"`).
#' @param h_bio,h_tech Biological and technical kernel bandwidths in micrometres
#'   (defaults 30 and 5).
#' @param contamination `"percell_hc"` (default) or `"none"`.
#' @param dispersion `"nb1"` (default) or `"nb2"`.
#' @param kernel_per_image If `TRUE`, neighbourhoods are built within each image.
#' @param image_re Second image-level random-effect block: `"none"`,
#'   `"intercept"`, `"slopes"`, or `"condition_slopes"`.
#' @param types Cell types in the order they should index the random-effect
#'   blocks. `NULL` (default) uses every observed type in alphabetical order;
#'   pass an explicit vector to reproduce a locked fit, whose block order the
#'   caller chose.
#' @param resp_term Responder interaction term prefix for condition cohorts; if
#'   `NULL` it is derived from `condition_col`.
#' @param verbose Whether to print progress.
#' @param ... Further arguments passed to the streaming fitter (`n_iter`,
#'   `threads`, `chunk_size`, `drop_sparse_neff`, `within_image`,
#'   `edge_correct`, `data_informed_tau`, `tau_shrinkage`, and the memory and
#'   approximation settings described below).
#'
#' @section Memory and approximation settings:
#'
#' Three pass-through arguments change how much memory a fit needs, or trade
#' exactness for time. Two of them are approximations that are ON by default.
#'
#' \describe{
#'   \item{`ambient_mode`}{`"cache"` (default) or `"stream"`. The contamination
#'     model needs an ambient field per cell and gene; `"cache"` materialises
#'     that n by G product, `"stream"` recomputes each chunk's columns from the
#'     weights and the counts. The numbers are identical. On a 1.2M cell,
#'     5,001 gene panel the cached product alone is over 5 GB, so `"stream"` is
#'     what makes a full transcriptome panel fit in memory at all.}
#'
#'   \item{`alpha_warmup`}{Default `6`. The per-gene dispersion MLE is re-fitted
#'     only on the first `alpha_warmup` iterations and on the last one; in
#'     between, alpha is frozen at its warmed-up value. Alpha typically settles
#'     within five iterations while the MLE is a large share of each iteration's
#'     cost, so this is on by default. **It is an approximation**: pass `Inf` to
#'     re-fit the dispersion on every iteration. Lowering it below the default
#'     does move results -- on the breast cancer cohort `4` changes the number
#'     of calls at `lfsr < 0.05`.}
#'
#'   \item{`alpha_max_n`}{Default `Inf`, meaning the dispersion MLE sees every
#'     cell. A finite value caps it at an even, deterministic subsample; the
#'     estimator, the Brent search and its tolerance are unchanged. The
#'     dispersion is one scalar per gene and its standard error falls as
#'     `1/sqrt(n)`, so most cells add little. **Validate any cap on your own
#'     cohort before trusting it.** A cap of 50,000 reproduced two cohorts of
#'     roughly 10^5 cells exactly, and on a 1.2M cell panel -- where the same
#'     cap is a far smaller fraction of the data -- it left the number of calls
#'     unchanged while swapping the identity of 72 of them. The binding
#'     quantity is not the absolute subsample.}
#' }
#' @return A [PACEFit] with the fitted model (reporting layers empty).
#' @examples
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' \donttest{
#' fit <- paceModel(spe, celltype_col = "cellType", verbose = FALSE)
#' fit
#' }
#' @rdname paceModel
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom Matrix rowSums colSums Matrix sparseMatrix bdiag
#' @export
setMethod(
  "paceModel", "SpatialExperiment",
  function(object, celltype_col, image_col = NULL, condition_col = NULL,
           assay_name = "counts", h_bio = 30, h_tech = 5,
           contamination = c("percell_hc", "none"),
           dispersion = c("nb1", "nb2"), kernel_per_image = FALSE,
           image_re = c("none", "intercept", "slopes", "condition_slopes"),
           types = NULL, resp_term = NULL, verbose = TRUE, ...) {
    contamination <- match.arg(contamination)
    dispersion    <- match.arg(dispersion)
    image_re      <- match.arg(image_re)

    prep <- .pace_prepare(object, celltype_col, image_col, condition_col,
                          assay_name, resp_term)
    res <- pace_fit_streaming(
      prep$Y, prep$df, types = types,
      celltype_col = celltype_col, image_col = prep$image_col,
      coord_cols = c("x", "y"), h_bio = h_bio, h_tech = h_tech,
      contamination = contamination, dispersion = dispersion,
      condition_col = condition_col, kernel_per_image = kernel_per_image,
      image_re = image_re, verbose = verbose, ...)

    methods::new(
      "PACEFit",
      fit       = res$fit,
      cellTypes = as.character(res$types),
      context   = list(df = res$df, X_fixed = res$X_fixed,
                       genes = colnames(res$fit$U)),
      params = list(h_bio = h_bio, h_tech = h_tech,
                    contamination = contamination, dispersion = dispersion,
                    celltype_col = celltype_col, image_col = prep$image_col,
                    condition_col = condition_col, resp_term = prep$resp_term,
                    kernel_per_image = kernel_per_image, image_re = image_re,
                    assay_name = assay_name, edge_correct = res$edge_correct))
  })

## ===========================================================================
## Stage 2: shrink the neighbour slopes
## ===========================================================================

#' Shrink the neighbour slopes
#'
#' Stabilises the fitted per-(gene, focal, neighbour) proximity slopes with
#' multivariate adaptive shrinkage, populating [neighbourSlopes()].
#'
#' One mash model is fitted per neighbour cell type, and every one of them is
#' fitted over the same genes. mash calibrates against the genes it is given,
#' so a common gene set keeps the lfsr comparable across neighbours; a gene
#' that cannot be used for one neighbour is dropped for all of them, with a
#' message.
#'
#' @param object A [PACEFit] from [paceModel()].
#' @param ... Further arguments passed to the shrinkage step, notably
#'   `null_correlation` (default `TRUE`: estimate the correlation between focal
#'   cell types under the null and pass it to mash as `V`; `FALSE` treats them
#'   as independent), `data_driven`, and `shrink_threads`.
#'
#'   `shrink_threads` (default 1) is how many R processes shrink the neighbour
#'   slices at once. The slices are independent and carry about 90% of the work,
#'   so this is where the time goes: on the full breast cancer cohort it takes
#'   the shrinkage from 55.4 to 21.1 seconds, and the called slopes agree with
#'   the serial path to 8.3e-13 with identical calls and no sign flips.
#'
#'   It is not, however, the serial computation. With `data_driven = TRUE`
#'   the parallel path cannot reproduce the serial one: `cov_pca()` draws its
#'   starting vectors from the RNG, and a serial run consumes that stream slice
#'   by slice. Above 1 each slice is seeded with its own index, so the parallel
#'   result is reproducible run to run but differs from the serial one -- by
#'   ~1e-12 at cohort scale, but by up to ~1e-6 on small slices, where there is
#'   less data to swamp the difference. Serial is the default because it is the
#'   stream the package's fixtures were made with, and because the parallel path
#'   raised a BiocParallel reducer error on one supported R version that could
#'   not be reproduced on any other. Raise it when a large cohort's shrinkage is
#'   the bottleneck. With `data_driven = FALSE` nothing draws from the stream and
#'   the two are identical.
#' @return The `PACEFit` with the shrunken neighbour slopes added.
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' fit <- paceShrink(fit)
#' head(neighbourSlopes(fit))
#' @rdname paceShrink
#' @export
setMethod("paceShrink", "PACEFit", function(object, ...) {
  slopes <- pace_shrink(object@fit, object@cellTypes,
                        resp_term = object@params$resp_term, ...)
  neighbourSlopes(object) <- slopes
  object
})

## ===========================================================================
## Stage 3: variance decomposition
## ===========================================================================

#' Variance decomposition
#'
#' Partitions per-gene expression variance into cell-type identity, spatial cell
#' state, contamination, and residual (plus a responder block for condition
#' cohorts), populating [varianceDecomposition()]. Needs the fitted object plus
#' the same `SpatialExperiment` used for [paceModel()] (to read the counts).
#'
#' The decomposition needs per-cell-type statistics of the fitted means, which
#' the fit stores; it never builds a `cells x genes` matrix. An older fit without
#' those statistics has them computed from its stored fitted means, or, if those
#' were dropped, rebuilt exactly from the fit and `spe` a block of genes at a time.
#'
#' @param object A [PACEFit] from [paceModel()].
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted.
#' @param ... Unused.
#' @return The `PACEFit` with the variance decomposition added.
#' @examples
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' fit <- paceDecompose(fit, spe)
#' head(varianceDecomposition(fit))
#' @rdname paceDecompose
#' @export
setMethod("paceDecompose", "PACEFit", function(object, spe, ...) {
  ## sparse cells x genes counts of the fitted genes, checked against the fit
  Y  <- .pace_counts_for_fit(object, spe)
  df <- object@context$df
  ## The decomposition needs only per-cell-type statistics of the fitted means:
  ## stored by the solver at fit time, taken from the matrices of an older fit
  ## that keeps them, or rebuilt a block of genes at a time for a stripped fit.
  stats <- .pace_fit_statistics(object, spe)
  dec <- pace_decompose(object@fit, df, Y, object@cellTypes,
                        object@context$X_fixed, resp_term = object@params$resp_term,
                        stats = stats)
  sf  <- .pace_single_frame(Y, df, object@params$celltype_col, dec,
                            has_condition = !is.null(object@params$condition_col))
  varianceDecomposition(object) <- list(perGene = sf, blocks = dec)
  object
})

## ===========================================================================
## Stage 4: per-pair driver tables
## ===========================================================================

#' Per-pair driver scores
#'
#' Ranks the genes mediating each focal-neighbour relationship by driver score,
#' populating [topDrivers()]. Requires [paceShrink()] and [paceDecompose()] to
#' have run first.
#'
#' The driver scores read each cell type's mean fitted mean, which the fit
#' stores. For an older fit saved without its fitted means or those statistics
#' they are rebuilt exactly from the fit and `spe`, as in [paceDecompose()], a
#' block of genes at a time; pass `spe` for such a fit.
#'
#' @param object A [PACEFit] with shrunken slopes and a decomposition.
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted. Needed
#'   only for an older fit that stores neither its fitted means nor their
#'   per-cell-type statistics.
#' @param pairs Optional list of focal-neighbour pairs to score; `NULL` scores
#'   all pairs.
#' @param ... Unused.
#' @return The `PACEFit` with the driver tables added.
#' @examples
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' fit <- paceDrivers(fit, spe)
#' names(topDrivers(fit))
#' @rdname paceDrivers
#' @export
setMethod("paceDrivers", "PACEFit", function(object, spe = NULL, pairs = NULL, ...) {
  if (nrow(object@neighbourSlopes) == 0L || length(object@varianceDecomposition) == 0L)
    stop("Run paceShrink() and paceDecompose() before paceDrivers().", call. = FALSE)
  if (!is.null(spe) && !methods::is(spe, "SpatialExperiment"))
    stop("`spe` must be the SpatialExperiment that was fitted; pass pairs by ",
         "name: paceDrivers(fit, spe, pairs = ...).", call. = FALSE)
  if (is.null(object@fit$stats) && is.null(object@fit$mu) && is.null(spe))
    stop("this fit was saved without its fitted means; pass the ",
         "SpatialExperiment it was fitted on: paceDrivers(fit, spe).", call. = FALSE)
  ## The scores need only each cell type's mean fitted mean (a statistic stored
  ## at fit time; rebuilt a block of genes at a time for an older stripped fit).
  mu_means <- .pace_fit_statistics(object, spe)$mu_mean
  args <- list(object@fit, object@neighbourSlopes,
               object@varianceDecomposition$blocks, object@cellTypes,
               mu_means = mu_means, pairs = pairs)
  ## condition cohorts score the responder interaction; pass the term and 0/1 indicator.
  if (!is.null(object@params$condition_col)) {
    args$resp_term  <- object@params$resp_term
    args$resp_dummy <- object@context$df$.resp_dummy
  }
  drivers <- tryCatch(
    do.call(pace_top_drivers, args),
    error = function(e) {
      warning("topDrivers could not be computed: ", conditionMessage(e), call. = FALSE)
      list()
    })
  topDrivers(object) <- drivers
  object
})

## ===========================================================================
## One-shot convenience: run the whole pipeline
## ===========================================================================

#' Fit a PACE model to a SpatialExperiment
#'
#' `paceFit()` runs the full pipeline in one call:
#' [paceModel()] -> [paceShrink()] -> [paceDecompose()] -> [paceDrivers()]. For
#' step-by-step control (inspecting the fitted model, re-running the downstream
#' without refitting), call those stages directly.
#'
#' @param object A [SpatialExperiment::SpatialExperiment].
#' @param pairs Optional list of focal-neighbour pairs for the driver tables.
#' @param ... Arguments passed to [paceModel()] (`celltype_col`, `image_col`,
#'   `contamination`, `dispersion`, ...).
#' @return A fully populated [PACEFit].
#' @examples
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' \donttest{
#' fit <- paceFit(spe, celltype_col = "cellType", verbose = FALSE)
#' head(neighbourSlopes(fit))
#' }
#' @rdname paceFit
#' @export
setMethod("paceFit", "SpatialExperiment", function(object, ..., pairs = NULL) {
  ## `verbose` reaches paceModel() through `...`, but the later stages are
  ## called positionally, so it has to be carried across by hand. Without this
  ## the shrinkage step reported its progress even under verbose = FALSE.
  dots <- list(...)
  verbose <- if (is.null(dots$verbose)) TRUE else isTRUE(dots$verbose)
  fit <- paceModel(object, ...)
  fit <- paceShrink(fit, verbose = verbose)
  fit <- paceDecompose(fit, object)
  paceDrivers(fit, object, pairs = pairs)
})

## ===========================================================================
## Fit-construction primitives (inspect the neighbourhood / contamination field)
## ===========================================================================

#' Build the neighbourhood kernels
#'
#' Returns the per-cell kernel-weighted neighbour abundances: the Gaussian
#' biological kernel `K_bio` (used for the proximity coefficients) and the
#' short-range exponential technical kernel `K_tech`. Useful for inspecting the
#' neighbourhood before or independently of a fit.
#'
#' @param object A [SpatialExperiment::SpatialExperiment].
#' @param celltype_col colData column with the cell-type annotation.
#' @param h_bio,h_tech Bandwidths in micrometres (defaults 30, 5).
#' @param eps Truncation radius; defaults to `3 * h_bio`.
#' @param ... Passed to the kernel builder.
#' @return A list with `K_bio` and `K_tech` (cells x cell types).
#' @examples
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' nb <- buildNeighbourhood(spe, celltype_col = "cellType")
#' dim(nb$K_bio)
#' @rdname buildNeighbourhood
#' @export
setMethod("buildNeighbourhood", "SpatialExperiment",
  function(object, celltype_col, h_bio = 30, h_tech = 5, eps = NULL, ...) {
    ct     <- as.character(SummarizedExperiment::colData(object)[[celltype_col]])
    types  <- sort(unique(ct))
    coords <- SpatialExperiment::spatialCoords(object)
    if (is.null(eps)) eps <- 3 * h_bio
    pace_neighbour_kernel(coords, ct, types, h_bio = h_bio, h_tech = h_tech,
                          eps = eps, ...)
  })

#' Build the ambient contamination field
#'
#' Returns the sparse cross-cell-type ambient weight matrix that defines the
#' per-cell contamination term (the short-range technical field a cell receives
#' from its heterotypic neighbours).
#'
#' @param object A [SpatialExperiment::SpatialExperiment].
#' @param celltype_col colData column with the cell-type annotation.
#' @param image_col Optional colData column grouping cells into images.
#' @param h_tech Technical bandwidth in micrometres (default 5).
#' @param assay_name Counts assay name (default `"counts"`).
#' @param ... Passed to the ambient-field builder.
#' @return The ambient-field object (sparse weight matrix and image index).
#' @examples
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' af <- ambientField(spe, celltype_col = "cellType")
#' names(af)
#' @rdname ambientField
#' @export
setMethod("ambientField", "SpatialExperiment",
  function(object, celltype_col, image_col = NULL, h_tech = 5,
           assay_name = "counts", ...) {
    ct     <- as.character(SummarizedExperiment::colData(object)[[celltype_col]])
    types  <- sort(unique(ct))
    coords <- SpatialExperiment::spatialCoords(object)
    ## cells x genes; kept sparse (only a few genes are densified for validation)
    Y      <- .pace_as_dgc(Matrix::t(SummarizedExperiment::assay(object, assay_name)))
    image  <- factor(if (!is.null(image_col))
                       as.character(SummarizedExperiment::colData(object)[[image_col]])
                     else rep("all", nrow(Y)))
    pace_ambient_field(coords, Y, ct, image, types, h_tech = h_tech, ...)
  })
