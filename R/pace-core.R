## pace_core.R -- cohort-agnostic PACE-MV method core (Stage 0 package extraction).
##
## A FAITHFUL lift of streaming/builders/build_bc_streaming.R: the computation and
## its order are identical, but every `Sys.getenv("R_*")` flag becomes an explicit
## function argument, and there is NO setwd / saveRDS / quit in the core. The
## solver itself (fit_pace_mvpql_streaming) and the reporting helpers
## (mvpql_to_results_multi, apply_mashr_shrinkage, mvpql_variance_decomposition_multi,
## build_random_design_multi, .compute_data_informed_weights) are assumed already
## sourced -- this file only repackages the BUILDER layer.
##
## Public surface (the future package API):
##   pace_neighbour_kernel()  -- Gaussian K_bio + exponential K_tech neighbour fields
##   pace_ambient_field()     -- sparse E^tech weight matrix W (streamed contamination)
##   pace_anchors()           -- homotypic-core negative-control anchors
##   pace_fit_streaming()     -- compose the above + streaming PQL fit  (THE method)
##   pace_shrink()            -- mash shrinkage of the neighbour slopes
##   pace_decompose()         -- per-gene variance decomposition (4-block)
##   pace_top_drivers()       -- per-pair driver-score tables
##
## Generic usage -- ANY spatial dataset, no per-cohort build script:
##   res <- pace_fit_streaming(Y, coldata,                 # counts + cell metadata
##                             celltype_col = "cellType",   # which columns hold the
##                             image_col    = "sample")     #   cell type + sample id
##   shr <- pace_shrink(res$fit, res$types)
##   dec <- pace_decompose(res$fit, res$df, res$Y, res$types, res$X_fixed,
##                         stats = res$fit$stats)
##   drv <- pace_top_drivers(res$fit, shr, dec, res$types,
##                           res$fit$stats$mu_mean)   # all focal x neighbour pairs
## Cell types, coordinates and image grouping are READ FROM THE DATA; the method
## defaults reproduce the manuscript recipe. (The Bioconductor wrapper paceFit(spe)
## reads these from a SpatialExperiment's conventions, collapsing the call to one
## argument.) Cohort-specific values (a fixed cell-type order, chosen pairs) are
## only ever needed to reproduce a *locked* fit -- see the reproduction test.

## ----------------------------------------------------------------------------
## Input helpers for the compiled neighbourhood code (src/core/).
## ----------------------------------------------------------------------------

## A single positive thread count.
.pace_thread_count <- function(threads) {
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L)
    stop("`threads` must be a single positive integer.", call. = FALSE)
  threads
}

## A radius or bandwidth: one finite number greater than 0.
.pace_check_positive <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
    stop("`", name, "` must be a single finite number greater than 0.", call. = FALSE)
  invisible(value)
}

## Coordinates as an n x 2 double matrix of finite values. PACE neighbourhoods
## are two-dimensional, so a third coordinate column is refused, not dropped.
.pace_coordinate_matrix <- function(coords) {
  coords <- as.matrix(coords)
  if (ncol(coords) != 2L)
    stop("coordinates must have exactly two columns (x, y); got ", ncol(coords), ".",
         call. = FALSE)
  storage.mode(coords) <- "double"
  if (!all(is.finite(coords)))
    stop("coordinates must be finite (no NA, NaN or Inf).", call. = FALSE)
  coords
}

## 0-based codes of `values` in `levels`; values not in `levels` get -1.
.pace_codes <- function(values, levels) {
  code <- match(values, levels) - 1L
  code[is.na(code)] <- -1L
  code
}

## Image codes for "search within each image": NA images get -1 (no image).
.pace_image_codes <- function(image) {
  image <- as.character(image)
  image_levels <- unique(image)
  image_levels <- image_levels[!is.na(image_levels)]
  list(code = .pace_codes(image, image_levels), n_images = length(image_levels))
}

## Counts as a double dgCMatrix, cells x genes (no copy when already one), with
## every stored value finite. Matrix and base matrices are converted directly;
## other matrix-like assays (DelayedArray, HDF5-backed) are realised through
## their own sparse coercion when they have one, and through as.matrix() otherwise.
.pace_as_dgc <- function(Y) {
  if (!methods::is(Y, "dgCMatrix")) {
    ## only densify when no sparse coercion exists; errors from an existing
    ## coercion (including memory errors) are raised, not hidden
    if (!methods::is(Y, "Matrix") && !is.matrix(Y))
      Y <- if (methods::canCoerce(Y, "dgCMatrix")) methods::as(Y, "dgCMatrix") else as.matrix(Y)
    if (!methods::is(Y, "dgCMatrix"))
      Y <- methods::as(methods::as(methods::as(Y, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }
  ## the same check as all(is.finite(Y@x)), without R's logical copy of it
  if (!pace_all_finite_cpp(Y@x))
    stop("counts must be finite (no NA, NaN or Inf).", call. = FALSE)
  Y
}

## The quadrature angles of the edge correction, theta_k = 2 pi k / n_angles.
.pace_edge_angles <- function(n_angles = 1000L) {
  theta <- seq(0, 2 * pi, length.out = n_angles + 1L)[-1L]
  list(cos = cos(theta), sin = sin(theta))
}

## ----------------------------------------------------------------------------
## Isotropic area-fraction edge correction.
## area_fraction(i) = area(disc(i, r) intersect image_rectangle) / (pi r^2),
## approximated by angular quadrature over `n_angles` directions:
##   r_eff(theta) = min(r, distance to the rectangle edge along theta)
##   area_fraction = (pi * sum_theta r_eff^2 / n_angles) / (pi * r^2)
## Computed in C++ (src/core/neighbourhood.cpp); cells at least r from every edge
## share one value, computed once.
## ----------------------------------------------------------------------------
pace_area_fraction <- function(coords, r,
                               xmin = NULL, xmax = NULL,
                               ymin = NULL, ymax = NULL,
                               n_angles = 1000L, threads = 1L) {
  .pace_check_positive(r, "r")
  coords <- .pace_coordinate_matrix(coords)
  if (!nrow(coords)) return(numeric(0))
  if (is.null(xmin)) xmin <- min(coords[, 1])
  if (is.null(xmax)) xmax <- max(coords[, 1])
  if (is.null(ymin)) ymin <- min(coords[, 2])
  if (is.null(ymax)) ymax <- max(coords[, 2])
  angles <- .pace_edge_angles(n_angles)
  pace_area_fraction_cpp(coords, r, xmin, xmax, ymin, ymax, angles$cos, angles$sin,
                         .pace_thread_count(threads))
}

## The original R edge-correction quadrature, vectorised over cells: used only to
## validate the compiled edge correction inside pace_ambient_field(). Same
## arithmetic as the R loop it replaced (rowSums() accumulates like sum()).
.pace_edge_fraction_reference <- function(cells, r, xmin, xmax, ymin, ymax, n_angles = 1000L) {
  theta <- seq(0, 2 * pi, length.out = n_angles + 1L)[-1L]
  cos_t <- cos(theta)
  sin_t <- sin(theta)
  x0 <- cells[, 1]
  y0 <- cells[, 2]
  ones <- rep(1, length(x0))
  d_right  <- outer(xmax - x0, ifelse(cos_t > 0, cos_t, NA), "/")
  d_left   <- outer(xmin - x0, ifelse(cos_t < 0, cos_t, NA), "/")
  d_top    <- outer(ymax - y0, ifelse(sin_t > 0, sin_t, NA), "/")
  d_bottom <- outer(ymin - y0, ifelse(sin_t < 0, sin_t, NA), "/")
  d_right[is.na(d_right)] <- Inf
  d_left[is.na(d_left)] <- Inf
  d_top[is.na(d_top)] <- Inf
  d_bottom[is.na(d_bottom)] <- Inf
  d_max <- pmin(d_right, d_left, d_top, d_bottom)
  r_eff <- pmin(r * outer(ones, rep(1, n_angles)), d_max)
  (pi * rowSums(r_eff^2) / n_angles) / (pi * r^2)
}

## ----------------------------------------------------------------------------
## Neighbour kernels: for every cell, the kernel-weighted abundance of each
## neighbour cell type within radius `eps`.
##   K_bio[i, c]  = sum_{j in type c, j != i} exp(-(d_ij / h_bio)^2)   (Gaussian)
##   K_tech[i, c] = sum_{j in type c, j != i} exp(-d_ij / h_tech)      (exponential)
## Neighbours are found and summed in C++ on a grid, without storing pairs, with
## the same distance arithmetic and summation order as the frNN + sparseMatrix
## implementation it replaces. `per_image = TRUE`: images are separate samples,
## so no cross-image neighbours; otherwise one physical section, global search.
## ----------------------------------------------------------------------------
pace_neighbour_kernel <- function(coords, celltype, types, h_bio, h_tech, eps,
                                  image = NULL, per_image = FALSE, threads = 1L) {
  .pace_check_positive(h_bio, "h_bio")
  .pace_check_positive(h_tech, "h_tech")
  .pace_check_positive(eps, "eps")
  coords <- .pace_coordinate_matrix(coords)
  n <- nrow(coords)
  neighbour_type <- .pace_codes(as.character(celltype), types)
  group <- rep(-1L, n)
  if (per_image) {
    stopifnot(!is.null(image))
    group <- .pace_image_codes(image)$code
  }
  kernels <- pace_neighbour_kernels_cpp(coords, neighbour_type, length(types), group,
                                        isTRUE(per_image), h_bio, h_tech, eps,
                                        .pace_thread_count(threads))
  K_tech <- kernels$K_tech
  K_bio  <- kernels$K_bio
  colnames(K_tech) <- paste0(types, "_near")
  colnames(K_bio)  <- types
  list(K_bio = K_bio, K_tech = K_tech)
}

## ----------------------------------------------------------------------------
## Drop (focal, neighbour) kernel columns with too little effective support.
## For each focal cell type, zero K_bio[focal cells, nb] when the centred
## column's effective sample size n_eff = sum(Kc^2) / max(Kc^2) < neff_min.
## ----------------------------------------------------------------------------
pace_drop_sparse_k <- function(K_bio, celltype, types, neff_min, verbose = TRUE) {
  dropped <- pace_drop_sparse_kernel_cpp(K_bio, .pace_codes(as.character(celltype), types),
                                         neff_min)
  if (verbose)
    message(sprintf("    drop_sparse_k: zeroed %d (focal, neighbour) pairs with n_eff < %g",
                    as.integer(dropped$n_dropped), neff_min))
  out <- dropped$kernel
  dimnames(out) <- dimnames(K_bio)
  out
}

## ----------------------------------------------------------------------------
## Within-(image, celltype) centring of each neighbour-kernel column, so the
## random slopes are estimated on the within-group deviation only.
## ----------------------------------------------------------------------------
pace_center_within_image <- function(K_bio, celltype, image, types) {
  image_levels <- unique(as.character(image))
  out <- pace_centre_within_groups_cpp(K_bio, .pace_codes(as.character(image), image_levels),
                                       .pace_codes(as.character(celltype), types),
                                       length(image_levels), length(types))
  dimnames(out) <- dimnames(K_bio)
  out
}

## ----------------------------------------------------------------------------
## Homotypic-core negative-control anchors (de-circularised contamination refs).
## A cell type's clean identity profile is estimated from its spatially-isolated
## (>= homo_frac same-type-neighbour) cells; a gene anchors type X when another
## type clearly owns it (owner_mean > owner_thresh) and X's core level is < a
## small fraction (core_thresh) of the owner's.
## Returns the memory-light form: an (n_types x G) 0/1 mask + a length-n celltype
## index (1..n_types), expanded per cell inside the solver.
## ----------------------------------------------------------------------------
## The per-type means and the same-type neighbour fractions are computed in C++
## from the sparse counts (exactly as colMeans() and mean() compute them); the
## anchor decision below stays in R.
pace_anchors <- function(coords, Y, celltype, image, types,
                         homo_frac = 0.5, owner_thresh = 0.1, core_thresh = 0.1,
                         verbose = TRUE, threads = 1L) {
  n <- nrow(Y)
  threads <- .pace_thread_count(threads)
  labels <- as.character(celltype)
  if (anyNA(labels))
    stop("cell-type labels must not be missing.", call. = FALSE)
  counts <- .pace_as_dgc(Y)
  type_code <- .pace_codes(labels, types)
  type_means <- pace_group_column_means_cpp(counts, type_code, length(types),
                                            detection = FALSE, n_threads = threads)
  dimnames(type_means) <- list(types, colnames(Y))

  ## same-type-neighbour fraction within 30 um, per image of at least 50 cells
  images <- .pace_image_codes(image)
  same_frac <- pace_same_type_fraction_cpp(.pace_coordinate_matrix(coords),
                                           .pace_codes(labels, unique(labels)),
                                           images$code,
                                           radius = 30, min_image_cells = 50L,
                                           n_threads = threads)
  core <- which(same_frac >= homo_frac)

  ## clean profile from core cells (raw type-mean fallback when a type is sparse)
  core_code <- rep(-1L, n)
  core_code[core] <- type_code[core]
  core_group_means <- pace_group_column_means_cpp(counts, core_code, length(types),
                                                  detection = FALSE, n_threads = threads)
  core_sizes <- tabulate(core_code[core_code >= 0L] + 1L, nbins = length(types))
  core_means <- type_means
  for (ti in seq_along(types)) {
    if (core_sizes[ti] >= 20) core_means[ti, ] <- core_group_means[ti, ]
  }
  anchored <- pace_anchor_mask_cpp(core_means, owner_thresh, core_thresh)
  mask <- anchored$mask
  dimnames(mask) <- list(types, colnames(Y))
  n_anchor <- anchored$n_anchor
  if (verbose)
    message(sprintf("    anchors: homotypic-core cells %d/%d; anchors/type median=%d range=[%d,%d]",
                    length(core), n, as.integer(stats::median(n_anchor)),
                    min(n_anchor), max(n_anchor)))
  list(mask = mask, idx = as.integer(celltype))
}

## ----------------------------------------------------------------------------
## Sparse E^tech ambient weight matrix W (n x n).
##   W[i, j] = exp(-d_ij / h_tech) / area_fraction(i)   for cross-celltype
##             neighbours j within rad = 3 * h_tech of i (self excluded).
## Streaming identity: W %*% Y[, g] == dense E_tech[, g]. Validated on a few
## evenly spaced genes before the fit when `validate = TRUE`.
## W and the edge fractions are built in C++ without storing neighbour pairs. The
## validation has two parts: every row is compared with E_tech rebuilt from
## dbscan::frNN pairs (using the compiled edge fractions), and up to 1000 cells are
## also compared with edge fractions recomputed by the original R quadrature
## (.pace_edge_fraction_reference()), so the edge correction is checked too.
## ----------------------------------------------------------------------------
pace_ambient_field <- function(coords, Y, celltype, image, types, h_tech,
                               edge_correct = TRUE, validate = TRUE, verbose = TRUE,
                               threads = 1L) {
  .pace_check_positive(h_tech, "h_tech")
  n <- nrow(Y)
  rad <- 3 * h_tech
  if (!is.finite(rad))
    stop("`h_tech` is too large: the ambient radius 3 * h_tech must be finite.", call. = FALSE)
  coords <- .pace_coordinate_matrix(coords)
  ct_all <- as.character(celltype)
  if (anyNA(ct_all))
    stop("cell-type labels must not be missing.", call. = FALSE)

  ## images as split() groups them: factor levels, NA in no image
  image_factor <- as.factor(image)
  image_code <- as.integer(image_factor) - 1L
  image_code[is.na(image_code)] <- -1L
  angles <- .pace_edge_angles()
  field <- pace_ambient_field_cpp(coords, .pace_codes(ct_all, unique(ct_all)), image_code,
                                  nlevels(image_factor), h_tech, isTRUE(edge_correct),
                                  angles$cos, angles$sin, .pace_thread_count(threads))
  af <- field$edge_fraction
  W <- methods::new("dgCMatrix", i = field$i, p = field$p, x = field$x, Dim = c(n, n))

  ## cheap correctness gate: W %*% Y == dense E_tech on a few spot-check genes.
  ## Genes are chosen deterministically (evenly spaced) so the check does not
  ## touch the global RNG state.
  if (validate) {
    spot_g <- unique(round(seq(1, ncol(Y), length.out = min(5L, ncol(Y)))))
    Y_spot <- as.matrix(Y[, spot_g, drop = FALSE])
    a_spot <- as.matrix(W %*% Y_spot)
    by_image <- split(seq_len(n), image)
    ## unnormalised E_tech sums from frNN pairs, per cell with a heterotypic neighbour
    raw_rows <- list()
    raw_sums <- list()
    for (si in seq_along(by_image)) {
      rows <- by_image[[si]]
      if (length(rows) < 2) next
      nn <- dbscan::frNN(coords[rows, , drop = FALSE], eps = rad)
      i_loc <- rep.int(seq_along(rows), lengths(nn$id))
      j_loc <- unlist(nn$id,   use.names = FALSE)
      d_loc <- unlist(nn$dist, use.names = FALSE)
      ct_im <- ct_all[rows]
      keep  <- ct_im[j_loc] != ct_im[i_loc]
      if (!any(keep)) next
      weighted <- exp(-d_loc[keep] / h_tech) * Y_spot[rows[j_loc[keep]], , drop = FALSE]
      sums <- rowsum(weighted, rows[i_loc[keep]])
      raw_rows[[length(raw_rows) + 1L]] <- as.integer(rownames(sums))
      raw_sums[[length(raw_sums) + 1L]] <- sums
    }
    target <- as.integer(unlist(raw_rows, use.names = FALSE))
    ## every row: frNN-based E_tech with the compiled edge fractions; cells without
    ## a heterotypic neighbour (all cells, when there is none) must have empty rows
    E_raw <- matrix(0, n, length(spot_g))
    E_all <- matrix(0, n, length(spot_g))
    if (length(target)) {
      E_raw[target, ] <- do.call(rbind, raw_sums)
      E_all[target, ] <- E_raw[target, , drop = FALSE] / af[target]
    }
    max_id <- if (n) max(abs(a_spot - E_all)) else 0
    if (length(target)) {
      ## up to 1000 evenly spaced cells, edge fractions recomputed independently
      pick <- unique(round(seq(1, length(target), length.out = min(1000L, length(target)))))
      checked <- target[pick]
      af_reference <- rep(1, length(checked))
      if (edge_correct) {
        checked_image <- as.character(image)[checked]
        for (im in unique(checked_image)) {
          in_image <- which(checked_image == im)
          rows <- by_image[[im]]
          af_reference[in_image] <- .pace_edge_fraction_reference(
            coords[checked[in_image], , drop = FALSE], rad,
            min(coords[rows, 1]), max(coords[rows, 1]),
            min(coords[rows, 2]), max(coords[rows, 2]))
        }
      }
      E_checked <- E_raw[checked, , drop = FALSE] / af_reference
      max_id <- max(max_id, abs(a_spot[checked, , drop = FALSE] - E_checked))
    }
    if (verbose)
      message(sprintf("    [W identity check] max|W%%*%%Y - dense E_tech| over %d genes = %.3e",
                      length(spot_g), max_id))
    if (max_id > 1e-8)
      stop("pace_ambient_field: W identity check FAILED (max abs diff ", max_id, ")")
  }

  list(W = W, image_idx = as.integer(image), n_images = nlevels(image))
}

## ----------------------------------------------------------------------------
## Fitted means and the contamination log-offset, rebuilt rather than stored.
## Used only for older, stripped fits that carry neither the n x G matrices nor
## the per-cell-type statistics (.pace_statistics_by_rebuild(), pace-stats.R).
##
## mu is a deterministic function of what the fit already holds: on the log
## scale eta = X B + Z U has rank at most p + q, so B, U, X_fixed and re_meta$Z
## are the factored form of it, and the dense n x G matrix is the expanded copy.
## Only the ambient field has to be recomputed, from the same counts and the
## same settings the fit used. Reconstruction is exact, not approximate.
##
## technical_offset_mat = log1p(mu_spill / mu_bio) must be rebuilt alongside mu:
## it is what gates the spillover block of the decomposition (`use_bleed`), so a
## fit missing it does not fail, it silently reports no spillover.
##
## The settings must come from the fit, never from the defaults of the exported
## ambientField(): a fit with another h_tech, image grouping, cell-type order or
## edge correction would otherwise be handed a different field and would return
## plausible, wrong numbers.
## ----------------------------------------------------------------------------
## What rebuilding mu needs beyond the fit: the counts (cells x fitted genes,
## kept sparse if the assay is) and the ambient field, or NULL for a fit
## without contamination.
.pace_mu_inputs <- function(object, spe) {
  f  <- object@fit
  df <- object@context$df
  p  <- object@params
  counts <- SummarizedExperiment::assay(spe, p$assay_name)
  if (ncol(counts) != nrow(df))
    stop("`spe` has ", ncol(counts), " cells but the fit has ", nrow(df),
         "; pass the same object used for paceModel().", call. = FALSE)
  Y <- .pace_as_dgc(Matrix::t(counts[object@context$genes, , drop = FALSE]))

  if (!identical(p$contamination, "percell_hc") || is.null(f$percell_bleed_rho))
    return(list(Y = Y, amb = NULL))

  if (is.null(p$edge_correct))
    stop("this fit predates the recording of `edge_correct` and its ambient ",
         "field cannot be rebuilt faithfully; refit, or decompose a fit that ",
         "retains `mu`.", call. = FALSE)

  amb <- pace_ambient_field(SpatialExperiment::spatialCoords(spe), Y,
                            df$celltype, df$imageID, object@cellTypes,
                            h_tech = p$h_tech, edge_correct = p$edge_correct,
                            validate = FALSE, verbose = FALSE)
  list(Y = Y, amb = amb)
}

## The two parts of the fitted mean for a block of genes (all cells x gene_idx):
## mu_bio, and mu_spill (NULL without contamination). The one place the rebuild
## formula lives. `gene_idx` selects genes (NULL: all of them); each gene column
## is computed independently, so any split of the genes gives the same values.
## These expressions must stay in step with the solver's own final pass
## (pace::final_pass_statistics), which a test ties them to.
.pace_mu_block <- function(object, inputs, gene_idx = NULL, threads = 1L) {
  f  <- object@fit
  df <- object@context$df
  if (is.null(gene_idx)) gene_idx <- seq_len(ncol(f$B))

  ## Count fits only. On the identity link the fitted mean is eta + offset,
  ## not its exponential, and a Gaussian fit has no df$nCount at all -- log(NULL)
  ## is numeric(0), so mu would silently collapse to a zero-length matrix and the
  ## decomposition built on it would fail somewhere far from here, if it failed
  ## at all. Refuse instead. The identity-link readouts are not written yet: the
  ## decomposition is a link-scale variance partition and whether that carries
  ## over is a real question, not a two-line change.
  if (identical(object@params$family, "gaussian")) {
    stop("This is a Gaussian (identity-link) fit from paceModelGaussian(). The mu ",
         "reconstruction, and the decomposition and driver readouts built on it, ",
         "are implemented for count fits only.", call. = FALSE)
  }

  ## mu_bio = exp(eta + offset), offset = log library size (pace_fit_streaming).
  eta    <- as.matrix(object@context$X_fixed %*% f$B[, gene_idx, drop = FALSE]) +
            as.matrix(f$re_meta$Z %*% f$U[, gene_idx, drop = FALSE])
  mu_bio <- pmax(exp(eta + log(df$nCount)), 1e-6)
  if (is.null(inputs$amb))
    return(list(mu_bio = mu_bio, mu_spill = NULL))

  ## Same expressions and floors as the solver's final pass -- and the same
  ## PRODUCT. The final pass gets its ambient field from pace::sparse_product_csc
  ## (see pace_sparse_product_cpp), so taking this one from R's `%*%` would put a
  ## second implementation back on the one path a test ties to the first, at a
  ## tolerance of 1e-10. Subsetting the genes first rather than passing a column
  ## range keeps this correct for a non-contiguous `gene_idx`, which this
  ## function's contract allows even though its only caller blocks contiguously.
  ambient  <- as.matrix(pace_sparse_product_cpp(
      .pace_as_dgc(inputs$amb$W), .pace_as_dgc(inputs$Y[, gene_idx, drop = FALSE]),
      1L, length(gene_idx), .pace_thread_count(threads)))
  mu_spill <- pmax(ambient * f$percell_bleed_rho, 0)
  list(mu_bio = mu_bio, mu_spill = mu_spill)
}

## ----------------------------------------------------------------------------
## THE method: compose kernels + ambient + anchors + data-informed tau and run
## the streaming PQL fit. All knobs that were `R_*` env flags in the builder are
## explicit arguments here. Returns the fit plus the working frame / X_fixed / Y
## needed by the reporting functions below.
## ----------------------------------------------------------------------------
pace_fit_streaming <- function(Y, df, types = NULL,
                               celltype_col, image_col, coord_cols = c("x", "y"),
                               h_bio = 30, h_tech = 5, eps = NULL,
                               contamination = c("percell_hc", "none"),
                               dispersion = c("nb1", "nb2"),              ## NB1 is canonical; matches paceModel()
                               condition_col = NULL,                      ## disease/condition column (optional)
                               kernel_per_image = FALSE,                  ## TRUE when images = separate samples
                               image_re = c("none", "intercept",          ## second RE block over images
                                            "slopes", "condition_slopes"),
                               drop_sparse_neff = 30,
                               within_image = TRUE,
                               edge_correct = TRUE,
                               data_informed_tau = TRUE,
                               det_min = 0.05, homo_frac = 0.5,
                               n_iter = 32L, threads = 4L, chunk_size = 128L,
                               tau_shrinkage = "adaptive",
                               alpha_warmup = 6, early_stop_tol = 2e-2, min_iter = 12L,
                               ## speed option, see fit_pace_mvpql_streaming()
                               ## The dispersion MLE runs a Brent search per gene over the
                               ## cells of that gene. alpha is ONE SCALAR per gene, so a
                               ## million cells is far past the point of diminishing returns,
                               ## and the search is the single largest term in the fit (45% of
                               ## wall on breast cancer). Subsampling evenly caps that work at
                               ## a fixed cost per gene instead of one proportional to n.
                               ## Inf keeps every cell, which is the historical behaviour.
                               alpha_max_n = Inf,
                               alpha_zero_collapse = TRUE,
                               alpha_fast_density = TRUE,
                               ## DEPRECATED: the fit keeps the per-cell-type statistics the
                               ## readouts need, and mu is rebuilt on demand from the fit and
                               ## the counts (.pace_mu_block()).
                               return_mu = FALSE,
                               ## Upper bound on the variance components, passed to the fitter.
                               ## A binding cap means the term is identified only by the ridge.
                               tau_max = 100,
                               ## "stream" avoids materialising the n x G ambient product; see
                               ## fit_pace_mvpql_streaming(). Identical numbers, less memory.
                               ambient_mode = c("cache", "stream"),
                               verbose = TRUE) {
  ambient_mode <- match.arg(ambient_mode)
  contamination <- match.arg(contamination)
  dispersion    <- match.arg(dispersion)
  image_re      <- match.arg(image_re)
  if (is.null(eps)) eps <- 3 * h_bio
  .pace_check_positive(h_bio, "h_bio")
  .pace_check_positive(h_tech, "h_tech")
  .pace_check_positive(eps, "eps")
  if (!is.finite(3 * h_tech))
    stop("`h_tech` is too large: the ambient radius 3 * h_tech must be finite.", call. = FALSE)
  use_etech <- contamination == "percell_hc"   ## E^tech ambient drives spillover
  has_cond  <- !is.null(condition_col)
  if (image_re == "condition_slopes" && !has_cond)
    stop("image_re = \"condition_slopes\" needs a `condition_col`.", call. = FALSE)

  ## ---- 1. working frame: raw labels -> factors, library size ----
  df <- as.data.frame(df)
  celltype_raw <- as.character(df[[celltype_col]])
  ## default: every observed cell type (a fixed order is only needed to reproduce
  ## a locked fit, in which case the caller passes `types` explicitly).
  if (is.null(types)) types <- sort(unique(celltype_raw))
  df$celltype <- factor(celltype_raw, levels = types)
  if (anyNA(df$celltype))
    stop(sum(is.na(df$celltype)), " cell(s) have a missing cell type or one not in `types`; ",
         "remove them or add their type to `types`.", call. = FALSE)
  df$imageID  <- factor(as.character(df[[image_col]]))
  ## The counts stay sparse from here on: every consumer (the kernels, the
  ## anchors, the data-informed weights, the solver and the decomposition) reads
  ## them column by column in C++.
  Y <- .pace_as_dgc(Y)
  ## optional max-per-celltype detection re-filter (no-op at det_min = 0.05)
  if (det_min > 0.05 + 1e-6) {
    codes <- .pace_codes(as.character(df$celltype), types)
    detection <- pace_group_column_means_cpp(Y, codes, length(types), detection = TRUE,
                                             n_threads = .pace_thread_count(threads))
    detection[tabulate(codes + 1L, nbins = length(types)) == 0L, ] <- 0
    Y <- Y[, which(apply(detection, 2, max) >= det_min), drop = FALSE]
  }
  df$nCount <- as.numeric(Matrix::rowSums(Y))
  Y_sparse <- Y
  coords <- as.matrix(df[, coord_cols])
  if (verbose)
    message(sprintf("pace_fit_streaming: %d cells x %d genes; %d images",
                    nrow(Y), ncol(Y), nlevels(df$imageID)))

  ## ---- 2. neighbour kernels (+ sparse-pair drop + within-image centring) ----
  ker <- pace_neighbour_kernel(coords, df$celltype, types, h_bio, h_tech, eps,
                               image = df$imageID, per_image = kernel_per_image,
                               threads = threads)
  K_bio  <- ker$K_bio
  K_tech <- ker$K_tech
  if (drop_sparse_neff > 0)
    K_bio <- pace_drop_sparse_k(K_bio, df$celltype, types, drop_sparse_neff, verbose)
  if (within_image)
    K_bio <- pace_center_within_image(K_bio, df$celltype, df$imageID, types)
  for (tc in types) df[[tc]] <- K_bio[, tc]   ## re_specs formula reads these by name

  ## ---- 3. fixed effects + RE spec ----
  ## E^tech path: X_fixed is intercept-only, or + a condition main effect when a
  ## disease contrast is supplied. The celltype RE block carries the neighbour
  ## slopes, crossed with the condition when present.
  offset_vec <- log(df$nCount)
  X_fixed <- if (use_etech && !has_cond) {
    matrix(1, nrow(df), 1, dimnames = list(NULL, "(Intercept)"))
  } else if (use_etech && has_cond) {
    stats::model.matrix(stats::as.formula(paste0("~ 1 + ", condition_col)), data = df)
  } else {
    cbind(`(Intercept)` = 1, K_tech)
  }
  types_rhs <- paste(types, collapse = " + ")
  celltype_formula <- if (has_cond)
    stats::as.formula(paste0("~ 1 + ", condition_col, " * (", types_rhs, ")"))
  else
    stats::as.formula(paste0("~ 1 + ", types_rhs))
  re_specs <- list(list(group_col = "celltype", formula = celltype_formula))

  ## optional SECOND RE block over images (e.g. patient-level neighbour slopes).
  ## Uses standardised kernel columns <type>_imgz; condition_slopes adds the
  ## condition interaction (the PAT_DZ design).
  if (image_re != "none") {
    ## A kernel column that is constant over all cells (e.g. zeroed for every
    ## focal by the sparse-pair drop) has zero SD: scale() returns NaN for every
    ## cell and model.matrix() then drops all rows. It carries no slope to
    ## estimate, so it is left out of the image-slope terms.
    ## scale() centres by the column mean and divides by the root mean square of
    ## the centred column; a column with no spread carries no slope to estimate.
    standardised <- lapply(types, function(tc) pace_standardise_cpp(df[[tc]]))
    names(standardised) <- types
    varying <- vapply(standardised, function(z) is.finite(z$sd) && z$sd > 0, logical(1))
    if (image_re != "intercept" && any(!varying) && verbose)
      message("    image_re: no image slope for constant kernel column(s) ",
              paste(types[!varying], collapse = ", "))
    ## sprintf(), not paste0(): paste0() turns an empty input into "_imgz"
    imgz <- sprintf("%s_imgz", types[varying])
    for (tc in types[varying]) df[[paste0(tc, "_imgz")]] <- standardised[[tc]]$values
    img_rhs <- switch(image_re,
      intercept        = "1",
      slopes           = paste(c("1", imgz), collapse = " + "),
      condition_slopes = paste(c("1", imgz, sprintf("%s:%s", condition_col, imgz)),
                               collapse = " + "))
    re_specs <- c(re_specs, list(list(
      group_col = "imageID",
      formula   = stats::as.formula(paste0("~ ", img_rhs)))))
  }

  ## ---- 4. anchors + sparse ambient field ----
  anchors <- if (contamination == "percell_hc")
    pace_anchors(coords, Y_sparse, df$celltype, df$imageID, types, homo_frac,
                 verbose = verbose, threads = threads)
  else NULL
  ambient_W <- NULL
  ambient_image_idx <- NULL
  ambient_n_images <- 0L
  if (use_etech) {
    amb <- pace_ambient_field(coords, Y_sparse, df$celltype, df$imageID, types, h_tech,
                              edge_correct = edge_correct, verbose = verbose,
                              threads = threads)
    ambient_W <- amb$W
    ambient_image_idx <- amb$image_idx
    ambient_n_images <- amb$n_images
  }

  ## ---- 5. data-informed tau weights ----
  ## The random-effect design is built ONCE here and handed to the solver, which
  ## used to rebuild the same object (Z is the largest thing in it).
  re_design <- build_random_design_multi(df, re_specs)
  data_informed_W <- NULL
  if (data_informed_tau) {
    data_informed_W <- .compute_data_informed_weights(
      re = re_design, Y = Y_sparse, df = df,
      focals = re_design$blocks[[1]]$group_levels, TYPES = types,
      celltype_col = "celltype", verbose = verbose, threads = threads)
  }

  ## ---- 6. PQL fit ----
  ## The streaming solver only implements the per-cell contamination (additive)
  ## path: it streams the ambient field as W %*% Y. With contamination = "none"
  ## there is no ambient field to stream, so we fall back to the dense oracle
  ## solver (the same solver the streaming path is byte-identical against;
  ## feasible for targeted panels). The percell_hc path is unchanged.
  if (contamination == "percell_hc") {
    fit <- fit_pace_mvpql_streaming(
      Y = Y_sparse, X_fixed = X_fixed, df = df, re_specs = re_specs, re = re_design,
      offset_vec = offset_vec, data_informed_W = data_informed_W,
      ambient_W = ambient_W, ambient_image_idx = ambient_image_idx,
      ambient_n_images = ambient_n_images,
      bleed_percell = TRUE,
      percell_anchor_mask = anchors$mask, percell_anchor_idx = anchors$idx,
      n_iter = as.integer(n_iter), tol = 5e-3,
      disp_model = dispersion,
      tau_shrinkage = tau_shrinkage,
      BPPARAM = BiocParallel::SerialParam(), n_threads = as.integer(threads),
      interior_precision = 1L, chunk_size = as.integer(chunk_size),
      alpha_warmup = alpha_warmup, alpha_max_n = alpha_max_n,
      alpha_zero_collapse = alpha_zero_collapse,
      alpha_fast_density = alpha_fast_density,
      early_stop_tol = early_stop_tol,
      min_iter = as.integer(min_iter),
      tau_max = tau_max,
      return_mu = return_mu, ambient_mode = ambient_mode, verbose = verbose)
  } else {
    ## contamination == "none": no ambient field; dense oracle (bleed_percell = FALSE).
    ## The dense oracle solver reads a dense Y; it is only reachable on targeted
    ## panels with contamination = "none".
    fit <- fit_pace_mvpql_joint_multi(
      Y = as.matrix(Y_sparse), X_fixed = X_fixed, df = df, re_specs = re_specs,
      offset_vec = offset_vec, data_informed_W = data_informed_W,
      n_iter = as.integer(n_iter), tol = 5e-3,
      disp_model = dispersion, tau_shrinkage = tau_shrinkage,
      BPPARAM = BiocParallel::SerialParam(), n_threads = as.integer(threads),
      interior_precision = 1L, chunk_size = as.integer(chunk_size),
      tau_max = tau_max, verbose = verbose)
  }

  list(fit = fit, df = df, X_fixed = X_fixed, Y = Y_sparse,
       K_tech = K_tech, K_bio = K_bio, types = types,
       ## Returned so the fit can record it: it is the one ambient-field input
       ## not otherwise recoverable from the fit, and .pace_mu() needs it to
       ## rebuild the same field.
       edge_correct = edge_correct)
}

## ----------------------------------------------------------------------------
## Reporting: mash shrinkage of the neighbour slopes.
## ----------------------------------------------------------------------------
pace_shrink <- function(fit, types, resp_term = NULL, shrink_threads = 1L,
                        verbose = TRUE, ...) {
  results_mv <- mvpql_to_results_multi(fit, keep_block = "celltype")
  apply_mashr_shrinkage(results = results_mv, focals = types,
                        neighbours = types, resp_term = resp_term,
                        shrink_threads = shrink_threads, verbose = verbose, ...)
}

## ----------------------------------------------------------------------------
## Reporting: per-gene variance decomposition with the 4-block percentage view.
## `stats` are the per-cell-type statistics (.pace_fit_statistics(), pace-stats.R);
## `Y` is the sparse cells x genes counts, read only for per-type count means.
## ----------------------------------------------------------------------------
pace_decompose <- function(fit, df, Y, types, X_fixed, resp_term = NULL,
                           stats = fit$stats, threads = 1L) {
  dec <- mvpql_variance_decomposition_stats(
    fit = fit, stats = stats, df = df, Y = Y, vars = types,
    X_fixed = X_fixed, resp_term = resp_term, focal_levels = types, threads = threads)
  g5 <- dec$gene_focal_5block
  four <- pace_four_block_shares_cpp(g5$celltype_offset_sq, g5$V_state_baseline,
                                     g5$V_state_responder, g5$V_spill, g5$V_disp)
  dec$gene_focal_4block <- tibble::tibble(
    focal = g5$focal,
    gene  = g5$gene,
    `Cell type %`     = four$pct_celltype,
    `Spatial state %` = four$pct_state,
    `Spillover %`     = four$pct_spill,
    `Residual %`      = four$pct_residual,
    celltype_offset_sq = g5$celltype_offset_sq,
    V_state = four$V_state,
    V_spill = g5$V_spill,
    V_disp  = g5$V_disp,
    Total   = four$Total,
    spec       = g5$spec,
    focal_mean = g5$focal_mean)
  dec
}

## ----------------------------------------------------------------------------
## Reporting: per-pair driver-score tables.
##   driver_score = b_shrunk^2 * spec^2 * focal_mean, lfsr < 0.05, ranked desc.
## `pairs` is a list of c(focal, neighbour); resp_term = NULL gives term == nb.
## ----------------------------------------------------------------------------
## `mu_means` is the cell types x genes matrix of mean fitted means
## (the `mu_mean` statistic, .pace_fit_statistics()), so the n x G `mu` need not be held.
pace_top_drivers <- function(fit, shrunken_long, dec, types, mu_means, pairs = NULL,
                             resp_term = NULL, resp_dummy = NULL) {
  ## default: every ordered focal != neighbour pair (a chosen subset is only a
  ## reporting convenience, not part of the method).
  if (is.null(pairs)) {
    pairs <- list()
    for (fc in types) for (nc in types) if (fc != nc) pairs <- c(pairs, list(c(fc, nc)))
  }
  g5 <- dec$gene_focal_5block
  Z_re <- .pace_as_dgc(fit$re_meta$Z)
  cells_by_ct <- lapply(types, function(c)
    pace_column_nonzero_rows_cpp(Z_re, match(paste0(c, "::(Intercept)"), colnames(Z_re))))
  names(cells_by_ct) <- types
  alpha_g <- pmax(fit$alpha, 0)
  gene_names_fit <- colnames(fit$U)
  has_resp <- !is.null(resp_term)

  ## One covariance per focal, over its own cells and all its neighbour columns:
  ## the diagonal is var(N_t) for each neighbour, and with the condition
  ## indicator as the row scale it is var(R N_t). Computing them per focal rather
  ## than per pair keeps the sparse reads to one pass per focal.
  pair_variances <- lapply(types, function(fc) {
    cells_c <- cells_by_ct[[fc]]
    columns <- match(paste0(fc, "::", types), colnames(Z_re))
    present <- which(!is.na(columns))
    var_n <- stats::setNames(rep(NA_real_, length(types)), types)
    var_rn <- var_n
    if (length(present) && length(cells_c) > 1L) {
      var_n[present] <- diag(pace_subset_covariance_cpp(Z_re, cells_c, columns[present], numeric(0)))
      if (has_resp) {
        col_R <- match(paste0(fc, "::", resp_term), colnames(Z_re))
        responder <- if (is.na(col_R)) {
          if (is.null(resp_dummy)) rep(0, length(cells_c)) else resp_dummy[cells_c]
        } else pace_subset_column_cpp(Z_re, cells_c, col_R)
        var_rn[present] <- diag(pace_subset_covariance_cpp(Z_re, cells_c, columns[present], responder))
      }
    }
    list(columns = columns, var_n = var_n, var_rn = var_rn)
  })
  names(pair_variances) <- types

  ## The long slopes table split once by (focal, neighbour, term), and the
  ## decomposition once by focal. Scanning both per pair meant 210 passes over a
  ## 70k-row table on a 15-type fit, which cost hundreds of MB of churn.
  slopes_by_pair <- split(shrunken_long,
                          paste(shrunken_long$focal, shrunken_long$neighbour, shrunken_long$term,
                                sep = "\r"))
  blocks_by_focal <- split(g5[, c("gene", "spec", "focal_mean")], as.character(g5$focal))
  rows_or_none <- function(index, key, template) {
    rows <- index[[key]]
    if (is.null(rows)) template[0, , drop = FALSE] else rows
  }

  out <- list()
  for (p in pairs) {
    fc <- p[1]
    nc <- p[2]
    pk <- paste(fc, nc, sep = "_")
    target_term <- if (is.null(resp_term)) nc else paste0(resp_term, ":", nc)
    s <- rows_or_none(slopes_by_pair, paste(fc, nc, target_term, sep = "\r"), shrunken_long) |>
      dplyr::distinct(gene, .keep_all = TRUE)
    fm <- rows_or_none(blocks_by_focal, fc, g5[, c("gene", "spec", "focal_mean")]) |>
      dplyr::distinct(gene, .keep_all = TRUE)
    if (!nrow(s) || !nrow(fm)) next

    nb_index <- match(nc, types)
    col_N <- if (is.na(nb_index)) NA_integer_ else pair_variances[[fc]]$columns[nb_index]
    if (is.na(col_N)) {
      out[[pk]] <- list(scores = s[0, ], status = "dropped (n_eff)",
                        expected_false_sign = 0, false_sign_rate = NA_real_)
      next
    }
    var_N <- pair_variances[[fc]]$var_n[[nb_index]]
    mu_bar_per_gene <- mu_means[fc, ]
    names(alpha_g) <- gene_names_fit
    names(mu_bar_per_gene) <- gene_names_fit

    ## inner_join(s, fm, by = "gene"): both are distinct by gene, so the join is
    ## one-to-one and match() does it without a hash table per pair.
    hit <- match(s$gene, fm$gene)
    s <- s[!is.na(hit), , drop = FALSE]
    hit <- hit[!is.na(hit)]
    base <- s |>
      dplyr::mutate(spec = fm$spec[hit], focal_mean = fm$focal_mean[hit],
                    mu_bar = mu_bar_per_gene[gene], alpha = alpha_g[gene])
    var_RN <- 0
    u_raw_vec <- numeric(0)
    if (has_resp) {
      ## condition x spatial: baseline slope V_S (from the raw BLUP u) plus the
      ## condition-interaction slope V_RxS (from the shrunken estimate).
      var_RN <- pair_variances[[fc]]$var_rn[[nb_index]]
      row_u <- match(colnames(Z_re)[col_N], rownames(fit$U))
      u_vec <- if (is.na(row_u))
        setNames(rep(0, length(gene_names_fit)), gene_names_fit)
      else setNames(as.numeric(fit$U[row_u, ]), gene_names_fit)
      u_raw_vec <- u_vec[base$gene]
    }
    scores <- pace_driver_scores_cpp(base$estimate_shrunk, base$spec, base$focal_mean,
                                     base$mu_bar, base$alpha, u_raw_vec, var_N, var_RN, has_resp)
    ## The vectorised expressions this replaces carried the gene names of their
    ## named inputs into the derived columns, and the stored tables of older fits
    ## have them, so each column keeps the names of its first named input.
    named_by <- function(values, ...) {
      for (source in list(...)) {
        if (!is.null(names(source))) return(stats::setNames(values, names(source)))
      }
      values
    }
    base <- base |>
      dplyr::mutate(
        driver_score       = named_by(scores$MCSD, spec, focal_mean),
        driver_score_spec4 = named_by(scores$MCSD4, spec, focal_mean),
        V_resid = named_by(scores$V_resid, alpha, mu_bar),
        V_S     = if (has_resp) named_by(scores$V_S, u_raw_vec) else scores$V_S,
        V_total = named_by(scores$V_total, u_raw_vec, alpha, mu_bar),
        R2_S    = named_by(scores$R2_S, u_raw_vec, alpha, mu_bar))
    if (!has_resp) {
      res <- base |>
        dplyr::filter(lfsr < 0.05) |>
        dplyr::arrange(dplyr::desc(driver_score)) |>
        dplyr::mutate(rank = dplyr::row_number()) |>
        dplyr::rename(b_clean = estimate_shrunk) |>
        dplyr::select(rank, gene, driver_score, driver_score_spec4, b_clean, spec, focal_mean,
                      R2_S, mu_bar, alpha, V_S, V_resid, V_total, lfsr, sd_shrunk)
    } else {
      base <- base |>
        dplyr::mutate(u_raw  = u_raw_vec,
                      V_RxS  = scores$V_RxS,
                      R2_RxS = named_by(scores$R2_RxS, u_raw_vec, alpha, mu_bar))
      res <- base |>
        dplyr::filter(lfsr < 0.05) |>
        dplyr::arrange(dplyr::desc(driver_score)) |>
        dplyr::mutate(rank = dplyr::row_number()) |>
        dplyr::rename(b_clean = estimate_shrunk) |>
        dplyr::select(rank, gene, driver_score, driver_score_spec4, b_clean, u_raw, spec, focal_mean,
                      R2_S, R2_RxS, mu_bar, alpha, V_S, V_RxS, V_resid, V_total,
                      lfsr, sd_shrunk)
    }
    status <- if (nrow(res) >= 3) "significant" else "honestly null"
    ## how many of this pair's calls are expected to have the wrong sign
    fsr <- expected_false_sign(res$lfsr)
    out[[pk]] <- list(scores = res, status = status,
                      expected_false_sign = fsr$expected_false_sign,
                      false_sign_rate = fsr$false_sign_rate)
  }
  out
}
