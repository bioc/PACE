## ============================================================
## helpers/pace_mvpql.R
## Multivariate PQL for PACE NB1 GLMM at panel scale.
##
## Joint model (after PQL linearisation):
##   Z = X B + Z_random U + E
##   vec(U) ~ N(0, diag(tau) %x% I_genes)   -- shared variance components
##   E_g ~ N(0, sigma2_g * diag(1/w_g))     -- per-gene PQL weights
##
## What's "multivariate":
##   - The variance-component vector tau (length = K_terms) is estimated by
##     EM / REML across all genes simultaneously.
##   - Per-gene marginal slopes (B[:,g]) and BLUPs (U[:,g]) are produced and
##     are drop-in compatible with the existing PACE downstream pipeline
##     (mashr, variance decomposition, MCSD).
##
## Per-gene NB1 dispersion alpha_g is updated via a Pearson moment estimator,
## so genes keep their own overdispersion despite sharing tau.
##
## What is approximate vs canonical glmmTMB NB1:
##   - PQL Gaussian-on-z working response (Breslow & Lin 1995); biased only
##     for very low mean genes (mu_bar < ~1).
##   - Variance components shared across genes (tau is a single K-vector).
##   - SEs come from the per-gene sandwich (Liang-Zeger), not full ML obs.
##     information.
## ============================================================


# Build the random-effect design matrix Z_random for `(1 + vars || celltype)`
#
# For each (celltype, term) pair we get one column. Cells contribute to the
# columns of their own celltype, with value 1 for the intercept term and the
# covariate value for slope terms.
#
# @param df data frame with celltype factor and the var columns
# @param vars character vector of slope-term names (covariates)
# @return list(Z, term_of_col, group_of_col, K_terms, K_groups)
build_random_design <- function(df, vars, group_col = "celltype") {
  if (!is.factor(df[[group_col]])) df[[group_col]] <- factor(df[[group_col]])
  groups <- levels(df[[group_col]])
  K_g <- length(groups)
  K_t <- length(vars) + 1L

  X_terms <- cbind("(Intercept)" = 1, as.matrix(df[, vars, drop = FALSE]))
  cell_groups <- as.integer(df[[group_col]])  # 1..K_g
  n <- nrow(df)

  ## column index for (term t, group g) is (t-1)*K_g + g
  col_idx <- function(t, g) (t - 1L) * K_g + g

  ii <- rep(seq_len(n), times = K_t)
  jj <- as.integer(rep(seq_len(K_t) - 1L, each = n) * K_g + rep(cell_groups, times = K_t))
  xx <- as.numeric(X_terms)
  Z  <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, K_t * K_g))

  term_of_col  <- rep(c("(Intercept)", vars), each = K_g)
  group_of_col <- rep(groups, times = K_t)
  colnames(Z)  <- paste0(group_of_col, "::", term_of_col)

  list(Z = Z,
       term_of_col  = term_of_col,
       group_of_col = group_of_col,
       term_levels  = c("(Intercept)", vars),
       group_levels = groups,
       K_terms = K_t, K_groups = K_g)
}


# Build a stacked random-effect design from multiple `(formula || group_col)`
# blocks. Each block contributes K_t_b * K_g_b columns to Z and to the
# tau-vector. Blocks share the same outer EM machinery but each has its
# own (K_t_b, K_g_b) shape -- e.g. melanoma uses
#   block "celltype" with terms = (1, Responder, neighbours, Responder:neighbours)
#   block "imageID"  with terms = (1).
#
# @param df cell-level data frame
# @param re_specs list of specs, each a `list(group_col, formula)`. The
#   formula's RHS columns become the per-group "terms" (one variance
#   component each). Use `~ 1 + Responder * (Tumour + Fibroblast + ...)`
#   for the canonical melanoma celltype block and `~ 1` for an
#   intercept-only image block.
# @return list with Z (sparse, columns concatenated across blocks),
#   term_of_col, group_of_col, block_of_col (length q), and per-block
#   metadata in `blocks` (list of K_t_b, K_g_b, term_levels, group_levels,
#   col_offsets[2] for Z indexing).
build_random_design_multi <- function(df, re_specs) {
  n <- nrow(df)
  blocks <- vector("list", length(re_specs))
  Z_list <- vector("list", length(re_specs))
  X_terms_list      <- vector("list", length(re_specs))
  cell_grp_list     <- vector("list", length(re_specs))
  cells_by_grp_list <- vector("list", length(re_specs))
  q_offset <- 0L
  for (b in seq_along(re_specs)) {
    spec <- re_specs[[b]]
    gcol <- spec$group_col
    if (!is.factor(df[[gcol]])) df[[gcol]] <- factor(df[[gcol]])
    groups <- levels(df[[gcol]])
    K_g_b  <- length(groups)
    cell_groups <- as.integer(df[[gcol]])

    ## model.matrix-based random-effect block.
    X_terms <- stats::model.matrix(spec$formula, df)
    ## model.matrix() applies na.omit and silently drops rows. The index vectors
    ## below are built from n, so a short X_terms makes sparseMatrix RECYCLE x
    ## rather than error, handing every later cell another cell's covariate.
    if (nrow(X_terms) != n) {
      stop(sprintf("re_spec '%s': model.matrix dropped %d of %d rows ",
                   gcol, n - nrow(X_terms), n),
           "(NA in the formula variables). ",
           "Remove or impute those cells before fitting.")
    }
    storage.mode(X_terms) <- "double"
    K_t_b   <- ncol(X_terms)
    term_levels <- colnames(X_terms)
    ## The design is assembled in the core. R used to build three vectors of
    ## n * K_terms elements -- at 1.2M cells and twenty terms roughly 400 MB --
    ## and hand them to Matrix::sparseMatrix(), which then sorted them into
    ## column order. Walking the groups in order and the cells ascending within
    ## a group emits the rows already sorted, so the compressed-column arrays
    ## are written once, with no triplet buffers and no sort.
    design <- pace_random_design_block_cpp(X_terms, cell_groups - 1L, K_g_b)
    Z_b <- design$Z
    colnames(Z_b) <- paste0(rep(groups, times = K_t_b), "::",
                            rep(term_levels, each = K_g_b))
    ## The same grouping the core already computed on its single pass, handed
    ## back as the list of ascending cell indices per group that split() gave.
    ## A group with no cells is integer(0), as it was before.
    starts <- design$group_start
    cbg <- lapply(seq_len(K_g_b), function(g) {
      if (starts[g + 1L] == starts[g]) integer(0)
      else design$cells_by_group[(starts[g] + 1L):starts[g + 1L]] + 1L
    })
    blocks[[b]] <- list(
      group_col    = gcol,
      term_levels  = term_levels,
      group_levels = groups,
      K_terms      = K_t_b,
      K_groups     = K_g_b,
      col_offset   = q_offset,
      n_cols       = K_t_b * K_g_b,
      term_of_col  = rep(term_levels, each = K_g_b),
      group_of_col = rep(groups,      times = K_t_b)
    )
    Z_list[[b]]       <- Z_b
    X_terms_list[[b]] <- X_terms
    cell_grp_list[[b]]     <- cell_groups
    cells_by_grp_list[[b]] <- cbg
    q_offset    <- q_offset + K_t_b * K_g_b
  }
  Z <- do.call(cbind, Z_list)
  block_of_col <- unlist(lapply(seq_along(blocks), function(b)
                                  rep(b, blocks[[b]]$n_cols)))
  term_of_col  <- unlist(lapply(blocks, `[[`, "term_of_col"))
  group_of_col <- unlist(lapply(blocks, `[[`, "group_of_col"))
  ## Block name for downstream filters (e.g. mashr only wants celltype block)
  block_name_of_col <- unlist(lapply(seq_along(blocks), function(b)
    rep(blocks[[b]]$group_col, blocks[[b]]$n_cols)))
  list(Z = Z,
       term_of_col       = term_of_col,
       group_of_col      = group_of_col,
       block_of_col      = block_of_col,
       block_name_of_col = block_name_of_col,
       blocks            = blocks,
       X_terms_list      = X_terms_list,
       cell_grp_list     = cell_grp_list,
       cells_by_grp_list = cells_by_grp_list,
       q                 = ncol(Z))
}


# Estimate the prior degrees-of-freedom d_0 for an inverse-chi-squared
# prior on per-gene variances, by marginal MLE on the cross-gene
# distribution of log(s^2). Smyth 2004, eqs 6-7. With per-gene df d_g = 1
# (one BLUP per (t, c) per gene), trigamma(0.5) = pi^2/2 ~= 4.93.
#
# Returns Inf when cross-gene variance is at or below the noise floor
# (i.e., apply full shrinkage to the panel mean).
.estimate_d0 <- function(log_s2, d_g = 1, d0_max = 5) {
  log_s2 <- log_s2[is.finite(log_s2)]
  .estimate_d0_from(stats::var(log_s2), length(log_s2), d_g = d_g, d0_max = d0_max)
}

# The second half of .estimate_d0(), given the cross-gene variance of log s2 and
# how many genes entered it. The variance itself is computed in C++
# (pace::tau_eb_summaries); the root-find stays here because it is R's own
# trigamma() and uniroot().
.estimate_d0_from <- function(v_obs, n_used, d_g = 1, d0_max = 5) {
  ## pace::estimate_d0: excess = var(log s2) - trigamma(d_g / 2), then the root
  ## of trigamma(x) = excess by Brent's method, doubled and capped.
  ##
  ## When excess <= 0 the strict marginal MLE says all the variability is
  ## sampling noise (d_0 = Inf). For BLUPs from a strongly shrunk fit that
  ## happens by construction in the early iterations, so d_0 is capped at
  ## d0_max instead and adaptive shrinkage always allows some per-gene
  ## variation; otherwise the fit never escapes "everything = panel".
  stopifnot(d_g == 1)
  pace_estimate_d0_cpp(as.numeric(v_obs), as.numeric(n_used), as.numeric(d0_max))
}


# Adaptive per-gene tau via EB shrinkage toward a per-(term, celltype)
# panel mean. Returns a q x G matrix of per-(t, c, g) variance estimates.
#
# For each (t, c):
#   panel_{t,c} = mean over genes of s2_g
#   d_0_{t,c}   = marginal-MLE prior df from cross-gene var(log s2_g)
#   tau_{t,c,g} = (d_0 * panel + s2_g) / (d_0 + 1)
# Adaptive per-gene tau via half-Cauchy prior on the SD (sigma = sqrt tau).
#
# Polson & Scott (2012) Bayesian Analysis 7: 887; Gelman 2006 Bayesian
# Analysis 1: 515. Auxiliary-variable representation:
#   sigma_{tcg} ~ Half-Cauchy(0, lambda_{t,c})
#   <=>  tau_{tcg} | a_{tcg} ~ IG(1/2, 1/a_{tcg})
#        a_{tcg} | lambda_{t,c} ~ IG(1/2, 1/lambda^2_{t,c})
#
# Variance-components analogue of vash (Lu & Stephens 2016, Bioinformatics
# 32: 3428), the GLMM-extension of limma::squeezeVar with adaptive prior.
#
# Heavy Cauchy tails on sigma let rare-celltype genes (e.g. MYLK in
# Myoepithelial) escape the panel-scale shrinkage that an inverse-chi-sq
# prior would over-impose (Gelman 2006).
#
# ECM updates per (t, c) (closed-form mode-based; see Wand 2014):
#   tau_g  <- (1/a_g + (u_g^2 + V_g)/2) / 2          [mode of IG(1, beta)]
#   a_g    <- (1/tau_g + 1/lambda^2)        / 2      [mode of IG(1, beta)]
#   lambda^2 <- 2 * mean_g(1/a_g)                    [closed-form M-step]
#
# Run a few iterations per outer PQL iter; warm-start from previous-iter values.

# Empirical-Bayes hierarchical shrinkage of per-(term, group) tau toward
# a per-term global mean, with shrinkage weight w_c = n_c / (n_c + lambda).
#
# Stabilises variance-component estimates for small clusters: for n_c >>
# lambda, w_c ~ 1 (no shrinkage); for n_c << lambda, w_c ~ 0 (full shrinkage
# to global mean). Default lambda is 0.5 * median(n_c), so the median-sized
# cluster gets ~67% local weight.
#
# @param tau_mat matrix(K_t, K_g) of per-(term, group) Schall estimates
# @param n_per_group integer vector length K_g of cells per group
# @param lambda_factor multiplier on median(n_per_group) controlling shrinkage
#   strength (default 0.5; larger -> more shrinkage)
# @return matrix same shape as tau_mat, hierarchically shrunken
## ---------------------------------------------------------------------------
## Data-informed prior weighting on the random-slope variance (R_DATA_INFORMED_TAU)
##
## Constructs a q × G multiplicative weight matrix W to scale tau_g_array per
## (random-effect column × gene).  Pairs where the data has little support for
## a (focal, neighbour, gene) slope get tau shrunk toward zero (i.e., the prior
## tightens at zero, the BLUP is shrunk to ≈ 0 at the PQL fit level).
##
## W[q, g] is constructed per random-effect column:
##   - For the celltype block, each column q encodes (focal_celltype, term).
##     - intercept / Responder / Responder:Cond terms:  W = detection_rate(focal, g)
##     - K_<nb> / Responder:K_<nb>             terms:  W = detection_rate(focal, g)
##                                                       × K_variance(focal, nb)
##   - For the imageID block: W = 1 (no data-weighting at image level)
##   - For the T_CD8_memory or other split blocks: same as celltype rules.
##
## Each per-(focal, term) weight slab is normalised so its 95th percentile = 1.
## This keeps well-supported pairs at their current PACE prior, while shrinking
## low-detection or low-K-variance pairs proportionally.
##
## Args:
##   re        the build_random_design_multi(df, re_specs) output
##   Y         n × G count matrix (used for detection rates)
##   df        cell metadata data frame (used for K_<nb> columns + celltype)
##   focals    vector of focal celltype names from re$blocks[[celltype]]
##   TYPES     vector of neighbour celltype names (= focals in symmetric design)
##
## The weight is detection_rate × normalised-K-variance, row-normalised so the
## most-informative gene per slope has W = 1 (linear; no exponent), with a fixed
## 1e-8 floor for numerical safety only.
##
## Returns: q × G matrix of multiplicative weights to apply as
##   tau_g_array <- tau_g_array * W
.compute_data_informed_weights <- function(re, Y, df, focals, TYPES,
                                            celltype_col = "celltype",
                                            verbose = TRUE, threads = 1L) {
  q <- ncol(re$Z); G <- ncol(Y)
  W <- matrix(1, q, G, dimnames = list(colnames(re$Z), colnames(Y)))
  if (!(celltype_col %in% names(df))) {
    if (verbose) cat(sprintf("  [data-informed tau] '%s' not in df; returning identity weights\n",
                              celltype_col))
    return(W)
  }
  ct_chr <- as.character(df[[celltype_col]])

  ## Detection rate per (focal, gene): mean of (Y > 0) over the focal's cells,
  ## computed in C++ from the sparse counts exactly as colMeans(Y[cells, ] > 0).
  ## A focal with no cells keeps rate 0.
  det_rate <- matrix(0, length(focals), G,
                     dimnames = list(focals, colnames(Y)))
  focal_code <- .pace_codes(ct_chr, focals)
  det_means <- pace_group_column_means_cpp(.pace_as_dgc(Y), focal_code, length(focals),
                                           detection = TRUE,
                                           n_threads = .pace_thread_count(threads))
  focal_sizes <- tabulate(focal_code[focal_code >= 0L] + 1L, nbins = length(focals))
  for (f_index in seq_along(focals)) {
    if (focal_sizes[f_index] > 0L) det_rate[f_index, ] <- det_means[f_index, ]
  }

  ## K-variance per (focal, neighbour): the diagonal of the focal's kernel-column
  ## covariance, which is stats::var over that focal's cells. A focal with fewer
  ## than five cells carries no usable variance and stays NA, as before.
  K_var <- matrix(NA_real_, length(focals), length(TYPES),
                  dimnames = list(focals, TYPES))
  present <- TYPES[TYPES %in% names(df)]
  if (length(present)) {
    covariance <- .pace_group_covariances(df[, present, drop = FALSE], ct_chr, focals, threads)
    for (f_index in seq_along(focals)) {
      if (focal_sizes[f_index] < 5L) next
      K_var[f_index, present] <- diag(matrix(covariance[, , f_index], length(present),
                                             length(present)))
    }
  }
  ## Normalise K_var per focal to [0, 1] so the weighting is scale-free across focals
  K_var_norm <- pace_normalise_rows_to_max_cpp(K_var)
  dimnames(K_var_norm) <- dimnames(K_var)

  ## Walk through each random-effect block and record, for every random-effect
  ## column, which focal's detection profile it takes and how its neighbour's
  ## K-variance scales it. A column of another block (imageID) keeps weight 1.
  focal_of_row <- integer(q)
  scale_of_row <- numeric(q)
  for (bi in seq_along(re$blocks)) {
    blk_i <- re$blocks[[bi]]
    is_celltype_blk <- identical(blk_i$group_col, celltype_col) ||
                       all(blk_i$group_levels %in% focals)
    if (!is_celltype_blk) next
    for (t in seq_len(blk_i$K_terms)) {
      term <- blk_i$term_levels[t]
      ## Determine if this term is a K-slope (i.e., references a neighbour celltype)
      ## Possible patterns:  "Tumour", "ResponderPD:Tumour", "(Intercept)", "ResponderPD"
      ## A term carries a K-slope if any of TYPES is a substring of the term name.
      K_nb_match <- TYPES[vapply(TYPES, function(t_nb)
        grepl(paste0("(^|:)", t_nb, "$"), term, fixed = FALSE),
        logical(1))]
      is_K_slope <- length(K_nb_match) > 0L

      for (c in seq_len(blk_i$K_groups)) {
        focal <- blk_i$group_levels[c]
        col_idx <- blk_i$col_offset + (t - 1L) * blk_i$K_groups + c
        focal_of_row[col_idx] <- match(focal, focals)
        ## intercept / Responder terms take the detection profile unscaled
        scale_of_row[col_idx] <- if (is_K_slope) K_var_norm[focal, K_nb_match[1L]] else 1
      }
    }
  }

  ## Per-term normalisation: divide each row by its maximum so the most-
  ## informative gene in that row has W = 1 (no extra shrinkage relative to
  ## PACE's own EB-estimated tau) and everything else scales linearly with
  ## its detection x K-variance. No exponent, no arbitrary floor -- just
  ## numerical safety against W = 0 in the per-gene WLS solve
  ## (pace::data_informed_weights).
  W <- pace_data_informed_weights_cpp(det_rate, focal_of_row, scale_of_row)
  dimnames(W) <- list(colnames(re$Z), colnames(Y))

  if (verbose) {
    cat(sprintf("  [data-informed tau] q x G weight matrix: dim %d x %d, range [%.3f, %.3f], median %.3f\n",
                nrow(W), ncol(W), min(W), max(W), median(W)))
  }
  W
}

.shrink_tau_hierarchical <- function(tau_mat, n_per_group, lambda_factor = 0.5) {
  out <- pace_tau_hierarchical_cpp(tau_mat, as.integer(n_per_group), lambda_factor)
  dimnames(out) <- dimnames(tau_mat)
  out
}

# Bound the variance components from above.
#
# Every tau write in the engines is a pmax against a floor. Without a ceiling a
# column that is identified only by the ridge climbs by mean(u^2) per iteration
# and the inflated posterior variance feeds the next update, so the term never
# settles. A binding cap is diagnostic, not cosmetic: it means the term is not
# identified by the data, so it is reported rather than applied silently.
.clamp_tau <- function(tau, tau_max, it = NA_integer_) {
  if (is.null(tau_max) || !is.finite(tau_max)) return(tau)
  capped <- pace_tau_clamp_cpp(tau, tau_max)
  n_bind <- capped$n_binding
  if (n_bind > 0L) {
    warning(sprintf("iter %s: %d variance component(s) reached tau_max = %.4g. ",
                    as.character(it), n_bind, tau_max),
            "A binding cap means the term is identified only by the ridge; ",
            "inspect the design rather than raising the cap.", call. = FALSE)
    tau <- capped$tau
    dimnames(tau) <- dimnames(capped$tau)
  }
  tau
}


# Force a double interior when the ridge is small relative to the diagonal it
# is added to.
#
# The ridge (1/tau) is accumulated inside the float kernel, where the smallest
# representable increment on a diagonal of order sum(w) is eps_float * sum(w).
# Below that the regularisation is quantised away and the terms it identifies
# come loose. The threshold scales as 1/n, so the same model and the same
# design lose their ridge purely by growing the cohort. Threshold (~84 float
# ulps) follows the correctness audit's recommendation.
## min()/max() of an all-NA or empty vector warn and return +/-Inf. That is the
## answer we want, but suppressing the warning would also hide any other, so
## take the empty case explicitly instead.
.finite_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) Inf else min(x)
}
.finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) -Inf else max(x)
}

.ridge_needs_double <- function(lam_chk, w_chk) {
  .ridge_needs_double_from(.finite_min(lam_chk), .finite_max(colSums(w_chk, na.rm = TRUE)))
}

# The same decision from the two statistics it reads, so a caller that builds the
# chunk's weights a few genes at a time can accumulate them without holding the
# whole chunk: both are per gene, so the decision is the chunk's either way.
.ridge_needs_double_from <- function(lam_min, w_max) {
  is.finite(lam_min) && is.finite(w_max) && lam_min < 1e-5 * w_max
}


.adaptive_tau_half_cauchy <- function(s2_mat, K_t, K_g, gene_names = NULL,
                                       lambda_sq_prev = NULL,
                                       a_prev         = NULL,
                                       n_em_iter      = 5L,
                                       tau_floor      = 1e-4) {
  stopifnot(nrow(s2_mat) == K_t * K_g)
  fitted <- pace_tau_half_cauchy_cpp(
    s2_mat, as.integer(n_em_iter), tau_floor,
    if (is.null(lambda_sq_prev)) numeric(0) else as.numeric(lambda_sq_prev),
    if (is.null(a_prev)) matrix(0, 0, 0) else as.matrix(a_prev))
  out <- fitted$tau
  dimnames(out) <- list(rownames(s2_mat), gene_names)
  attr(out, "lambda_sq") <- fitted$lambda_sq
  attr(out, "a")         <- fitted$a
  attr(out, "panel")     <- fitted$panel
  out
}

.adaptive_tau_eb <- function(s2_mat, K_t, K_g, gene_names = NULL,
                              d0_min = 1, tau_floor = NULL,
                              p_fixed = 0L) {
  q <- nrow(s2_mat)
  stopifnot(q == K_t * K_g)
  ## Wolfinger-O'Connell (1993) / Schall (1991) REML correction on the variance-
  ## component EM update. PACE's PQL EM update used mean(u^2 + V), the ML moment
  ## estimator, which underestimates tau for sparse counts (Lin-Breslow 1996
  ## JASA 91:1007). The inflation q / (q - p_fixed) is ~1.002 at the gene-pool
  ## level but 1.20-1.40 per block, so it is applied to both the per-gene s2 and
  ## the panel mean the shrinkage targets. Off unless R_REML_TAU is set.
  use_reml <- nzchar(Sys.getenv("R_REML_TAU", unset = ""))
  reml_fac <- if (use_reml && K_g > p_fixed) pmin(K_g / pmax(K_g - p_fixed, 1L), 2) else 1
  summaries <- pace_tau_eb_summaries_cpp(s2_mat, reml_fac)
  ## Cap the minimum d0: very weak shrinkage produces wild per-gene tau that
  ## destabilise the per-gene WLS solve. d0_min = 1 means at least 50/50 between
  ## the gene and the panel.
  d0 <- vapply(seq_len(q), function(k)
    max(.estimate_d0_from(summaries$log_variance[k], summaries$n_log_finite[k]), d0_min),
    numeric(1))
  out <- pace_tau_eb_apply_cpp(s2_mat, reml_fac, summaries$panel, summaries$panel_median, d0,
                               if (is.null(tau_floor)) NA_real_ else tau_floor)
  dimnames(out) <- list(rownames(s2_mat), gene_names)
  attr(out, "panel") <- summaries$panel
  attr(out, "d0")    <- d0
  out
}

# NB1 and NB2 MLEs for the dispersion alpha given fitted mu and counts y.
#
# Reliable replacement for the Pearson moment estimator, which collapses to ~0
# under PQL when the working response makes mu ~= y by construction. Both are
# pace::dispersion_mle: Brent's minimiser over log alpha, with the negative
# binomial density in its mean parameterisation (NB1 takes size = mu / alpha,
# NB2 the constant size = 1 / alpha).
#
# `max_n` subsamples evenly for speed. NB1 dispersion is well determined from
# ~10k cells, but alpha differences propagate through the PQL weights into the
# mashr and MCSD ranking, so the default Inf is what preserves the manuscript.
#
# `zero_collapse` sums the zero-count cells in closed form (exact algebra,
# reassociated, so alpha moves in its last digits); see the design note.
.alpha_nb1_mle <- function(y, mu, max_n = Inf, zero_collapse = FALSE,
                           fast_density = FALSE) {
  pace_dispersion_mle_cpp(as.numeric(y), as.numeric(mu), nb2 = FALSE,
                          zero_collapse = zero_collapse, max_cells = max_n,
                          fast_density = fast_density)
}

.alpha_nb2_mle <- function(y, mu, max_n = Inf) {
  pace_dispersion_mle_cpp(as.numeric(y), as.numeric(mu), nb2 = TRUE,
                          zero_collapse = FALSE, max_cells = max_n,
                          fast_density = FALSE)
}


# One penalised WLS solve for one gene
#
# Solves (D' W D + Lambda) [beta; u] = D' W z where D = [X | Z_random] and
# Lambda is block-diagonal with zeros on the X block and 1/tau (broadcast
# across celltypes within each term) on the Z block.
.solve_one_gene <- function(X, Z, w, z, lam_diag) {
  ## D' W D as 2x2 blocks
  Xw <- X * w
  Zw <- Z * w  # element-wise scale of sparse rows
  XtWX <- crossprod(X, Xw)
  XtWZ <- crossprod(X, Zw)
  ZtWZ <- crossprod(Z, Zw)
  ZtWZ <- as(ZtWZ, "dgCMatrix")
  ## Add penalty on Z block diagonal
  diag(ZtWZ) <- diag(ZtWZ) + lam_diag

  XtWz <- crossprod(X, w * z)
  ZtWz <- crossprod(Z, w * z)

  ## Stack into one symmetric system
  ## A = [XtWX  XtWZ]
  ##     [t(XtWZ) ZtWZ + Lambda]
  ## b = [XtWz; ZtWz]
  p <- ncol(X); q <- ncol(Z)
  A <- Matrix::bdiag(Matrix::Matrix(XtWX), ZtWZ)
  ## Insert XtWZ blocks
  A[seq_len(p), p + seq_len(q)] <- as(XtWZ, "dgCMatrix")
  A[p + seq_len(q), seq_len(p)] <- as(t(XtWZ), "dgCMatrix")
  b <- c(as.numeric(XtWz), as.numeric(ZtWz))

  sol <- as.numeric(Matrix::solve(A, b))
  list(beta = sol[seq_len(p)], u = sol[p + seq_len(q)], A = A, b = b)
}


# Fast per-gene WLS solve exploiting Z's celltype-block-diagonal structure
#
# For build_random_design()'s Z, every cell has nonzeros only in columns
# belonging to its own celltype. So Z'W Z is block-diagonal across celltypes
# (each block is K_t x K_t). Constructing the system directly via dense
# per-celltype crossprods skips the sparse Matrix machinery (Z*w row scale,
# sparse crossprod, bdiag, slot assignment) that dominated the per-gene
# cost at panel scale -- profile showed ~40% of fit time in those Matrix
# internals.
#
# Returns the same beta/u as .solve_one_gene plus the full diag(A^-1)
# (computed via dense Cholesky once).
#
# @param X cells x p fixed-effect design (dense)
# @param X_terms cells x K_t random-effect term values (dense, includes
#   intercept column)
# @param cells_by_ct list of length K_g; cells_by_ct[[c]] is the integer
#   vector of cell row-indices belonging to celltype c (in the order Z's
#   columns are laid out by build_random_design)
# @param w cell weights (length n)
# @param z working response (length n)
# @param lam_diag length-q penalty (1/tau broadcast)
# @param K_t number of random-effect terms per group (incl. intercept)
# @param K_g number of groups (celltypes)
.solve_one_gene_block <- function(X, X_terms, cells_by_ct, w, z, lam_diag,
                                   K_t, K_g) {
  p <- ncol(X); q <- K_t * K_g
  ## Fixed-block crossprods (dense, full panel)
  Xw       <- X * w
  XtWX     <- crossprod(X, Xw)
  XtWz_vec <- as.numeric(crossprod(X, w * z))
  ## Pre-allocate dense ZtWZ (block-diagonal, K_g blocks of K_t x K_t),
  ## XtWZ (p x q), and ZtWz (length q)
  ZtWZ <- matrix(0, q, q)
  XtWZ <- matrix(0, p, q)
  ZtWz <- numeric(q)
  ## Column index in q for (term t, group c) is (t-1)*K_g + c -- matches
  ## the column order in build_random_design().
  for (c in seq_len(K_g)) {
    idx <- cells_by_ct[[c]]
    if (!length(idx)) next
    Xt_c  <- X_terms[idx, , drop = FALSE]    # n_c x K_t
    w_c   <- w[idx]
    Wxt_c <- Xt_c * w_c                       # n_c x K_t
    cols  <- ((seq_len(K_t) - 1L) * K_g + c)
    ZtWZ[cols, cols] <- crossprod(Xt_c, Wxt_c)            # K_t x K_t
    XtWZ[, cols]     <- crossprod(X[idx, , drop = FALSE], Wxt_c)  # p x K_t
    ZtWz[cols]       <- crossprod(Xt_c, w_c * z[idx])
  }
  diag(ZtWZ) <- diag(ZtWZ) + lam_diag
  ## Assemble dense A = [XtWX, XtWZ; XtWZ^T, ZtWZ + Lambda]
  A <- rbind(cbind(XtWX,  XtWZ),
             cbind(t(XtWZ), ZtWZ))
  b <- c(XtWz_vec, ZtWz)
  ## Cholesky-based solve + diag(A^-1) extraction
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(R)) {
    return(list(beta = rep(NA_real_, p), u = rep(NA_real_, q),
                Ainv_diag = rep(NA_real_, p + q)))
  }
  sol  <- backsolve(R, backsolve(R, b, transpose = TRUE))
  Ainv <- chol2inv(R)
  list(beta      = sol[seq_len(p)],
       u         = sol[p + seq_len(q)],
       Ainv_diag = diag(Ainv))
}


# Multi-block analogue of .solve_one_gene_block.
#
# For each block b, ZtWZ's within-block sub-matrix is block-diagonal across
# that block's groups (same trick as single-block). Cross-block sub-matrices
# (b1 vs b2) are NOT block-diagonal -- each cell contributes to exactly one
# (group_b1, group_b2) intersection -- but they're cheap to build by
# bucketing cells in each group of b1 by their b2 group.
.solve_one_gene_multiblock <- function(X, X_terms_list, cell_grp_list,
                                        cells_by_grp_list,
                                        w, z, lam_diag, blocks) {
  p   <- ncol(X)
  B_n <- length(blocks)
  q   <- sum(vapply(blocks, function(b) b$n_cols, integer(1)))

  Xw       <- X * w
  XtWX     <- crossprod(X, Xw)
  XtWz_vec <- as.numeric(crossprod(X, w * z))

  ZtWZ <- matrix(0, q, q)
  XtWZ <- matrix(0, p, q)
  ZtWz <- numeric(q)

  ## Within-block: block-diagonal by group
  for (b in seq_len(B_n)) {
    blk   <- blocks[[b]]
    K_t_b <- blk$K_terms; K_g_b <- blk$K_groups
    Xt_b  <- X_terms_list[[b]]
    cbg   <- cells_by_grp_list[[b]]
    for (g in seq_len(K_g_b)) {
      idx <- cbg[[g]]
      if (!length(idx)) next
      Xt_g  <- Xt_b[idx, , drop = FALSE]
      w_g   <- w[idx]
      Wxt_g <- Xt_g * w_g
      cols_in_b <- blk$col_offset + ((seq_len(K_t_b) - 1L) * K_g_b + g)
      ZtWZ[cols_in_b, cols_in_b] <- crossprod(Xt_g, Wxt_g)
      XtWZ[, cols_in_b]          <- crossprod(X[idx, , drop = FALSE], Wxt_g)
      ZtWz[cols_in_b]            <- crossprod(Xt_g, w_g * z[idx])
    }
  }

  ## Cross-block (b1 < b2): bucket cells in each group of b1 by their
  ## b2-group, then per (g1, g2) intersection do one dense crossprod.
  if (B_n >= 2L) {
    for (b1 in seq_len(B_n - 1L)) {
      blk1   <- blocks[[b1]]
      K_t_1  <- blk1$K_terms; K_g_1 <- blk1$K_groups
      Xt_1   <- X_terms_list[[b1]]
      cbg_1  <- cells_by_grp_list[[b1]]
      for (b2 in (b1 + 1L):B_n) {
        blk2   <- blocks[[b2]]
        K_t_2  <- blk2$K_terms; K_g_2 <- blk2$K_groups
        Xt_2   <- X_terms_list[[b2]]
        grp2   <- cell_grp_list[[b2]]
        for (g1 in seq_len(K_g_1)) {
          idx1 <- cbg_1[[g1]]
          if (!length(idx1)) next
          ## Bucket idx1 by their group in block b2.
          g2_of_idx1 <- grp2[idx1]
          splits     <- split(idx1, g2_of_idx1)
          for (g2_str in names(splits)) {
            g2  <- as.integer(g2_str)
            idx <- splits[[g2_str]]
            if (!length(idx)) next
            Xt_1_sub <- Xt_1[idx, , drop = FALSE]
            Xt_2_sub <- Xt_2[idx, , drop = FALSE]
            w_sub    <- w[idx]
            block_xy <- crossprod(Xt_1_sub, Xt_2_sub * w_sub)  # K_t_1 x K_t_2
            cols1 <- blk1$col_offset + ((seq_len(K_t_1) - 1L) * K_g_1 + g1)
            cols2 <- blk2$col_offset + ((seq_len(K_t_2) - 1L) * K_g_2 + g2)
            ZtWZ[cols1, cols2] <- block_xy
            ZtWZ[cols2, cols1] <- t(block_xy)
          }
        }
      }
    }
  }

  diag(ZtWZ) <- diag(ZtWZ) + lam_diag

  A <- rbind(cbind(XtWX,  XtWZ),
             cbind(t(XtWZ), ZtWZ))
  b_vec <- c(XtWz_vec, ZtWz)

  R <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(R)) {
    return(list(beta = rep(NA_real_, p), u = rep(NA_real_, q),
                Ainv_diag = rep(NA_real_, p + q)))
  }
  sol  <- backsolve(R, backsolve(R, b_vec, transpose = TRUE))
  Ainv <- chol2inv(R)
  list(beta      = sol[seq_len(p)],
       u         = sol[p + seq_len(q)],
       Ainv_diag = diag(Ainv))
}


# Vectorised multi-block solve over a chunk of genes.
#
# Trick: the per-celltype crossprods that build ZtWZ / XtWZ / ZtWz are linear
# in `w` and `z`. So per chunk of genes we can compute one big dgemm that
# produces every gene's per-celltype block at once -- amortising R-level
# loop overhead and BLAS dispatch across the chunk. Per-gene Cholesky and
# inverse stay per-gene (parallel via BPPARAM).
#
# Per-celltype crossprod trick (within-block, single gene):
#   ZtWZ_block_c[t1, t2] = sum_{i in c} X_terms[i, t1] * X_terms[i, t2] * w[i]
# Vectorised over all G_chunk genes:
#   M_pair[i, kk] = X_terms[idx_c, t1[kk]] * X_terms[idx_c, t2[kk]]
#   crossprod(M_pair, w[idx_c, gene_idx]) -> K_t² × G_chunk
#   reshape to (K_t, K_t, G_chunk).
# Same trick for XtWZ (pairs of (p, t) indices) and for cross-block sub-
# matrices (pairs of (t in b1, t in b2) indices over the (g1, g2) intersection).
#
# @param X         cells × p fixed-effect design (dense)
# @param X_terms_list, cell_grp_list, cells_by_grp_list  see build_random_design_multi
# @param w, z      full n × G PQL weights and working response
# @param lam_diag  full q × G  per-gene penalty
# @param blocks    re_meta$blocks (per-block layout meta)
# @param gene_idx  integer vector of which genes to process this chunk
# @param BPPARAM   BiocParallel back-end for the per-gene Cholesky+solve
# @param n_threads optional integer; OpenMP threads for the C++ kernel.
#   If NULL, derived from `bpworkers(BPPARAM)` (back-compat). Pass an
#   explicit value to decouple kernel threading from the bplapply fork
#   (e.g. BPPARAM = SerialParam() with n_threads = 4 for OMP-only).
# @param interior_precision 0 (default, double) or 1 (float interior).
#   See core/gene_solve.hpp for what the float interior costs and buys.
# @return list of per-gene results, one per element of gene_idx
.solve_genes_chunk_multiblock <- function(X, X_terms_list, cell_grp_list,
                                            cells_by_grp_list,
                                            w, z, lam_diag, blocks,
                                            gene_idx,
                                            BPPARAM = BiocParallel::SerialParam(),
                                            n_threads = NULL,
                                            interior_precision = 0L) {
  G_chunk <- length(gene_idx)
  p   <- ncol(X); B_n <- length(blocks)
  q   <- sum(vapply(blocks, function(b) b$n_cols, integer(1)))

  ## The solve is pace::solve_genes_chunk: the within-block tensors, the
  ## cross-block tensors over the cells two groups share, and the per-gene
  ## Cholesky, in one call per chunk. Every loop writes its own output slot, so
  ## the result does not depend on the thread count. The reference it is read
  ## against is the explicit-Z solve in tests/testthat/test-gene-solve.R.
  if (is.null(n_threads))
    n_threads <- tryCatch(max(1L, BiocParallel::bpworkers(BPPARAM)),
                          error = function(e) 1L)

  ## Take the columns only when they really are a subset. Every caller reaches
  ## here with the identity -- the chunk loops subset w, z and lam_diag before
  ## the call and then pass seq_along() of the chunk -- and subsetting by the
  ## identity still copies the whole matrix. At 1.2M cells and the default
  ## 128-gene chunk that is 1.25 GB apiece, three times over, to rebuild the
  ## matrices we were handed. The binding borrows these through a const span
  ## and the core copies into its own storage, so passing them through
  ## unsubsetted cannot be seen by the caller.
  columns_are_identity <- identical(gene_idx, seq_len(ncol(w)))
  if (!columns_are_identity) {
    w        <- w[, gene_idx, drop = FALSE]
    z        <- z[, gene_idx, drop = FALSE]
    lam_diag <- lam_diag[, gene_idx, drop = FALSE]
  }

  res <- pace_solve_genes_chunk_cpp(
    x_fixed  = X,
    w        = w,
    z        = z,
    lam_diag = lam_diag,
    q_total  = q,
    blocks   = blocks,
    terms_list = lapply(X_terms_list, function(m) {
                   storage.mode(m) <- "double"; m }),
    cells_by_group_list = cells_by_grp_list,
    cell_group_list     = cell_grp_list,
    single_precision    = as.integer(interior_precision) != 0L,
    n_threads           = .pace_thread_count(n_threads))

  ## The three matrices go back as matrices. This used to explode them into
  ## G_chunk lists of three vectors, which the caller then looped over to write
  ## column by column -- six per-gene allocations, and two full copies of data
  ## the core had already laid out contiguously, to end up with the same three
  ## matrices. The caller assigns whole blocks instead.
  res
}



# Multivariate PQL fit
#
# @param Y cells × genes count matrix (integer)
# @param X_fixed cells × p design matrix for fixed effects (must include
#   intercept column)
# @param df data frame with celltype + vars columns (used to build Z_random)
# @param vars character vector of neighbour-count covariates (random slopes)
# @param offset_vec length-cells log(nCount), or NULL
# @param n_iter PQL iterations (default 3)
# @param tol convergence on max delta-eta per gene
# @param tau_init initial variance components (length K_terms = 1 + length(vars))
# @param BPPARAM BiocParallel back-end for the per-gene WLS solve loop
#   inside each outer iter (default SerialParam). Each gene's solve is
#   independent given the current (z, w, tau), so this parallelises cleanly;
#   the cross-gene tau-EM step still runs serially on the gathered results.
# @param tau_shrinkage One of:
#   * "shared" (legacy): single tau_{t,c} per (term, celltype), shared
#     across all genes.
#   * "adaptive" (default): per-gene tau_{t,c,g} with empirical-Bayes
#     shrinkage toward the per-(t,c) panel mean. Prior strength estimated
#     by marginal MLE on the cross-gene distribution of log(BLUP^2 + V),
#     following Smyth 2004 (limma-trend) generalised to NB1-GLMM variance
#     components. Lets high-variance rare-celltype genes (e.g. MYLK in
#     Myoepithelial) escape the panel-mean shrinkage that hurt them under
#     "shared", at the cost of one extra log/uniroot per (t,c) per iter.
# @return list:
#   - B (p × g) per-gene fixed effects
#   - U (q × g) per-gene BLUPs (rows ordered as build_random_design)
#   - se_B (p × g) per-gene SE of fixed effects (sandwich)
#   - se_U (q × g) per-gene SE of BLUPs (sandwich, conditional on tau)
#   - alpha (g) per-gene NB1 dispersion (Pearson moment estimator)
#   - tau (K_terms) shared variance components
#   - mu (cells × g) fitted means
#   - re_meta from build_random_design (term_of_col, group_of_col, etc.)
#   - n_iter, converged, secs, history (per-iter diagnostics)
fit_pace_mvpql <- function(Y, X_fixed, df, vars,
                           offset_vec = NULL,
                           n_iter = 3, tol = 1e-3,
                           tau_init = NULL,
                           tau_shrinkage = c("shared", "adaptive", "half_cauchy", "hierarchical"),
                           BPPARAM = BiocParallel::SerialParam(),
                           alpha_max_n = Inf,
                           verbose = TRUE) {
  tau_shrinkage <- match.arg(tau_shrinkage)

  if (is.null(offset_vec)) offset_vec <- rep(0, nrow(Y))
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"

  n <- nrow(Y); g_n <- ncol(Y); p <- ncol(X_fixed)
  re <- build_random_design(df, vars)
  Z  <- re$Z; q <- ncol(Z)
  K_t <- re$K_terms; K_g <- re$K_groups
  if (verbose) cat(sprintf("  [mvpql] n=%d  g=%d  p_fixed=%d  q_random=%d (= %d terms x %d groups)\n",
                           n, g_n, p, q, K_t, K_g))

  ## Precompute structural helpers for the fast block-by-celltype solve.
  ## X_terms: dense n x K_t matrix of random-effect term values (intercept +
  ## slope covariates) -- the nonzero VALUES that Z holds, just in dense
  ## form for cheap per-celltype crossprods.
  ## cells_by_ct: list of cell-row indices per celltype (matches the column
  ## ordering build_random_design uses: column (t,c) at index (t-1)*K_g + c).
  X_terms_dense <- cbind("(Intercept)" = 1, as.matrix(df[, vars, drop = FALSE]))
  storage.mode(X_terms_dense) <- "double"
  if (!is.factor(df$celltype)) df$celltype <- factor(df$celltype)
  cells_by_ct <- split(seq_len(n), as.integer(df$celltype))
  ## split() omits empty levels; restore so cells_by_ct has length K_g
  if (length(cells_by_ct) != K_g) {
    out_idx <- vector("list", K_g)
    nm <- match(names(cells_by_ct), seq_len(K_g))
    for (k in seq_along(cells_by_ct)) out_idx[[nm[k]]] <- cells_by_ct[[k]]
    out_idx[vapply(out_idx, is.null, logical(1))] <- list(integer(0))
    cells_by_ct <- out_idx
  }

  ## Initial mu, alpha, tau (per (term, celltype) matrix; scalar/vector input
  ## broadcast across celltypes for backward compatibility)
  mu    <- pmax(Y, 0.5)
  alpha <- rep(1, g_n)
  if (is.null(tau_init)) {
    tau_mat <- matrix(1, nrow = K_t, ncol = K_g,
                       dimnames = list(re$term_levels, re$group_levels))
  } else if (is.matrix(tau_init)) {
    stopifnot(nrow(tau_init) == K_t, ncol(tau_init) == K_g)
    tau_mat <- tau_init
  } else {
    stopifnot(length(tau_init) == K_t)
    tau_mat <- matrix(tau_init, nrow = K_t, ncol = K_g, byrow = FALSE,
                       dimnames = list(re$term_levels, re$group_levels))
  }

  hist <- list(tau = list(), alpha = list(), rel_delta = numeric())

  prev_eta <- log(mu) - offset_vec
  converged <- FALSE

  ## Allocate B, U, re_var, se_* once
  B <- matrix(0, p, g_n)
  U <- matrix(0, q, g_n)
  re_var <- matrix(0, q, g_n)
  se_B <- matrix(NA_real_, p, g_n)
  se_U <- matrix(NA_real_, q, g_n)

  ## Helper: REML M-step for per-(term, celltype) tau:
  ## tau_{t,c} = mean over genes of (u_{tcg}^2 + posterior_var_{tcg})
  ## Each parameter informed by all G genes' BLUPs for that (t, c) combo --
  ## plenty of data to estimate K_t * K_g = 90 parameters from G * K_g
  ## BLUPs per parameter (= G when celltype is fixed).
  .em_tau <- function(U_, V_) {
    out <- matrix(NA_real_, K_t, K_g,
                   dimnames = list(re$term_levels, re$group_levels))
    for (t in seq_len(K_t)) {
      for (c in seq_len(K_g)) {
        col <- (t - 1L) * K_g + c
        out[t, c] <- pmax(mean(U_[col, ]^2 + V_[col, ], na.rm = TRUE), 1e-6)
      }
    }
    out
  }

  ## n_per_group for hierarchical shrinkage (single-block)
  n_per_group_main <- vapply(cells_by_ct, length, 0L)

  prev_alpha <- alpha
  tau_outer_history <- list(tau_mat)

  ## Per-gene tau matrix (q x G). On iter 1 every gene sees the same panel
  ## tau_mat broadcast; from iter 2 we shrink per-gene toward the panel.
  tau_g_array <- matrix(rep(as.numeric(t(tau_mat)), times = g_n),
                         nrow = q, ncol = g_n,
                         dimnames = list(paste0(re$group_of_col, "::",
                                                  re$term_of_col),
                                          NULL))
  if (!is.null(colnames(Y))) colnames(tau_g_array) <- colnames(Y)
  d0_per_tc      <- rep(NA_real_, q)
  lambda_sq_prev <- NULL                    ## carried across outer iters
  a_prev         <- NULL                    ## (half-Cauchy auxiliary state)
  for (it in seq_len(n_iter)) {
    t_it <- Sys.time()
    eta <- log(mu) - offset_vec
    z   <- eta + (Y - mu) / mu
    w   <- sweep(mu, 2, (1 + alpha), "/")

    ## --- One WLS solve per gene with per-gene (or shared) tau ---
    ## Per-gene lam_diag = 1/tau_g_array[, gi]. Under "shared" mode every
    ## column of tau_g_array is identical to the broadcast tau_mat.
    lam_diag_mat <- 1 / tau_g_array  # q x G

    per_gene <- BiocParallel::bplapply(seq_len(g_n), function(gi) {
      lam_g <- lam_diag_mat[, gi]
      sol <- .solve_one_gene_block(X_fixed, X_terms_dense, cells_by_ct,
                                    w[, gi], z[, gi], lam_g, K_t, K_g)
      list(beta    = sol$beta,
           u       = sol$u,
           re_var  = pmax(sol$Ainv_diag[p + seq_len(q)], 0),
           ## Always compute SEs so early-exit at convergence still
           ## yields valid se_B/se_U (mashr filter rejects NA SEs).
           se_beta = sqrt(pmax(sol$Ainv_diag[seq_len(p)], 0)),
           se_u    = sqrt(pmax(sol$Ainv_diag[p + seq_len(q)], 0)))
    }, BPPARAM = BPPARAM)

    for (gi in seq_len(g_n)) {
      B[, gi]      <- per_gene[[gi]]$beta
      U[, gi]      <- per_gene[[gi]]$u
      re_var[, gi] <- per_gene[[gi]]$re_var
      se_B[, gi]   <- per_gene[[gi]]$se_beta
      se_U[, gi]   <- per_gene[[gi]]$se_u
    }

    ## --- Update mu and per-gene alpha (NB1 MLE, with damping for stability) ---
    eta_new <- as.matrix(X_fixed %*% B) + as.matrix(Z %*% U)
    mu      <- pmax(exp(eta_new + offset_vec), 1e-6)

    ## Per-gene alpha is embarrassingly parallel; bplapply across the same
    ## worker pool used for Stage 3 (free since that pool is otherwise idle here).
    alpha_list <- BiocParallel::bplapply(seq_len(g_n),
                    function(gi) .alpha_nb1_mle(Y[, gi], mu[, gi],
                                                  max_n = alpha_max_n),
                    BPPARAM = BPPARAM)
    alpha <- unlist(alpha_list, use.names = FALSE)
    alpha[!is.finite(alpha)] <- prev_alpha[!is.finite(alpha)]
    alpha <- pmin(pmax(alpha, 1e-4), 50)
    prev_alpha <- alpha

    ## --- REML EM step on per-(term, celltype) tau ---
    tau_mat <- .em_tau(U, re_var)
    if (tau_shrinkage == "hierarchical") {
      tau_mat <- .shrink_tau_hierarchical(tau_mat, n_per_group_main)
    }
    tau_outer_history[[length(tau_outer_history) + 1L]] <- tau_mat

    ## --- Per-gene tau update ---
    s2_mat <- U^2 + re_var   # per-gene "data" estimate of variance
    rownames(s2_mat) <- rownames(tau_g_array)
    if (tau_shrinkage == "shared" || tau_shrinkage == "hierarchical") {
      tau_g_array[] <- rep(as.numeric(t(tau_mat)), times = g_n)
    } else if (tau_shrinkage == "adaptive") {
      tau_g_array <- .adaptive_tau_eb(s2_mat, K_t, K_g,
                                       gene_names = colnames(Y))
      d0_per_tc <- attr(tau_g_array, "d0")
    } else if (tau_shrinkage == "half_cauchy") {
      tau_g_array <- .adaptive_tau_half_cauchy(
        s2_mat, K_t, K_g,
        gene_names     = colnames(Y),
        lambda_sq_prev = lambda_sq_prev,
        a_prev         = a_prev,
        n_em_iter      = if (it <= 2L) 8L else 5L
      )
      lambda_sq_prev <- attr(tau_g_array, "lambda_sq")
      a_prev         <- attr(tau_g_array, "a")
    }

    ## Outer convergence on eta
    rel_delta <- max(abs(eta_new - prev_eta) /
                       pmax(abs(prev_eta), 1e-3), na.rm = TRUE)
    prev_eta <- eta_new

    hist$tau[[it]]   <- tau_mat
    hist$alpha[[it]] <- alpha
    hist$rel_delta[it] <- rel_delta

    if (verbose) {
      tau_med_per_term <- apply(tau_mat, 1, stats::median)
      d0_summary <- if (tau_shrinkage == "adaptive") {
        d0_finite <- d0_per_tc[is.finite(d0_per_tc)]
        if (length(d0_finite) > 0) {
          sprintf("  d0[med,p25,p75]=[%.1f,%.1f,%.1f] inf=%d/%d",
                  stats::median(d0_finite),
                  stats::quantile(d0_finite, 0.25),
                  stats::quantile(d0_finite, 0.75),
                  sum(!is.finite(d0_per_tc)), length(d0_per_tc))
        } else " (all d0=Inf -> full shrinkage)"
      } else if (tau_shrinkage == "half_cauchy" && !is.null(lambda_sq_prev)) {
        ## Per-gene tau dispersion: ratio of 90th to 10th percentile across
        ## genes for a representative slope term, indicating how much
        ## per-gene differentiation the half-Cauchy is allowing.
        spread_log <- log(apply(tau_g_array, 1,
                                 function(x) stats::quantile(x, 0.9)) /
                          apply(tau_g_array, 1,
                                 function(x) stats::quantile(x, 0.1)))
        sprintf("  lambda^2[med]=%.3g  per-gene tau spread (log p90/p10) [med,p75]=[%.2f,%.2f]",
                stats::median(lambda_sq_prev, na.rm = TRUE),
                stats::median(spread_log,    na.rm = TRUE),
                stats::quantile(spread_log, 0.75, na.rm = TRUE))
      } else ""
      cat(sprintf("  [mvpql] iter %d  rel_delta=%.3g  alpha[med,p25,p75]=[%.2f,%.2f,%.2f]  tau_med_per_term=[%s]%s  (%.1fs)\n",
                  it, rel_delta,
                  median(alpha),
                  stats::quantile(alpha, 0.25, na.rm=TRUE),
                  stats::quantile(alpha, 0.75, na.rm=TRUE),
                  paste(sprintf("%.3f", tau_med_per_term), collapse=","),
                  d0_summary,
                  as.numeric(difftime(Sys.time(), t_it, units = "secs"))))
    }
    if (it >= 2 && is.finite(rel_delta) && rel_delta < tol) {
      converged <- TRUE
      break
    }
    ## Release per-iter scratch before next iteration to keep peak RSS down
    rm(per_gene); gc(verbose = FALSE)
  }

  rownames(B) <- colnames(X_fixed)
  rownames(U) <- colnames(Z)
  rownames(se_B) <- colnames(X_fixed)
  rownames(se_U) <- colnames(Z)
  if (!is.null(colnames(Y))) {
    colnames(B)    <- colnames(Y)
    colnames(U)    <- colnames(Y)
    colnames(se_B) <- colnames(Y)
    colnames(se_U) <- colnames(Y)
  }

  list(B = B, U = U, se_B = se_B, se_U = se_U,
       alpha = alpha,
       tau_mat        = tau_mat,                          # K_t × K_g panel
       tau_g_array    = tau_g_array,                      # q × G per-gene
       d0_per_tc      = d0_per_tc,                        # length q (adaptive only)
       lambda_sq      = lambda_sq_prev,                   # length q (half_cauchy only)
       tau            = apply(tau_mat, 1, stats::median), # back-compat
       tau_shrinkage  = tau_shrinkage,
       mu      = mu,
       re_meta = re, n_iter = it, converged = converged, history = hist)
}


# Multivariate PQL fit -- multi-block random effects.
#
# Generalises fit_pace_mvpql() to arbitrary `(formula || group_col)` random
# effects. Each block has its own per-(term, group) variance matrix; the EM
# update operates per block, the WLS solve sees the stacked Z. Used for the
# canonical melanoma model:
#   gene ~ 1 + offset + Responder + spillover_near +
#          (1 + Responder * (vars) || celltype) +
#          (1 || imageID)
#
# @param re_specs list of `list(group_col, formula)` blocks. Each formula's
#   model.matrix columns become per-group random terms.
# @return Same structure as fit_pace_mvpql() but with `re_meta` holding a
#   multi-block design (build_random_design_multi output) and `tau_blocks`
#   a list of per-block tau matrices.
fit_pace_mvpql_multi <- function(Y, X_fixed, df, re_specs,
                                  offset_vec = NULL,
                                  n_iter = 6, tol = 5e-3,
                                  tau_shrinkage = c("shared","adaptive","half_cauchy","hierarchical"),
                                  BPPARAM = BiocParallel::SerialParam(),
                                  chunk_size = 256L,
                                  alpha_max_n = Inf,
                                  sample_weight = NULL,
                                  n_threads = NULL,
                                  interior_precision = 0L,
                                  ## Upper bound on the variance components; a binding cap means the
                                  ## term is identified only by the ridge. See fit_pace_mvpql_streaming.
                                  tau_max = 100,
                                  verbose = TRUE) {
  tau_shrinkage <- match.arg(tau_shrinkage)
  if (is.null(offset_vec)) offset_vec <- rep(0, nrow(Y))
  ## Resolve OpenMP thread count for the C++ kernel. Decoupled from BPPARAM
  ## so the caller can use SerialParam (no fork) while still getting
  ## multi-threaded solve inside pace::solve_genes_chunk.
  if (is.null(n_threads))
    n_threads <- tryCatch(max(1L, BiocParallel::bpworkers(BPPARAM)),
                           error = function(e) 1L)
  n_threads <- max(1L, as.integer(n_threads))
  Y <- as.matrix(Y); storage.mode(Y) <- "double"
  n <- nrow(Y); g_n <- ncol(Y); p <- ncol(X_fixed)
  ## Optional per-cell sample weight: multiplies the WLS weight w in the
  ## solve. Use to balance cluster contributions (e.g. cohort_w_i =
  ## 1 / n_cells_in_image_celltype to give each (image, celltype) equal
  ## total weight). Default NULL = no reweighting.
  if (!is.null(sample_weight)) {
    if (length(sample_weight) != n) stop("sample_weight must have length n")
    sample_weight <- as.numeric(sample_weight)
    if (any(!is.finite(sample_weight)) || any(sample_weight < 0))
      stop("sample_weight must be non-negative finite")
  }

  re <- build_random_design_multi(df, re_specs)
  Z  <- re$Z; q <- ncol(Z)
  if (verbose) {
    blk_str <- paste(vapply(re$blocks, function(b)
      sprintf("%s[%dx%d]", b$group_col, b$K_terms, b$K_groups),
      character(1)), collapse = "+")
    cat(sprintf("  [mvpql.multi] n=%d  g=%d  p_fixed=%d  q_random=%d (= %s)\n",
                n, g_n, p, q, blk_str))
  }

  ## Per-block panel tau matrix (term x group)
  tau_blocks <- lapply(re$blocks, function(b) {
    matrix(1, b$K_terms, b$K_groups,
           dimnames = list(b$term_levels, b$group_levels))
  })
  ## tau_g_array: q x G; broadcast tau_blocks into the q-length per-gene vector
  build_tau_vec <- function() {
    out <- numeric(q)
    for (bi in seq_along(re$blocks)) {
      blk <- re$blocks[[bi]]
      ## tau_blocks[[bi]] is K_t_b x K_g_b; col-major flatten matches the
      ## column order in Z (term-major: term1 group1..K, term2 group1..K, ...)
      slice <- as.numeric(t(tau_blocks[[bi]]))
      out[(blk$col_offset + 1L):(blk$col_offset + blk$n_cols)] <- slice
    }
    out
  }
  tau_g_array <- matrix(build_tau_vec(), nrow = q, ncol = g_n)
  rownames(tau_g_array) <- colnames(Z)
  if (!is.null(colnames(Y))) colnames(tau_g_array) <- colnames(Y)

  mu    <- pmax(Y, 0.5)
  alpha <- rep(1, g_n)
  prev_eta   <- log(mu) - offset_vec
  prev_alpha <- alpha
  hist <- list(tau_blocks = list(), alpha = list(), rel_delta = numeric(),
               n_nonfinite = integer(), tau_max_seen = numeric())
  converged <- FALSE

  B <- matrix(0, p, g_n); U <- matrix(0, q, g_n)
  re_var <- matrix(0, q, g_n)
  se_B <- matrix(NA_real_, p, g_n); se_U <- matrix(NA_real_, q, g_n)

  ## Per-block EM update of tau (REML-style: mean of u^2 + V across genes
  ## within each (term, group) of the block).
  .em_tau_blocks <- function(U_, V_) {
    out <- vector("list", length(re$blocks))
    for (bi in seq_along(re$blocks)) {
      blk <- re$blocks[[bi]]
      m <- matrix(NA_real_, blk$K_terms, blk$K_groups,
                  dimnames = list(blk$term_levels, blk$group_levels))
      for (t in seq_len(blk$K_terms)) {
        for (c in seq_len(blk$K_groups)) {
          col <- blk$col_offset + (t - 1L) * blk$K_groups + c
          m[t, c] <- pmax(mean(U_[col, ]^2 + V_[col, ], na.rm = TRUE), 1e-6)
        }
      }
      out[[bi]] <- m
    }
    out
  }

  ## Set once the convergence criterion has been met on a float interior; the
  ## next pass then runs in double so the SEs we keep are not float-quantised.
  final_pass <- FALSE
  for (it in seq_len(n_iter)) {
    t_it <- Sys.time()
    last_iter <- (it == n_iter) || final_pass
    lam_diag_mat <- 1 / tau_g_array

    ## Process genes in chunks. To control peak memory, compute the
    ## working response (z) and weights (w) lazily per chunk rather than
    ## allocating the full cells x genes matrices in the parent process.
    ## At DKD scale (110k cells x ~3000 genes) the full matrices would be
    ## ~7 GB each; per-chunk slices are ~150 MB and freed between chunks.
    ## Fork-time copy-on-write tax on workers drops correspondingly.
    chk_starts <- seq.int(1L, g_n, by = max(1L, as.integer(chunk_size)))
    for (cs in chk_starts) {
      gene_idx_chk <- cs:min(cs + chunk_size - 1L, g_n)
      mu_chk    <- mu[, gene_idx_chk, drop = FALSE]
      eta_chk   <- log(mu_chk) - offset_vec
      z_chk     <- eta_chk + (Y[, gene_idx_chk, drop = FALSE] - mu_chk) / mu_chk
      w_chk     <- sweep(mu_chk, 2, (1 + alpha[gene_idx_chk]), "/")
      if (!is.null(sample_weight)) w_chk <- w_chk * sample_weight
      lam_chk   <- lam_diag_mat[, gene_idx_chk, drop = FALSE]
      rm(eta_chk, mu_chk)
      ## Float-interior schedule: float for non-final iters, double for the
      ## last iter so the SE diagonal inherits double precision. Caller can
      ## force double everywhere with interior_precision = 0.
      iter_precision <- if (last_iter) 0L else as.integer(interior_precision)
      if (iter_precision != 0L && .ridge_needs_double(lam_chk, w_chk))
        iter_precision <- 0L
      per_gene_chk <- .solve_genes_chunk_multiblock(
        X_fixed, re$X_terms_list, re$cell_grp_list, re$cells_by_grp_list,
        w_chk, z_chk, lam_chk, re$blocks,
        gene_idx = seq_along(gene_idx_chk),
        n_threads = n_threads,
        interior_precision = iter_precision,
        BPPARAM = BPPARAM)
      ## Whole blocks: the solve returns its three matrices in gene order.
      ## Always compute SEs so early-exit at convergence still yields valid
      ## se_B/se_U (mashr filter rejects NA SEs), matching the single-block
      ## sibling above. Non-final iterations may run a float interior, so
      ## these are refreshed in double on the last iteration.
      random_rows <- p + seq_len(q)
      B[, gene_idx_chk]      <- per_gene_chk$B
      U[, gene_idx_chk]      <- per_gene_chk$U
      re_var[, gene_idx_chk] <- pmax(per_gene_chk$Ainv_diag[random_rows, , drop = FALSE], 0)
      se_B[, gene_idx_chk]   <- sqrt(pmax(per_gene_chk$Ainv_diag[seq_len(p), , drop = FALSE], 0))
      se_U[, gene_idx_chk]   <- sqrt(pmax(per_gene_chk$Ainv_diag[random_rows, , drop = FALSE], 0))
      rm(per_gene_chk, z_chk, w_chk, lam_chk)
    }

    eta_new <- as.matrix(X_fixed %*% B) + as.matrix(Z %*% U)
    mu      <- pmax(exp(eta_new + offset_vec), 1e-6)
    ## Per-gene alpha is embarrassingly parallel; bplapply across the same
    ## worker pool used for Stage 3 (free since that pool is otherwise idle here).
    alpha_list <- BiocParallel::bplapply(seq_len(g_n),
                    function(gi) .alpha_nb1_mle(Y[, gi], mu[, gi],
                                                  max_n = alpha_max_n),
                    BPPARAM = BPPARAM)
    alpha <- unlist(alpha_list, use.names = FALSE)
    alpha[!is.finite(alpha)] <- prev_alpha[!is.finite(alpha)]
    alpha <- pmin(pmax(alpha, 1e-4), 50)
    prev_alpha <- alpha

    tau_blocks <- .em_tau_blocks(U, re_var)
    if (tau_shrinkage == "hierarchical") {
      for (bi in seq_along(tau_blocks)) {
        n_c_b <- vapply(re$cells_by_grp_list[[bi]], length, 0L)
        tau_blocks[[bi]] <- .shrink_tau_hierarchical(tau_blocks[[bi]], n_c_b)
      }
    }
    if (tau_shrinkage == "shared" || tau_shrinkage == "hierarchical") {
      tau_g_array[] <- rep(build_tau_vec(), times = g_n)
    } else {
      ## Per-block adaptive shrinkage: run .adaptive_tau_eb on each block's
      ## sub-matrix of s2_mat.
      s2_mat <- U^2 + re_var
      rownames(s2_mat) <- rownames(tau_g_array)
      for (bi in seq_along(re$blocks)) {
        blk <- re$blocks[[bi]]
        rng <- (blk$col_offset + 1L):(blk$col_offset + blk$n_cols)
        s2_b <- s2_mat[rng, , drop = FALSE]
        tau_b_g <- if (tau_shrinkage == "adaptive") {
          .adaptive_tau_eb(s2_b, blk$K_terms, blk$K_groups,
                           gene_names = colnames(Y))
        } else {
          .adaptive_tau_half_cauchy(s2_b, blk$K_terms, blk$K_groups,
                                     gene_names = colnames(Y))
        }
        tau_g_array[rng, ] <- tau_b_g
      }
    }
    tau_g_array <- .clamp_tau(tau_g_array, tau_max, it)
    hist$tau_max_seen[it] <- max(tau_g_array, na.rm = TRUE)

    rd_mat    <- abs(eta_new - prev_eta) / pmax(abs(prev_eta), 1e-3)
    ## na.rm below hides every gene whose linear predictor has gone non-finite,
    ## so count them and refuse to declare convergence while any remain.
    n_nonfinite <- sum(!is.finite(rd_mat))
    rel_delta   <- max(rd_mat, na.rm = TRUE)
    rm(rd_mat)
    prev_eta <- eta_new
    hist$n_nonfinite[it]  <- n_nonfinite
    hist$tau_blocks[[it]] <- tau_blocks
    hist$alpha[[it]]      <- alpha
    hist$rel_delta[it]    <- rel_delta

    if (verbose) {
      tau_med_per_block <- vapply(tau_blocks, function(m)
                                    stats::median(m), numeric(1))
      cat(sprintf("  [mvpql.multi] iter %d  rel_delta=%.3g  alpha[med]=%.2f  tau_med_per_block=[%s] (max %.3g)  (%.1fs)\n",
                  it, rel_delta, median(alpha),
                  paste(sprintf("%.3f", tau_med_per_block), collapse=","),
                  hist$tau_max_seen[it],
                  as.numeric(difftime(Sys.time(), t_it, units = "secs"))))
    }
    if (n_nonfinite > 0L) {
      warning(sprintf(
        "iter %d: %d gene-cell linear predictors are non-finite; these are excluded from the convergence metric and the fit will not be reported as converged.",
        it, n_nonfinite), call. = FALSE)
    }
    if (it >= 2 && is.finite(rel_delta) && rel_delta < tol && n_nonfinite == 0L) {
      if (final_pass || as.integer(interior_precision) == 0L) {
        converged <- TRUE; break
      }
      ## Converged on a float interior: take one more pass in double so the
      ## retained SEs come from an unquantised ridge, then stop.
      final_pass <- TRUE
    }
    ## Release per-iter scratch before next iteration to keep peak RSS down
    if (exists("per_gene", inherits = FALSE)) {
      rm(per_gene); gc(verbose = FALSE)
    }
  }

  rownames(B)    <- colnames(X_fixed); rownames(U)    <- colnames(Z)
  rownames(se_B) <- colnames(X_fixed); rownames(se_U) <- colnames(Z)
  if (!is.null(colnames(Y))) {
    colnames(B) <- colnames(U) <- colnames(se_B) <- colnames(se_U) <- colnames(Y)
  }

  list(B = B, U = U, se_B = se_B, se_U = se_U,
       alpha          = alpha,
       tau_blocks     = tau_blocks,
       tau_g_array    = tau_g_array,
       tau_shrinkage  = tau_shrinkage,
       mu             = mu,
       re_meta        = re,
       n_iter         = it, converged = converged, history = hist)
}


# Multivariate, per-gene Nakagawa-style variance decomposition for MV-PQL.
#
# Produces a per-gene + per-(gene, focal) variance decomposition without
# refitting -- entirely from the joint MV-PQL posterior. Each gene gets its
# OWN effective variance components (\eqn{\hat\tau_{t,g}}), but those
# components are computed using both the shrunken BLUPs AND their posterior
# variances (the diagonal of \eqn{(D^\top W_g D + \Lambda)^{-1}}), which
# carry the multivariate information from the shared-tau prior. So the
# result is "multivariate-in-prior, per-gene-in-output".
#
# Per-gene panel-level decomposition (link scale):
#   \deqn{\hat\tau_{t,g} = \frac{1}{K_g} \sum_c (u_{tcg}^2 + \mathrm{Var}(u_{tcg} \mid y_g))}
#   \deqn{\sigma^2_{\mathrm{fix}}(g) = \mathrm{Var}_i(X_i \beta_g)}
#   \deqn{\sigma^2_{\mathrm{disp}}(g) = \log(1 + 1/\bar\mu_g + \alpha_g)} (NB1)
#   \deqn{\mathrm{share}_t(g) = \hat\tau_{t,g} / \mathrm{total}(g)}
#
# Per-(gene, focal) decomposition uses the same link-scale block sums as
# coef_ss_link but corrects each ranef term with the posterior variance:
#   \deqn{\mathrm{Cell type}~ \mathrm{SS}(g, c) = (u_{(\mathrm{Int}), c, g}^2 + V_{(\mathrm{Int}), c, g}) \cdot n_c}
#   \deqn{\mathrm{Spatial state}~ \mathrm{SS}(g, c) = \sum_{t \in \mathrm{vars}} (u_{tcg}^2 + V_{tcg}) \cdot \sum_{i \in c} N_{i,t}^2}
#
# @param fit Output of fit_pace_mvpql()
# @param df cell-level df (must contain celltype + vars columns)
# @param Y cells x genes count matrix used in the fit (for per-cell mu_bar
#   and per-focal mu_bar)
# @param vars character vector of slope-term names (= the random-slope
#   covariates passed to fit_pace_mvpql)
# @param X_fixed cells x p fixed-effect design matrix (the same one passed
#   to fit_pace_mvpql); used to compute sigma^2_fix(g) and the spillover
#   block contribution
# @param focal_levels which celltypes to evaluate per-focal decomposition for
#   (default: all levels of df$celltype)
# @return list with:
#   - gene_total: per-gene panel-level decomposition tibble
#   - gene_focal: per-(gene, focal) decomposition tibble
#   - agg_mean / agg_median: per-focal aggregate share tables (matching
#     canonical nakagawa_decomp_full.rds schema)
mvpql_variance_decomposition <- function(fit, df, Y, vars, X_fixed,
                                          focal_levels = NULL) {
  re <- fit$re_meta
  K_t <- re$K_terms; K_g <- re$K_groups
  groups <- re$group_levels; terms <- re$term_levels
  gene_names <- colnames(fit$U)
  G <- length(gene_names)
  if (is.null(focal_levels))
    focal_levels <- as.character(unique(df$celltype))

  ## Helper: column index in U/se_U for (term t, group g)
  col_idx <- function(t, g) (t - 1L) * K_g + g

  ## ---- Pre-compute per-cell helpers ----
  cell_ct <- as.character(df$celltype)
  ct_to_idx <- match(cell_ct, groups)

  ## Pre-compute Σ_{i in c} N_{i,t}^2 for each (term, focal_celltype)
  ## (plus n_c for the intercept block).
  ssN <- array(0, dim = c(K_t, K_g),
               dimnames = list(terms, groups))
  for (g_idx in seq_along(groups)) {
    cells_in_g <- which(ct_to_idx == g_idx)
    if (!length(cells_in_g)) next
    ssN[1L, g_idx] <- length(cells_in_g)  ## intercept SS multiplier = n_c
    for (t in seq_len(length(vars))) {
      vname <- vars[t]
      vals  <- df[[vname]][cells_in_g]
      ssN[t + 1L, g_idx] <- sum(vals^2, na.rm = TRUE)
    }
  }

  ## ---- Per-gene panel-level decomposition ----
  per_gene <- vector("list", G)
  for (gi in seq_len(G)) {
    g <- gene_names[gi]

    ## sigma^2_fix on the link scale (variance of X_i' beta_g across cells)
    eta_fix <- as.numeric(X_fixed %*% fit$B[, gi])
    sigma2_fix <- stats::var(eta_fix, na.rm = TRUE)

    ## EM-corrected per-term variance components (multivariate posterior)
    tau_g <- numeric(K_t); names(tau_g) <- terms
    for (t in seq_len(K_t)) {
      cols_t <- ((t - 1L) * K_g + 1L):(t * K_g)
      U_tg   <- fit$U[cols_t, gi]
      V_tg   <- fit$se_U[cols_t, gi]^2
      tau_g[t] <- mean(U_tg^2 + V_tg, na.rm = TRUE)
    }

    ## NB1 link-scale dispersion (Nakagawa et al. 2017)
    mu_g     <- fit$mu[, gi]
    mu_bar   <- mean(mu_g, na.rm = TRUE)
    alpha_g  <- unname(fit$alpha[gi])
    sigma2_disp <- log(1 + 1 / max(mu_bar, 1e-9) + max(alpha_g, 0))

    sigma2_intercept <- tau_g[1L]
    sigma2_state     <- sum(tau_g[-1L])
    total            <- sigma2_fix + sigma2_intercept + sigma2_state + sigma2_disp

    per_gene[[gi]] <- tibble::tibble(
      gene             = g,
      sigma2_fix       = sigma2_fix,
      sigma2_intercept = sigma2_intercept,
      sigma2_state     = sigma2_state,
      sigma2_disp      = sigma2_disp,
      sigma2_total     = total,
      `Cell type %`    = sigma2_intercept / total * 100,
      `Spatial state %`= sigma2_state    / total * 100,
      `Spillover %`    = sigma2_fix      / total * 100,
      `Residual %`     = sigma2_disp     / total * 100,
      R2_marg          = sigma2_fix      / total,
      R2_cond          = (sigma2_fix + sigma2_intercept + sigma2_state) / total,
      mu_bar           = mu_bar,
      alpha            = alpha_g
    )
  }
  gene_total <- dplyr::bind_rows(per_gene)

  ## ---- Per-(gene, focal) decomposition (link-scale SS shares) ----
  ## Match canonical block_sum_coef_ss conventions:
  ##   Cell type   = ranef intercept BLUP^2 (+ posterior var) * n_c
  ##   Spatial state = sum_t ranef slope BLUP^2 (+ posterior var) * Sigma N_{i,t}^2
  ##   Spillover   = sum_t fixed-effect beta_t^2 * Sigma N_{i,t}^2  (per-term diag)
  ##                 NOT (X beta)^2 -- that mixes intercept + cross terms which
  ##                 inflates spillover and breaks comparability with canonical.
  ##   Residuals   = working-residual SS + (intercept fixed-effect SS) +
  ##                 anything not categorized

  ## Identify which fixed-effect coefficients are spillover (_near) vs intercept
  fix_names <- rownames(fit$B)
  spill_idx <- which(grepl("_near$|spill", fix_names, ignore.case = TRUE))
  int_idx   <- which(fix_names == "(Intercept)")
  ## Unified spillover (2026-05-18+): the single celltype x celltype matrix
  ## phi_mat carries all spillover via fit$technical_offset_mat (n x G, the
  ## per-cell-per-gene bleed log-offset Σ_c' phi[c,c'] κ(g,c') N_c'(cell)).
  ## When present & nonzero this REPLACES the (now-removed) _near covariate
  ## as the Spillover block; legacy fits fall back to the _near path.
  use_bleed <- !is.null(fit$technical_offset_mat) &&
               any(fit$technical_offset_mat != 0, na.rm = TRUE)
  ## Per-term Sigma X^2 per focal celltype (n_c for intercept, ssN for the
  ## spillover terms which are *_near columns also present in df)
  fixed_term_ss_per_focal <- function(c_idx) {
    setNames(numeric(length(fix_names)), fix_names)
  }
  ## We can compute spillover term SS once per (focal, term): sum_{i in c} X_{i, term}^2
  spill_ss_per_focal <- matrix(0, nrow = length(spill_idx), ncol = length(groups),
                                dimnames = list(fix_names[spill_idx], groups))
  for (c_idx in seq_along(groups)) {
    cells_in_g <- which(ct_to_idx == c_idx)
    if (!length(cells_in_g)) next
    Xc <- X_fixed[cells_in_g, spill_idx, drop = FALSE]
    spill_ss_per_focal[, c_idx] <- colSums(Xc^2, na.rm = TRUE)
  }

  rows <- vector("list", G * length(focal_levels))
  ix <- 0L
  for (c_name in focal_levels) {
    c_idx <- match(c_name, groups)
    cells_in_c <- which(ct_to_idx == c_idx)
    if (!length(cells_in_c)) next
    n_c <- length(cells_in_c)

    for (gi in seq_len(G)) {
      g <- gene_names[gi]

      ## Cell type SS: (BLUP_int_c^2 + V_int_c) * n_c
      u_int <- fit$U[col_idx(1L, c_idx), gi]
      v_int <- fit$se_U[col_idx(1L, c_idx), gi]^2
      ss_celltype <- (u_int^2 + v_int) * n_c

      ## Spatial state SS: sum over slope terms of (BLUP^2 + V) * Σ N_{i,t}^2
      ss_state <- 0
      for (t in seq_along(vars)) {
        ut <- fit$U[col_idx(t + 1L, c_idx), gi]
        vt <- fit$se_U[col_idx(t + 1L, c_idx), gi]^2
        ss_state <- ss_state + (ut^2 + vt) * ssN[t + 1L, c_idx]
      }

      ## Spillover SS: from the phi bleed offset (centered SS within focal)
      ## when unified; else legacy per-term _near diag.
      if (use_bleed) {
        bo <- fit$technical_offset_mat[cells_in_c, gi]
        ss_spill <- sum((bo - mean(bo, na.rm = TRUE))^2, na.rm = TRUE)
      } else {
        beta_spill <- fit$B[spill_idx, gi]
        ss_spill   <- sum(beta_spill^2 * spill_ss_per_focal[, c_idx], na.rm = TRUE)
      }

      ## Residual SS: working residuals + intercept fixed contribution
      mu_c <- fit$mu[cells_in_c, gi]
      y_c  <- as.numeric(Y[cells_in_c, gi])
      r    <- (y_c - mu_c) / pmax(mu_c, 1e-9)
      ss_working_resid <- sum(r^2, na.rm = TRUE)
      ss_int  <- if (length(int_idx)) (fit$B[int_idx, gi])^2 * n_c else 0
      ss_resid <- ss_working_resid + ss_int

      ss_total <- ss_celltype + ss_state + ss_spill + ss_resid
      if (ss_total <= 0) next

      rows[[ix <- ix + 1L]] <- tibble::tibble(
        gene  = g,
        focal = c_name,
        n_focal = n_c,
        SS_celltype  = ss_celltype,
        SS_state     = ss_state,
        SS_spillover = ss_spill,
        SS_residual  = ss_resid,
        SS_total     = ss_total,
        `Cell type %`    = ss_celltype  / ss_total * 100,
        `Spatial state %`= ss_state     / ss_total * 100,
        `Spillover %`    = ss_spill     / ss_total * 100,
        `Residual %`     = ss_resid     / ss_total * 100
      )
    }
  }
  gene_focal <- dplyr::bind_rows(rows[seq_len(ix)])

  ## ---- Per-(gene, focal) Nakagawa-style within-focal variance partition ----
  ## Asks: for cells of focal celltype c, what fraction of the link-scale
  ## VARIANCE of eta is explained by spatial-state (slope BLUPs x neighbour
  ## counts), spillover (fixed-effect slopes x neighbour counts), and the
  ## NB1 link-scale residual.
  ##
  ## The celltype intercept BLUP is constant within focal cells, so it has
  ## zero within-focal variance; we report its squared link-scale offset
  ## separately as `celltype_offset_sq` so the user can see how much
  ## celltype-baseline shift exists for each gene.
  ##
  ## Contamination annotation: per (gene, focal) we also compute Tan-style
  ## specificity = focal_mean / (focal_mean + max_other_mean) and the raw
  ## ratio focal/max_other. Genes with ratio < 1 are more highly expressed
  ## in some other celltype than in the focal -- these are likely showing
  ## up due to neighbour-cell bleed-through into the focal segmentation
  ## (the canonical MCSD-varpart pipeline downweights them via spec^2).
  ## Annotated as `is_contaminated`; we leave them in the table for
  ## transparency but flag for downstream filtering.

  ## Pre-compute per-(gene, celltype) mean expression for specificity.
  ct_means <- matrix(0, length(groups), G,
                     dimnames = list(groups, gene_names))
  for (c_idx in seq_along(groups)) {
    cells_c <- which(ct_to_idx == c_idx)
    if (length(cells_c) > 0) {
      ct_means[c_idx, ] <- colMeans(Y[cells_c, , drop = FALSE], na.rm = TRUE)
    }
  }

  rows_naka <- vector("list", G * length(focal_levels))
  ix2 <- 0L
  for (c_name in focal_levels) {
    c_idx <- match(c_name, groups)
    cells_in_c <- which(ct_to_idx == c_idx)
    if (length(cells_in_c) < 5L) next
    n_c <- length(cells_in_c)

    ## per-focal neighbour-count design matrix
    N_c <- as.matrix(df[cells_in_c, vars, drop = FALSE])

    ## fixed-effect spillover beta indices (skip intercept)
    beta_spill_names <- fix_names[spill_idx]

    ## Per-focal specificity for ALL genes (vectorised once per focal)
    ## focal_mean / (focal_mean + max_other_mean), Tan style
    fmean      <- ct_means[c_idx, ]
    other_max  <- apply(ct_means[-c_idx, , drop = FALSE], 2, max, na.rm = TRUE)
    spec_focal <- fmean / pmax(fmean + other_max, 1e-9)
    foratio    <- fmean / pmax(other_max, 1e-9)

    for (gi in seq_len(G)) {
      g <- gene_names[gi]

      ## Per-focal slope BLUPs (length = length(vars))
      slope_blups <- vapply(seq_along(vars),
                             function(t) fit$U[col_idx(t + 1L, c_idx), gi],
                             numeric(1))
      slope_post  <- vapply(seq_along(vars),
                             function(t) fit$se_U[col_idx(t + 1L, c_idx), gi]^2,
                             numeric(1))

      ## eta_state contribution per focal cell (link scale)
      eta_state <- as.numeric(N_c %*% slope_blups)
      ## EM-corrected within-focal variance: Var(eta_state) + Σ_t V_t * Var(N_{i,t})
      var_N <- apply(N_c, 2, stats::var, na.rm = TRUE)
      V_state <- stats::var(eta_state, na.rm = TRUE) + sum(slope_post * var_N, na.rm = TRUE)

      ## eta_spillover: compute against the actual *_near design columns
      ## via X_fixed, NOT the kernel-weighted *_neighbour columns. Using
      ## the wrong design matrix here was masked when both columns had the
      ## same scale (canonical binary) but blows up under kernel features.
      if (use_bleed) {
        V_spill <- stats::var(fit$technical_offset_mat[cells_in_c, gi],
                              na.rm = TRUE)
      } else if (length(spill_idx)) {
        X_spill_c <- X_fixed[cells_in_c, spill_idx, drop = FALSE]
        eta_spill <- as.numeric(X_spill_c %*% fit$B[spill_idx, gi])
        V_spill   <- stats::var(eta_spill, na.rm = TRUE)
      } else {
        V_spill <- 0
      }

      ## NB1 link-scale residual within focal cells
      mu_c    <- fit$mu[cells_in_c, gi]
      mu_bar_c <- mean(mu_c, na.rm = TRUE)
      alpha_g <- unname(fit$alpha[gi])
      V_disp  <- log(1 + 1 / max(mu_bar_c, 1e-9) + max(alpha_g, 0))

      ## Celltype offset (constant within focal) reported separately
      ct_offset_sq <- fit$U[col_idx(1L, c_idx), gi]^2 +
                       fit$se_U[col_idx(1L, c_idx), gi]^2

      tot <- V_state + V_spill + V_disp
      if (!is.finite(tot) || tot <= 0) next

      ## All three percentages are RAW within-focal variance shares (sum
      ## to 100). spec is reported as a separate column for downstream
      ## filtering / weighting; spec^2 weighting is applied in the MCSD
      ## scoring pipeline (mcsd_robust_no_resp_varpart) NOT inside the
      ## variance decomposition. Keeps the decomposition table internally
      ## consistent (percentages sum to 100).
      spec_g <- if (is.finite(spec_focal[gi])) spec_focal[gi] else 0

      rows_naka[[ix2 <- ix2 + 1L]] <- tibble::tibble(
        gene = g, focal = c_name, n_focal = n_c,
        V_state = V_state, V_spill = V_spill, V_disp = V_disp,
        V_total_naka = tot,
        celltype_offset_sq = ct_offset_sq,
        focal_mean        = fmean[gi],
        max_other_mean    = other_max[gi],
        spec              = spec_g,
        focal_other_ratio = foratio[gi],
        is_contaminated   = !is.finite(foratio[gi]) || foratio[gi] < 1,
        ## Raw within-focal variance shares -- sum to 100%.
        `Spatial state %` = V_state / tot * 100,
        `Spillover %`     = V_spill / tot * 100,
        `Residual %`      = V_disp  / tot * 100
      )
    }
  }
  gene_focal_nakagawa <- dplyr::bind_rows(rows_naka[seq_len(ix2)])

  agg_focal_naka_mean <- gene_focal_nakagawa |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      `Spatial state %` = mean(`Spatial state %`, na.rm = TRUE),
      `Spillover %`     = mean(`Spillover %`,     na.rm = TRUE),
      `Residual %`      = mean(`Residual %`,      na.rm = TRUE),
      median_celltype_offset_sq = stats::median(celltype_offset_sq, na.rm = TRUE),
      n_genes           = dplyr::n(),
      .groups = "drop"
    )
  agg_focal_naka_median <- gene_focal_nakagawa |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      `Spatial state %` = stats::median(`Spatial state %`, na.rm = TRUE),
      `Spillover %`     = stats::median(`Spillover %`,     na.rm = TRUE),
      `Residual %`      = stats::median(`Residual %`,      na.rm = TRUE),
      .groups = "drop"
    )

  ## ---- 4-block hybrid decomposition (matches canonical schema) ----
  ## Cell type baseline = (BLUP_int_c^2 + V) treated as a link-scale variance
  ## contribution to the gene's total variability budget. The other three
  ## components are within-focal variances (Nakagawa-style). All four sum to
  ## 100% per (gene, focal), giving the canonical 4-bar stacked layout.
  gene_focal_4block <- gene_focal_nakagawa |>
    dplyr::mutate(
      total_4 = celltype_offset_sq + V_state + V_spill + V_disp,
      `Cell type %`     = celltype_offset_sq / total_4 * 100,
      `Spatial state %` = V_state            / total_4 * 100,
      `Spillover %`     = V_spill            / total_4 * 100,
      `Residual %`      = V_disp             / total_4 * 100
    ) |>
    dplyr::select(gene, focal, n_focal,
                  `Cell type %`, `Spatial state %`,
                  `Spillover %`, `Residual %`,
                  total_4)

  agg_focal_4block_mean <- gene_focal_4block |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      `Cell type %`     = mean(`Cell type %`,     na.rm = TRUE),
      `Spatial state %` = mean(`Spatial state %`, na.rm = TRUE),
      `Spillover %`     = mean(`Spillover %`,     na.rm = TRUE),
      `Residual %`      = mean(`Residual %`,      na.rm = TRUE),
      n_genes           = dplyr::n(),
      .groups = "drop"
    )
  agg_focal_4block_median <- gene_focal_4block |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      `Cell type %`     = stats::median(`Cell type %`,     na.rm = TRUE),
      `Spatial state %` = stats::median(`Spatial state %`, na.rm = TRUE),
      `Spillover %`     = stats::median(`Spillover %`,     na.rm = TRUE),
      `Residual %`      = stats::median(`Residual %`,      na.rm = TRUE),
      .groups = "drop"
    )

  ## ---- Per-focal aggregates ----
  agg_mean <- gene_focal |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      `Cell type %`     = mean(`Cell type %`,    na.rm = TRUE),
      `Spatial state %` = mean(`Spatial state %`, na.rm = TRUE),
      `Spillover %`     = mean(`Spillover %`,    na.rm = TRUE),
      `Residual %`      = mean(`Residual %`,     na.rm = TRUE),
      n_genes           = dplyr::n(),
      .groups = "drop"
    )
  agg_median <- gene_focal |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      `Cell type %`     = stats::median(`Cell type %`,    na.rm = TRUE),
      `Spatial state %` = stats::median(`Spatial state %`, na.rm = TRUE),
      `Spillover %`     = stats::median(`Spillover %`,    na.rm = TRUE),
      `Residual %`      = stats::median(`Residual %`,     na.rm = TRUE),
      .groups = "drop"
    )

  list(
    gene_total              = gene_total,
    gene_focal              = gene_focal,
    agg_mean                = agg_mean,
    agg_median              = agg_median,
    gene_focal_nakagawa     = gene_focal_nakagawa,
    agg_focal_naka_mean     = agg_focal_naka_mean,
    agg_focal_naka_median   = agg_focal_naka_median,
    gene_focal_4block       = gene_focal_4block,
    agg_focal_4block_mean   = agg_focal_4block_mean,
    agg_focal_4block_median = agg_focal_4block_median
  )
}


# Wrap one gene's mv-PQL output as a fit_light object compatible with the
# canonical PACE variance-decomposition consumers (01-core-ss.R::coef_ss_link,
# 03-multigene.R::run_block_decomp_no_resp).
#
# Output schema mirrors what extractRandomEffects() builds in 00-extract.R::
#   - used_rows: row indices of df used in the fit (we use all rows by
#     default since MV-PQL was fit on the full df)
#   - fixef: named numeric vector of fixed-effect coefficients
#   - ranef: named list of (group -> celltype x term BLUP matrix). Only
#     celltype is populated for the no-responder MV-PQL fit.
#   - residuals$working: per-cell working residuals = (y - mu)/mu
#   - sigma_disp: per-gene NB1 dispersion alpha
#   - mu_summary: list(mean = mean(mu), by_celltype = tapply(mu, celltype, mean))
#   - pdHess = TRUE (no boundary-singular issues by construction)
#
# @param fit Output of fit_pace_mvpql()
# @param df Original cell-level df (must include the celltype factor + the
#   neighbour-count covariates used in fitting).
# @param Y cells x genes integer count matrix passed into fit_pace_mvpql()
#   (for residual computation; can be reconstructed from df[, gene_names])
# @param gene gene name (must be a column of fit$U)
# @return fit_light object
mvpql_to_fit_light <- function(fit, df, Y, gene) {
  if (!gene %in% colnames(fit$U))
    stop(sprintf("gene '%s' not in MV-PQL fit", gene))
  re <- fit$re_meta

  ## fit$mu and fit$alpha may be unnamed; resolve by position via U colnames.
  gene_idx <- match(gene, colnames(fit$U))

  ## Reshape per-gene U vector (length q = K_t * K_g) into celltype x term
  ## matrix matching the lme4::ranef(fit)$cond[[group]] schema.
  K_g <- re$K_groups; K_t <- re$K_terms
  ranef_mat <- matrix(NA_real_, K_g, K_t,
                      dimnames = list(re$group_levels, re$term_levels))
  for (t in seq_len(K_t)) {
    for (g in seq_len(K_g)) {
      col <- (t - 1L) * K_g + g
      ranef_mat[g, t] <- fit$U[col, gene_idx]
    }
  }

  fixef_g <- as.numeric(fit$B[, gene_idx])
  names(fixef_g) <- rownames(fit$B)

  mu_g <- as.numeric(fit$mu[, gene_idx])
  y_g  <- as.numeric(Y[, gene])
  ## Working residual on link scale (matches glmmTMB residuals(type="working"))
  resid_working <- (y_g - mu_g) / pmax(mu_g, 1e-9)

  mu_summary <- list(
    mean        = mean(mu_g, na.rm = TRUE),
    by_celltype = if ("celltype" %in% names(df)) {
                    tapply(mu_g, as.character(df[["celltype"]]),
                           mean, na.rm = TRUE)
                  } else NULL
  )

  fl <- list(
    gene       = gene,
    used_rows  = as.character(seq_len(nrow(df))),
    n_obs      = nrow(df),
    fixef      = fixef_g,
    ranef      = list(celltype = ranef_mat),
    pdHess     = TRUE,
    sigma_disp = unname(fit$alpha[gene_idx]),
    varcorr    = NULL,
    mu_summary = mu_summary,
    residuals  = list(working = resid_working)
  )
  class(fl) <- c("fit_light", "list")
  fl
}


# Build the canonical fits_light list (one fit_light per gene) from a
# single MV-PQL fit. Drops directly into run_block_decomp_no_resp(fits, df).
mvpql_to_fits_light <- function(fit, df, Y) {
  ## Force df rownames to "1"..nrow(df) so used_rows lookups work
  rownames(df) <- as.character(seq_len(nrow(df)))
  out <- lapply(colnames(fit$U), function(g) mvpql_to_fit_light(fit, df, Y, g))
  names(out) <- colnames(fit$U)
  out
}


# Wrap mv-PQL output as a `results`-shaped named list (one entry per gene)
# that the canonical PACE downstream consumers (mashr_pipeline, MCSD)
# accept directly. Each entry is `list(ran_vals, fit_light, gene, error)`.
#
# The fit_light here is intentionally minimal -- only `pdHess = TRUE` so
# that `is_well_fit()` from 00-extract.R returns TRUE -- because the
# canonical variance-decomposition consumer reads `VarCorr$cond`,
# `sigma_disp`, and `mu_summary`, which mv-PQL doesn't currently produce
# per-gene (tau is shared across genes by design).
#
# @param fit Output of fit_pace_mvpql()
# @param group_col Name to put in ran_vals$group (default "celltype")
# @return Named list compatible with helpers/mashr_pipeline.R::apply_mashr_shrinkage
mvpql_to_results <- function(fit, group_col = "celltype") {
  rv_long <- mvpql_to_ran_vals(fit, group_col = group_col)
  by_gene <- split(rv_long, rv_long$gene)
  out <- lapply(names(by_gene), function(g) {
    list(
      ran_vals  = by_gene[[g]],
      fit_light = list(pdHess = TRUE),
      gene      = g,
      error     = FALSE
    )
  })
  names(out) <- names(by_gene)
  out
}


# Multi-block ran_vals: emit one row per (gene, block, group, level, term).
# For mashr downstream, filter to the celltype block.
mvpql_to_ran_vals_multi <- function(fit) {
  re <- fit$re_meta
  g_n <- ncol(fit$U)
  gene_names <- colnames(fit$U)
  if (is.null(gene_names)) gene_names <- paste0("g", seq_len(g_n))
  out <- vector("list", g_n)
  for (gi in seq_len(g_n)) {
    est <- fit$U[, gi]; se <- fit$se_U[, gi]
    sca <- est / se
    out[[gi]] <- tibble::tibble(
      effect    = "ran_vals",
      component = "cond",
      group     = re$block_name_of_col,
      level     = re$group_of_col,
      term      = re$term_of_col,
      estimate  = est,
      std.error = se,
      gene      = gene_names[gi],
      lower     = est - 2 * se,
      upper     = est + 2 * se,
      scaled_estimate = sca,
      pval      = 2 * stats::pnorm(-abs(sca))
    )
  }
  do.call(rbind, out)
}


# Multi-block `results` shaped object compatible with mashr_pipeline.
# Filters to the celltype block by default; mashr only consumes celltype REs.
mvpql_to_results_multi <- function(fit, keep_block = "celltype") {
  rv <- mvpql_to_ran_vals_multi(fit)
  rv <- rv[rv$group == keep_block, , drop = FALSE]
  rv$group <- "celltype"  # canonical key expected downstream
  by_gene <- split(rv, rv$gene)
  out <- lapply(names(by_gene), function(g)
    list(ran_vals = by_gene[[g]], fit_light = list(pdHess = TRUE),
         gene = g, error = FALSE))
  names(out) <- names(by_gene)
  out
}


# Multi-block variance decomposition. Uses the celltype-block random
# slopes for the spatial-state contribution, the imageID-block random
# intercept for an "image %" component, and the spillover_near columns
# of X_fixed for spillover. Per-(gene, focal) within-focal variance
# partition (Nakagawa style).
#
# @param fit Output of fit_pace_mvpql_multi()
# @param df cell-level df
# @param Y cells x genes count matrix
# @param vars character vector of neighbour-count covariates (the names
#   of the spatial slopes within the celltype block; used to compute the
#   eta_state variance from the relevant U columns)
# @param X_fixed cells x p fixed-effect design (must include the
#   spillover_near columns and may include Responder)
# @param resp_term name of the Responder term in the celltype-block
#   model.matrix (e.g. "ResponderPD"); set NULL if no Responder.
# @param focal_levels which celltypes to evaluate (default all in df$celltype)
# @param disp_model character; "nb1" (default, matches PACE-MV's NB1 dispersion
#   per `project_pace_overview.md`) uses Leckie et al. 2020 latent-scale
#   formula `log(1 + (1+alpha)/mu_bar)`; "nb2" uses the older Nakagawa-Johnson-
#   Schielzeth 2017 formula `log(1 + 1/mu_bar + alpha)`. The two coincide
#   asymptotically; at small mu_bar (sparse focals) the NB1 formula gives a
#   larger V_disp denominator, deflating the higher-level percentages.
# @param weight_by_spec_sq logical; if TRUE, additionally compute a per-focal
#   aggregate where each gene's percentage is weighted by gene-celltype
#   specificity squared (matches the MCSD spec_exponent=2 weighting). Returned
#   as `agg_focal_5block_specw` in the output list. This is the correct
#   manuscript-headline statistic when genes are panel-targeted and many are
#   non-marker contaminants for the focal celltype.
mvpql_variance_decomposition_multi <- function(fit, df, Y, vars, X_fixed,
                                                 resp_term = NULL,
                                                 focal_levels = NULL,
                                                 disp_model = c("nb1", "nb2"),
                                                 weight_by_spec_sq = TRUE) {
  disp_model <- match.arg(disp_model)
  re <- fit$re_meta
  ## Locate the celltype and image blocks
  blk_idx_ct  <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "celltype")
  blk_idx_img <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "imageID")
  if (length(blk_idx_ct) != 1L)
    stop("Expected exactly one celltype RE block")
  blk_ct  <- re$blocks[[blk_idx_ct]]
  groups  <- blk_ct$group_levels
  if (is.null(focal_levels)) focal_levels <- intersect(as.character(unique(df$celltype)), groups)
  gene_names <- colnames(fit$U); G <- length(gene_names)

  ## Identify column indices for the celltype block, broken out by term
  ct_col <- function(t_idx, g_idx) {
    blk_ct$col_offset + (t_idx - 1L) * blk_ct$K_groups + g_idx
  }
  ## Map term name -> term index within celltype block
  term2t <- setNames(seq_along(blk_ct$term_levels), blk_ct$term_levels)
  ## Spatial-slope term names within the celltype block: each `vars` entry
  ## is one term; if Responder is in the model, also Responder:vars (we
  ## sum BOTH into the within-focal eta_state variance because the model
  ## predicts gene = ... + (no_resp_slope + resp_slope * Responder) * neighbour.
  ## For the per-focal partition we use the actual fitted u values * the
  ## actual neighbour value, which already incorporates Responder.)
  has_resp <- !is.null(resp_term) && any(grepl(paste0("^", resp_term, ":"),
                                                blk_ct$term_levels))
  ## Spillover columns in X_fixed
  fix_names <- rownames(fit$B)
  spill_idx <- which(grepl("_near$|spill", fix_names, ignore.case = TRUE))
  ## Unified spillover (2026-05-18+): phi bleed offset replaces _near covariate
  use_bleed <- !is.null(fit$technical_offset_mat) &&
               any(fit$technical_offset_mat != 0, na.rm = TRUE)
  ## Image RE columns (intercept-only), if any
  img_intercept_cols <- if (length(blk_idx_img)) {
    blk_img <- re$blocks[[blk_idx_img]]
    int_t <- which(blk_img$term_levels == "(Intercept)")
    if (!length(int_t)) integer(0)
    else blk_img$col_offset + (int_t - 1L) * blk_img$K_groups + seq_len(blk_img$K_groups)
  } else integer(0)

  rows_naka <- vector("list", G * length(focal_levels))
  ix2 <- 0L
  ct_means <- matrix(0, length(groups), G,
                     dimnames = list(groups, gene_names))
  for (c_idx in seq_along(groups)) {
    cells_c <- which(as.character(df$celltype) == groups[c_idx])
    if (length(cells_c) > 0)
      ct_means[c_idx, ] <- colMeans(Y[cells_c, , drop = FALSE], na.rm = TRUE)
  }
  for (c_name in focal_levels) {
    c_idx <- match(c_name, groups)
    cells_in_c <- which(as.character(df$celltype) == c_name)
    if (length(cells_in_c) < 5L) next
    n_c <- length(cells_in_c)
    N_c <- as.matrix(df[cells_in_c, vars, drop = FALSE])
    fmean      <- ct_means[c_idx, ]
    other_max  <- apply(ct_means[-c_idx, , drop = FALSE], 2, max, na.rm = TRUE)
    spec_focal <- fmean / pmax(fmean + other_max, 1e-9)
    foratio    <- fmean / pmax(other_max, 1e-9)

    ## Per-focal slope BLUPs for the no-Responder neighbour terms
    slope_term_idx <- vapply(vars, function(v) term2t[[v]], integer(1))

    for (gi in seq_len(G)) {
      g <- gene_names[gi]
      slope_blups <- vapply(slope_term_idx, function(ti)
        fit$U[ct_col(ti, c_idx), gi], numeric(1))
      slope_post  <- vapply(slope_term_idx, function(ti)
        fit$se_U[ct_col(ti, c_idx), gi]^2, numeric(1))

      ## Baseline spatial state: variance of (N_c %*% slope_blups) within focal
      eta_state <- as.numeric(N_c %*% slope_blups)
      var_N     <- apply(N_c, 2, stats::var, na.rm = TRUE)
      V_state_baseline <- stats::var(eta_state, na.rm = TRUE) +
                          sum(slope_post * var_N, na.rm = TRUE)

      ## Responder spatial state: variance of (R_c * N_c %*% resp_blups)
      ## within focal cells. R_c = 1 for PD cells, 0 for SD cells; the
      ## ResponderPD:neighbour interaction only fires in PD cells, so this
      ## block isolates the PD-specific spatial differential.
      V_state_responder <- 0
      if (has_resp) {
        ## FAIL LOUD: has_resp requires the condition indicator; absence used to
        ## silently leave V_state_responder = 0 (the .resp_dummy bug).
        if (!".resp_dummy" %in% colnames(df))
          stop("variance decomposition: resp_term given and interaction terms present, but df$.resp_dummy (the condition 0/1 indicator) is missing. The builder must set it.")
        resp_term_names <- paste0(resp_term, ":", vars)
        matched <- vapply(resp_term_names, function(nm) {
          v <- term2t[[nm]]; if (is.null(v)) NA_integer_ else as.integer(v)
        }, integer(1))
        keep_v   <- which(!is.na(matched))     # positions in `vars` that HAVE an interaction term
        resp_idx <- matched[keep_v]
        if (length(resp_idx)) {
          rblups <- vapply(resp_idx, function(ti)
            fit$U[ct_col(ti, c_idx), gi], numeric(1))
          rpost  <- vapply(resp_idx, function(ti)
            fit$se_U[ct_col(ti, c_idx), gi]^2, numeric(1))
          R_c <- df[[".resp_dummy"]][cells_in_c]   # exact match (avoid $ partial-match)
          N_r <- N_c[, vars[keep_v], drop = FALSE] * R_c   # the ACTUAL surviving neighbours (was vars[seq_along] — misaligned for subset interactions)
          eta_resp <- as.numeric(N_r %*% rblups)
          var_Nr   <- apply(N_r, 2, stats::var, na.rm = TRUE)
          V_state_responder <- stats::var(eta_resp, na.rm = TRUE) +
                                sum(rpost * var_Nr, na.rm = TRUE)
        }
      }

      ## Spillover: phi bleed offset (unified) else legacy _near fixed effects
      V_spill <- if (use_bleed) {
        stats::var(fit$technical_offset_mat[cells_in_c, gi], na.rm = TRUE)
      } else if (length(spill_idx)) {
        X_spill_c <- X_fixed[cells_in_c, spill_idx, drop = FALSE]
        beta_spill <- fit$B[spill_idx, gi]
        stats::var(as.numeric(X_spill_c %*% beta_spill), na.rm = TRUE)
      } else 0

      mu_c     <- fit$mu[cells_in_c, gi]
      mu_bar_c <- mean(mu_c, na.rm = TRUE)
      alpha_g  <- unname(fit$alpha[gi])
      ## Latent-scale residual variance for the NB GLM (Leckie et al. 2020
      ## Psych Methods 25:787-801). Default "nb1" matches PACE-MV's
      ## NB1 dispersion (Var(Y) = mu * (1 + alpha)). The legacy "nb2" form
      ## was Nakagawa et al. 2017's NB2 expression; it underestimates V_disp
      ## by approximately log((1+alpha)/(1+alpha/mu)) at small mu_bar, which
      ## inflates the higher-level percentages of the partition.
      V_disp <- if (disp_model == "nb1") {
        log(1 + (1 + max(alpha_g, 0)) / max(mu_bar_c, 1e-9))
      } else {
        log(1 + 1 / max(mu_bar_c, 1e-9) + max(alpha_g, 0))
      }

      ct_offset_sq <- fit$U[ct_col(term2t[["(Intercept)"]], c_idx), gi]^2 +
                       fit$se_U[ct_col(term2t[["(Intercept)"]], c_idx), gi]^2

      ## 5-block total: Cell type / Spatial state / Responder spatial state /
      ## Spillover / Residual. The image RE block stays in the model so it
      ## absorbs cross-patient variance away from the celltype slopes, but
      ## it is intentionally not surfaced in this decomposition (which
      ## matches the canonical PACE 5-block schema).
      tot_5 <- ct_offset_sq + V_state_baseline + V_state_responder +
               V_spill + V_disp
      if (!is.finite(tot_5) || tot_5 <= 0) next
      spec_g <- if (is.finite(spec_focal[gi])) spec_focal[gi] else 0

      rows_naka[[ix2 <- ix2 + 1L]] <- tibble::tibble(
        gene = g, focal = c_name, n_focal = n_c,
        V_state_baseline  = V_state_baseline,
        V_state_responder = V_state_responder,
        V_spill = V_spill, V_disp = V_disp,
        celltype_offset_sq = ct_offset_sq,
        focal_mean = fmean[gi], max_other_mean = other_max[gi],
        spec = spec_g, focal_other_ratio = foratio[gi],
        is_contaminated = !is.finite(foratio[gi]) || foratio[gi] < 1,
        `Cell type %`               = ct_offset_sq      / tot_5 * 100,
        `Spatial state %`           = V_state_baseline  / tot_5 * 100,
        `Responder spatial state %` = V_state_responder / tot_5 * 100,
        `Spillover %`               = V_spill           / tot_5 * 100,
        `Residual %`                = V_disp            / tot_5 * 100
      )
    }
  }
  gene_focal_5block <- dplyr::bind_rows(rows_naka[seq_len(ix2)])

  agg_focal_5block_mean <- gene_focal_5block |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      `Cell type %`               = mean(`Cell type %`,               na.rm = TRUE),
      `Spatial state %`           = mean(`Spatial state %`,           na.rm = TRUE),
      `Responder spatial state %` = mean(`Responder spatial state %`, na.rm = TRUE),
      `Spillover %`               = mean(`Spillover %`,               na.rm = TRUE),
      `Residual %`                = mean(`Residual %`,                na.rm = TRUE),
      n_genes = dplyr::n(), .groups = "drop"
    )

  ## Spec²-weighted aggregate (matches MCSD spec_exponent=2 weighting).
  ## Each gene's per-focal percentage is weighted by spec²; non-marker
  ## genes (low spec) are downweighted. Same denominator as MCSD. This is
  ## the manuscript-defensible aggregate when the panel contains many genes
  ## that are not specific to the focal celltype.
  agg_focal_5block_specw <- if (isTRUE(weight_by_spec_sq)) {
    gene_focal_5block |>
      dplyr::group_by(focal) |>
      dplyr::summarise(
        spec_w_sum = sum(spec^2, na.rm = TRUE),
        `Cell type %`               = sum(`Cell type %`               * spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        `Spatial state %`           = sum(`Spatial state %`           * spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        `Responder spatial state %` = sum(`Responder spatial state %` * spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        `Spillover %`               = sum(`Spillover %`               * spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        `Residual %`                = sum(`Residual %`                * spec^2, na.rm = TRUE) / pmax(spec_w_sum, 1e-12),
        n_genes = dplyr::n(),
        n_specific_genes = sum(!is_contaminated, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::select(-spec_w_sum)
  } else NULL

  ## Pooled (SS-weighted) aggregate per Goldstein-Browne-Rasbash 2002.
  ## Standard multilevel VPC literature reports the OVERALL share of
  ## explained variance: sum(numerator)/sum(denominator) across genes,
  ## NOT mean(numerator/denominator). The latter biases toward genes with
  ## small absolute total variance (where 1pp of % is meaningless). Pooled
  ## is the per-gene-magnitude-weighted statistic and is what actually
  ## quantifies "how much of the focal's outcome variance is in this block".
  gene_focal_5block <- gene_focal_5block |>
    dplyr::mutate(Total =
      V_state_baseline + V_state_responder + V_spill + V_disp + celltype_offset_sq)

  agg_focal_5block_pooled <- gene_focal_5block |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      total_SS = sum(Total, na.rm = TRUE),
      `Cell type %`               = 100 * sum(celltype_offset_sq, na.rm = TRUE) / pmax(sum(Total, na.rm = TRUE), 1e-12),
      `Spatial state %`           = 100 * sum(V_state_baseline,   na.rm = TRUE) / pmax(sum(Total, na.rm = TRUE), 1e-12),
      `Responder spatial state %` = 100 * sum(V_state_responder,  na.rm = TRUE) / pmax(sum(Total, na.rm = TRUE), 1e-12),
      `Spillover %`               = 100 * sum(V_spill,            na.rm = TRUE) / pmax(sum(Total, na.rm = TRUE), 1e-12),
      `Residual %`                = 100 * sum(V_disp,             na.rm = TRUE) / pmax(sum(Total, na.rm = TRUE), 1e-12),
      n_genes = dplyr::n(),
      .groups = "drop"
    )

  list(
    gene_focal_5block       = gene_focal_5block,
    agg_focal_5block_mean   = agg_focal_5block_mean,
    agg_focal_5block_specw  = agg_focal_5block_specw,
    agg_focal_5block_pooled = agg_focal_5block_pooled,
    disp_model              = disp_model
  )
}


# Leave-one-patient-out (LOPO) cluster-jackknife of the variance decomposition.
#
# Cameron, Gelbach & Miller (2008) Rev Econ Stat: at small cluster counts,
# cluster-robust standard errors and naive aggregates are biased. The
# jackknife provides bias-corrected point estimates and proper cluster-aware
# uncertainty bands.
#
# Use case: PACE-MV's per-focal variance decomposition averages across all
# panel genes and across all focal cells. At sparse focals (Mel B_Cell n=99,
# T_Cell n=278) a single patient's cells can dominate the per-focal % via
# their contribution to var(eta_resp) and the (image, celltype, term) BLUP
# magnitudes. Without LOPO, the headline % is fragile: e.g. Mel T_Cell
# Resp:Spatial = 33.6% drops to 20.1% if patient 32151 is excluded.
#
# Cheap LOPO (this function): holds the fitted BLUPs fixed, recomputes the
# per-cell variance partition on each LOPO subset of cells. Tests sensitivity
# of the AGGREGATION to individual patients, not the BLUP estimates.
#
# Proper LOPO (not implemented; cost = n_patients refits): refits the model
# without each patient. Tests sensitivity of the BLUPs themselves. Reserve for
# final manuscript figures if reviewers push.
#
# @param fit, df, Y, vars, X_fixed, resp_term, focal_levels, disp_model,
#   weight_by_spec_sq same as mvpql_variance_decomposition_multi
# @param patient_col column name in df identifying patients/clusters (e.g.
#   "imageID" for Mel/DKD, "patient_id" for CRC).
# @return list with:
#   - `summary` per-focal table with full / LOPO mean / SD / min / max / max-min
#   - `worst_offender` per-focal table identifying the patient driving the
#     largest single-leave-out drop in Resp:Spatial %
#   - `lopo_long` row per (LOPO patient, focal): the full agg_focal_5block_mean
#     for each LOPO subset, useful for additional diagnostics
#   - `full` the full-data agg_focal_5block_mean for reference
mvpql_variance_decomposition_lopo <- function(fit, df, Y, vars, X_fixed,
                                                patient_col,
                                                resp_term = NULL,
                                                focal_levels = NULL,
                                                disp_model = c("nb1", "nb2"),
                                                weight_by_spec_sq = FALSE,
                                                BPPARAM = NULL,
                                                verbose = TRUE) {
  if (!patient_col %in% colnames(df))
    stop(sprintf("patient_col '%s' not in df", patient_col))
  if (is.null(BPPARAM)) BPPARAM <- BiocParallel::SerialParam()

  full <- mvpql_variance_decomposition_multi(
    fit = fit, df = df, Y = Y, vars = vars, X_fixed = X_fixed,
    resp_term = resp_term, focal_levels = focal_levels,
    disp_model = disp_model,
    weight_by_spec_sq = weight_by_spec_sq)$agg_focal_5block_mean

  patients <- as.character(unique(df[[patient_col]]))
  if (verbose)
    message(sprintf("[LOPO] %d %s units; %d focals; cheap-mode (no refit); BPPARAM=%s",
                    length(patients), patient_col,
                    if (is.null(focal_levels))
                      length(intersect(unique(as.character(df$celltype)),
                                       fit$re_meta$blocks[[1]]$group_levels))
                    else length(focal_levels),
                    class(BPPARAM)[1L]))

  one_patient <- function(p) {
    keep <- as.character(df[[patient_col]]) != p
    df_p <- df[keep, , drop = FALSE]
    Y_p  <- Y[keep, , drop = FALSE]
    X_p  <- X_fixed[keep, , drop = FALSE]
    fit_p <- fit
    fit_p$mu <- fit_p$mu[keep, , drop = FALSE]
    ## The consumer indexes technical_offset_mat by positions in the SUBSET
    ## frame, so it has to be subset alongside mu; leaving it at full length
    ## silently offsets every cell after the dropped patient's first row.
    if (!is.null(fit_p$technical_offset_mat))
      fit_p$technical_offset_mat <-
        fit_p$technical_offset_mat[keep, , drop = FALSE]
    dec_p <- mvpql_variance_decomposition_multi(
      fit = fit_p, df = df_p, Y = Y_p, vars = vars, X_fixed = X_p,
      resp_term = resp_term, focal_levels = focal_levels,
      disp_model = disp_model,
      weight_by_spec_sq = weight_by_spec_sq)$agg_focal_5block_mean
    dec_p$drop_patient <- p
    dec_p
  }

  results <- BiocParallel::bplapply(patients, one_patient, BPPARAM = BPPARAM)

  lopo_long <- dplyr::bind_rows(results)

  summary <- lopo_long |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      n_lopo = dplyr::n(),
      mean   = mean(`Responder spatial state %`),
      sd     = stats::sd(`Responder spatial state %`),
      min    = min(`Responder spatial state %`),
      max    = max(`Responder spatial state %`),
      max_minus_min = max - min,
      .groups = "drop") |>
    dplyr::left_join(
      full |> dplyr::select(focal, full = `Responder spatial state %`),
      by = "focal") |>
    dplyr::arrange(dplyr::desc(full))

  worst_offender <- lopo_long |>
    dplyr::group_by(focal) |>
    dplyr::slice_min(`Responder spatial state %`, n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(focal, drop_patient,
                     resp_pct_drop = `Responder spatial state %`) |>
    dplyr::left_join(
      full |> dplyr::select(focal, full = `Responder spatial state %`),
      by = "focal") |>
    dplyr::mutate(delta = full - resp_pct_drop) |>
    dplyr::arrange(dplyr::desc(delta))

  list(summary = summary,
       worst_offender = worst_offender,
       lopo_long = lopo_long,
       full = full)
}


# FAST leave-one-patient-out variance decomposition via incremental sums.
#
# Same statistical output as `mvpql_variance_decomposition_lopo()` but ~10-30x
# faster. The key is that variance at a LOPO subset is closed-form from the
# full sums and per-patient subtotals:
#   var(x[!p]) = (Sxx - Sxx_p) / (n - n_p) - ((Sx - Sx_p) / (n - n_p))^2
# Compute Sx, Sxx ONCE per (focal, gene); Sx_p, Sxx_p once per (focal,
# gene, patient); per LOPO patient just subtract and divide.
#
# DKD: 40 min serial -> ~3 min via incremental.
# CRC: 15 min -> ~1 min.
# Mel: 5 min -> ~20 sec.
#
# Numerically equivalent to the slow version up to floating-point roundoff
# (~1e-10 typical max diff on the per-focal Resp:Spatial %).
#
# @param variance_estimator one of:
#   - "within_patient" (default 2026-05-10): one-way ANOVA decomposition
#     by patient; numerator = SSW = sum_p (Sxx_p - Sx_p²/n_p), df_W =
#     sum_p (n_p - 1) = n_total - K_patients. Robust to between-patient
#     drift, which can inflate the total-variance estimator at sparse
#     focals (Mel B_Cell n=99 drops from 26% to 5% under within-patient).
#     This is the manuscript-defensible choice when Mundlak hybrid +
#     patient random slopes are claimed to absorb between-patient signal.
#   - "total": classical sample variance over all focal cells (legacy).
#     Conflates within- and between-patient variation.
#
# Output structure identical to mvpql_variance_decomposition_lopo:
# list(summary, worst_offender, lopo_long, full).
mvpql_variance_decomposition_lopo_fast <- function(fit, df, Y, vars, X_fixed,
                                                     patient_col,
                                                     resp_term = NULL,
                                                     focal_levels = NULL,
                                                     disp_model = c("nb1", "nb2"),
                                                     variance_estimator = c("within_patient", "total"),
                                                     verbose = TRUE) {
  disp_model <- match.arg(disp_model)
  variance_estimator <- match.arg(variance_estimator)
  if (!patient_col %in% colnames(df))
    stop(sprintf("patient_col '%s' not in df", patient_col))

  re <- fit$re_meta
  blk_idx_ct <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "celltype")
  if (length(blk_idx_ct) != 1L) stop("Expected exactly one celltype RE block")
  blk_ct <- re$blocks[[blk_idx_ct]]
  groups <- blk_ct$group_levels
  if (is.null(focal_levels)) focal_levels <- intersect(as.character(unique(df$celltype)), groups)
  gene_names <- colnames(fit$U); G <- length(gene_names)
  ct_col <- function(t_idx, g_idx) blk_ct$col_offset + (t_idx - 1L) * blk_ct$K_groups + g_idx
  term2t <- setNames(seq_along(blk_ct$term_levels), blk_ct$term_levels)
  has_resp <- !is.null(resp_term) && any(grepl(paste0("^", resp_term, ":"), blk_ct$term_levels))
  fix_names <- rownames(fit$B)
  spill_idx <- which(grepl("_near$|spill", fix_names, ignore.case = TRUE))
  ## Unified spillover (2026-05-18+): phi bleed offset replaces _near covariate
  use_bleed <- !is.null(fit$technical_offset_mat) &&
               any(fit$technical_offset_mat != 0, na.rm = TRUE)

  ## Per-celltype mean for spec calculation
  ct_means <- matrix(0, length(groups), G, dimnames = list(groups, gene_names))
  for (c_idx in seq_along(groups)) {
    cells_c <- which(as.character(df$celltype) == groups[c_idx])
    if (length(cells_c) > 0)
      ct_means[c_idx, ] <- colMeans(Y[cells_c, , drop = FALSE], na.rm = TRUE)
  }

  ## Helper: incremental variance for a vector x given full sums and per-row
  ## subset to remove. Returns var of x[setdiff].
  ##   var(x[keep]) = (Sxx - Sxx_p) / (n - n_p) - mean²
  inc_var <- function(Sxx_full, Sx_full, n_full, Sxx_p, Sx_p, n_p) {
    n_keep <- n_full - n_p
    Sxx_k <- Sxx_full - Sxx_p
    Sx_k  <- Sx_full  - Sx_p
    mean_k <- Sx_k / n_keep
    pmax((Sxx_k - n_keep * mean_k^2) / pmax(n_keep - 1, 1), 0)
  }

  patient_vec <- as.character(df[[patient_col]])
  patients <- unique(patient_vec)
  n_pat <- length(patients)
  pat_idx <- match(patient_vec, patients)

  if (verbose)
    message(sprintf("[LOPO-fast] %d %s units; %d focals; %d genes; incremental-variance mode",
                    n_pat, patient_col, length(focal_levels), G))

  ## Build the per-focal output, one table per (gene, focal, [LOPO-patient or full])
  ## that can later be aggregated.
  lopo_long_list <- vector("list", length(focal_levels))
  full_list      <- vector("list", length(focal_levels))

  for (fi in seq_along(focal_levels)) {
    c_name <- focal_levels[fi]; c_idx <- match(c_name, groups)
    cells_in_c <- which(as.character(df$celltype) == c_name)
    if (length(cells_in_c) < 5L) next
    n_c <- length(cells_in_c)
    pat_c <- pat_idx[cells_in_c]
    N_c <- as.matrix(df[cells_in_c, vars, drop = FALSE])
    if (has_resp && !".resp_dummy" %in% colnames(df))
      stop("variance decomposition (lopo): resp_term given but df$.resp_dummy is missing.")
    R_c <- if (".resp_dummy" %in% colnames(df)) df[[".resp_dummy"]][cells_in_c] else rep(0, n_c)

    ## Precompute the per-(focal cell, gene) eta matrices.
    slope_term_idx <- vapply(vars, function(v) term2t[[v]], integer(1))
    slope_blups_g <- fit$U[ct_col(slope_term_idx, c_idx), , drop = FALSE]   # K_t × G
    slope_post_g  <- fit$se_U[ct_col(slope_term_idx, c_idx), , drop = FALSE]^2
    eta_state_mat <- N_c %*% slope_blups_g                                   # n_c × G

    if (has_resp) {
      resp_term_names <- paste0(resp_term, ":", vars)
      matched <- vapply(resp_term_names, function(nm) {
        v <- term2t[[nm]]; if (is.null(v)) NA_integer_ else as.integer(v)
      }, integer(1))
      keep_v   <- which(!is.na(matched))
      resp_idx <- matched[keep_v]
      rblups_g <- fit$U[ct_col(resp_idx, c_idx), , drop = FALSE]
      rpost_g  <- fit$se_U[ct_col(resp_idx, c_idx), , drop = FALSE]^2
      N_r_c    <- N_c[, vars[keep_v], drop = FALSE] * R_c   # surviving neighbours (was vars[seq_along] — misaligned)
      eta_resp_mat <- N_r_c %*% rblups_g                                      # n_c × G
    } else {
      eta_resp_mat <- matrix(0, n_c, G)
      rpost_g <- matrix(0, length(slope_term_idx), G)
      N_r_c   <- N_c * 0
    }

    eta_spill_mat <- if (use_bleed) {
      as.matrix(fit$technical_offset_mat[cells_in_c, , drop = FALSE])             # n_c × G
    } else if (length(spill_idx)) {
      X_spill_c   <- X_fixed[cells_in_c, spill_idx, drop = FALSE]
      beta_spill  <- fit$B[spill_idx, , drop = FALSE]
      X_spill_c %*% beta_spill                                                # n_c × G
    } else matrix(0, n_c, G)

    mu_c_mat <- as.matrix(fit$mu[cells_in_c, , drop = FALSE])                 # n_c × G

    ## Spec / contam / cell-type-offset (constant across LOPO subsets)
    fmean      <- ct_means[c_idx, ]
    other_max  <- apply(ct_means[-c_idx, , drop = FALSE], 2, max, na.rm = TRUE)
    other_max[!is.finite(other_max)] <- 0
    spec_focal <- fmean / pmax(fmean + other_max, 1e-9)
    foratio    <- fmean / pmax(other_max, 1e-9)
    is_contam  <- !is.finite(foratio) | foratio < 1
    ct_offset_sq <- fit$U[ct_col(term2t[["(Intercept)"]], c_idx), ]^2 +
                     fit$se_U[ct_col(term2t[["(Intercept)"]], c_idx), ]^2

    alpha_g <- pmax(unname(fit$alpha), 0)

    ## ---- Full-data per-(focal, gene) sums ----
    Sx_eta_state  <- colSums(eta_state_mat)
    Sxx_eta_state <- colSums(eta_state_mat^2)
    Sx_eta_resp   <- colSums(eta_resp_mat)
    Sxx_eta_resp  <- colSums(eta_resp_mat^2)
    Sx_eta_spill  <- colSums(eta_spill_mat)
    Sxx_eta_spill <- colSums(eta_spill_mat^2)
    Sx_mu         <- colSums(mu_c_mat)
    ## For per-(t)-vars sums of N_c: K_t × 2
    Sx_N   <- colSums(N_c)         # K_t
    Sxx_N  <- colSums(N_c^2)       # K_t
    Sx_Nr  <- colSums(N_r_c)
    Sxx_Nr <- colSums(N_r_c^2)
    Sx_post_var_term  <- function(Sxx_N_, Sx_N_, n_, post_g) {
      ## var per term, weighted by post: sum_t post_g[t,g] * var_t
      var_per_term <- pmax((Sxx_N_ - n_ * (Sx_N_ / n_)^2) / pmax(n_ - 1, 1), 0)
      as.numeric(crossprod(post_g, var_per_term))
    }

    ## ---- Per-patient subtotals (n_pat × G), once per focal ----
    Sx_p_eta_state  <- matrix(0, n_pat, G)
    Sxx_p_eta_state <- matrix(0, n_pat, G)
    Sx_p_eta_resp   <- matrix(0, n_pat, G)
    Sxx_p_eta_resp  <- matrix(0, n_pat, G)
    Sx_p_eta_spill  <- matrix(0, n_pat, G)
    Sxx_p_eta_spill <- matrix(0, n_pat, G)
    Sx_p_mu         <- matrix(0, n_pat, G)
    n_p_focal       <- integer(n_pat)
    Sx_p_N   <- matrix(0, n_pat, length(slope_term_idx))
    Sxx_p_N  <- matrix(0, n_pat, length(slope_term_idx))
    Sx_p_Nr  <- matrix(0, n_pat, length(slope_term_idx))
    Sxx_p_Nr <- matrix(0, n_pat, length(slope_term_idx))
    for (pi_local in seq_len(n_pat)) {
      idx_p <- which(pat_c == pi_local)
      n_p_focal[pi_local] <- length(idx_p)
      if (!length(idx_p)) next
      Sx_p_eta_state[pi_local, ]  <- colSums(eta_state_mat[idx_p, , drop = FALSE])
      Sxx_p_eta_state[pi_local, ] <- colSums(eta_state_mat[idx_p, , drop = FALSE]^2)
      Sx_p_eta_resp[pi_local, ]   <- colSums(eta_resp_mat[idx_p, , drop = FALSE])
      Sxx_p_eta_resp[pi_local, ]  <- colSums(eta_resp_mat[idx_p, , drop = FALSE]^2)
      Sx_p_eta_spill[pi_local, ]  <- colSums(eta_spill_mat[idx_p, , drop = FALSE])
      Sxx_p_eta_spill[pi_local, ] <- colSums(eta_spill_mat[idx_p, , drop = FALSE]^2)
      Sx_p_mu[pi_local, ]         <- colSums(mu_c_mat[idx_p, , drop = FALSE])
      Sx_p_N[pi_local, ]   <- colSums(N_c[idx_p, , drop = FALSE])
      Sxx_p_N[pi_local, ]  <- colSums(N_c[idx_p, , drop = FALSE]^2)
      Sx_p_Nr[pi_local, ]  <- colSums(N_r_c[idx_p, , drop = FALSE])
      Sxx_p_Nr[pi_local, ] <- colSums(N_r_c[idx_p, , drop = FALSE]^2)
    }

    ## ---- Within-patient SS contributions (only patients with n_p >= 2) ----
    ## For variance_estimator = "within_patient": SSW = sum_p (Sxx_p - Sx_p^2/n_p)
    ## Only patients with n_p >= 2 contribute (need >=2 cells to estimate var).
    qq <- which(n_p_focal >= 2L)
    if (length(qq) > 0L) {
      np_q <- n_p_focal[qq]
      ssw_eta_state_q <- pmax(Sxx_p_eta_state[qq, , drop=FALSE] -
                               Sx_p_eta_state[qq, , drop=FALSE]^2 / np_q, 0)
      ssw_eta_resp_q  <- pmax(Sxx_p_eta_resp[qq, , drop=FALSE]  -
                               Sx_p_eta_resp[qq, , drop=FALSE]^2  / np_q, 0)
      ssw_eta_spill_q <- pmax(Sxx_p_eta_spill[qq, , drop=FALSE] -
                               Sx_p_eta_spill[qq, , drop=FALSE]^2 / np_q, 0)
      ssw_N_q  <- pmax(Sxx_p_N[qq, , drop=FALSE]  - Sx_p_N[qq, , drop=FALSE]^2  / np_q, 0)
      ssw_Nr_q <- pmax(Sxx_p_Nr[qq, , drop=FALSE] - Sx_p_Nr[qq, , drop=FALSE]^2 / np_q, 0)
      ssw_eta_state_full <- colSums(ssw_eta_state_q)
      ssw_eta_resp_full  <- colSums(ssw_eta_resp_q)
      ssw_eta_spill_full <- colSums(ssw_eta_spill_q)
      ssw_N_full  <- colSums(ssw_N_q)
      ssw_Nr_full <- colSums(ssw_Nr_q)
      df_w_full <- sum(np_q - 1L)
    } else {
      ssw_eta_state_full <- ssw_eta_resp_full <- ssw_eta_spill_full <- numeric(G)
      ssw_N_full <- ssw_Nr_full <- numeric(length(slope_term_idx))
      df_w_full <- 0L
    }

    ## Full decomp — pick total or within depending on variance_estimator
    if (variance_estimator == "total") {
      var_eta_state_full <- inc_var(Sxx_eta_state, Sx_eta_state, n_c, 0, 0, 0)
      var_eta_resp_full  <- inc_var(Sxx_eta_resp,  Sx_eta_resp,  n_c, 0, 0, 0)
      var_eta_spill_full <- inc_var(Sxx_eta_spill, Sx_eta_spill, n_c, 0, 0, 0)
      sum_post_var_N_full   <- Sx_post_var_term(Sxx_N,  Sx_N,  n_c, slope_post_g)
      sum_rpost_var_Nr_full <- Sx_post_var_term(Sxx_Nr, Sx_Nr, n_c, rpost_g)
    } else {  ## "within_patient"
      d <- max(df_w_full, 1L)
      var_eta_state_full <- ssw_eta_state_full / d
      var_eta_resp_full  <- ssw_eta_resp_full  / d
      var_eta_spill_full <- ssw_eta_spill_full / d
      var_per_term_N_full  <- ssw_N_full  / d
      var_per_term_Nr_full <- ssw_Nr_full / d
      sum_post_var_N_full   <- as.numeric(crossprod(slope_post_g, var_per_term_N_full))
      sum_rpost_var_Nr_full <- as.numeric(crossprod(rpost_g,      var_per_term_Nr_full))
    }
    mu_bar_full <- Sx_mu / n_c
    V_state_baseline_full  <- var_eta_state_full + sum_post_var_N_full
    V_state_responder_full <- var_eta_resp_full  + sum_rpost_var_Nr_full
    V_spill_full           <- var_eta_spill_full
    V_disp_full <- if (disp_model == "nb1") {
      log(1 + (1 + alpha_g) / pmax(mu_bar_full, 1e-9))
    } else {
      log(1 + 1 / pmax(mu_bar_full, 1e-9) + alpha_g)
    }
    tot_full <- ct_offset_sq + V_state_baseline_full + V_state_responder_full +
                V_spill_full + V_disp_full

    full_focal <- tibble::tibble(
      gene = gene_names, focal = c_name, n_focal = n_c,
      `Cell type %`               = 100 * ct_offset_sq            / tot_full,
      `Spatial state %`           = 100 * V_state_baseline_full   / tot_full,
      `Responder spatial state %` = 100 * V_state_responder_full  / tot_full,
      `Spillover %`               = 100 * V_spill_full            / tot_full,
      `Residual %`                = 100 * V_disp_full             / tot_full,
      spec = spec_focal, is_contaminated = is_contam)
    full_list[[fi]] <- full_focal

    ## ---- Per-LOPO incremental decomp ----
    lopo_focal <- vector("list", n_pat)
    for (pi_local in seq_len(n_pat)) {
      n_p <- n_p_focal[pi_local]
      n_keep <- n_c - n_p
      if (n_keep < 5L) next

      if (variance_estimator == "total") {
        var_es <- inc_var(Sxx_eta_state, Sx_eta_state, n_c,
                           Sxx_p_eta_state[pi_local, ], Sx_p_eta_state[pi_local, ], n_p)
        var_er <- inc_var(Sxx_eta_resp,  Sx_eta_resp,  n_c,
                           Sxx_p_eta_resp[pi_local, ],  Sx_p_eta_resp[pi_local, ],  n_p)
        var_sp <- inc_var(Sxx_eta_spill, Sx_eta_spill, n_c,
                           Sxx_p_eta_spill[pi_local, ], Sx_p_eta_spill[pi_local, ], n_p)
        var_per_term_N <- pmax(((Sxx_N - Sxx_p_N[pi_local, ]) -
                                (Sx_N - Sx_p_N[pi_local, ])^2 / n_keep) / pmax(n_keep - 1, 1),
                               0)
        var_per_term_Nr <- pmax(((Sxx_Nr - Sxx_p_Nr[pi_local, ]) -
                                  (Sx_Nr - Sx_p_Nr[pi_local, ])^2 / n_keep) / pmax(n_keep - 1, 1),
                                0)
      } else {  ## within_patient: SSW_LOPO = SSW_full - SSW_p, df_W -= (n_p - 1) if n_p >= 2
        if (n_p >= 2L) {
          q_pos <- match(pi_local, qq)
          ssw_es_drop <- ssw_eta_state_full - ssw_eta_state_q[q_pos, ]
          ssw_er_drop <- ssw_eta_resp_full  - ssw_eta_resp_q[q_pos, ]
          ssw_sp_drop <- ssw_eta_spill_full - ssw_eta_spill_q[q_pos, ]
          ssw_N_drop  <- ssw_N_full  - ssw_N_q[q_pos, ]
          ssw_Nr_drop <- ssw_Nr_full - ssw_Nr_q[q_pos, ]
          d_drop      <- max(df_w_full - (n_p - 1L), 1L)
        } else {
          ## n_p < 2: patient contributed nothing to SSW; LOPO == full
          ssw_es_drop <- ssw_eta_state_full
          ssw_er_drop <- ssw_eta_resp_full
          ssw_sp_drop <- ssw_eta_spill_full
          ssw_N_drop  <- ssw_N_full
          ssw_Nr_drop <- ssw_Nr_full
          d_drop      <- max(df_w_full, 1L)
        }
        var_es <- ssw_es_drop / d_drop
        var_er <- ssw_er_drop / d_drop
        var_sp <- ssw_sp_drop / d_drop
        var_per_term_N  <- ssw_N_drop  / d_drop
        var_per_term_Nr <- ssw_Nr_drop / d_drop
      }
      sum_post_var_N_p   <- as.numeric(crossprod(slope_post_g, var_per_term_N))
      sum_rpost_var_Nr_p <- as.numeric(crossprod(rpost_g,      var_per_term_Nr))
      mu_bar_p <- (Sx_mu - Sx_p_mu[pi_local, ]) / n_keep
      V_state_baseline_p  <- var_es + sum_post_var_N_p
      V_state_responder_p <- var_er + sum_rpost_var_Nr_p
      V_spill_p           <- var_sp
      V_disp_p <- if (disp_model == "nb1") {
        log(1 + (1 + alpha_g) / pmax(mu_bar_p, 1e-9))
      } else {
        log(1 + 1 / pmax(mu_bar_p, 1e-9) + alpha_g)
      }
      tot_p <- ct_offset_sq + V_state_baseline_p + V_state_responder_p +
               V_spill_p + V_disp_p
      lopo_focal[[pi_local]] <- tibble::tibble(
        gene = gene_names, focal = c_name, n_focal = n_keep,
        drop_patient = patients[pi_local],
        `Cell type %`               = 100 * ct_offset_sq          / tot_p,
        `Spatial state %`           = 100 * V_state_baseline_p    / tot_p,
        `Responder spatial state %` = 100 * V_state_responder_p   / tot_p,
        `Spillover %`               = 100 * V_spill_p             / tot_p,
        `Residual %`                = 100 * V_disp_p              / tot_p,
        spec = spec_focal, is_contaminated = is_contam)
    }
    lopo_long_list[[fi]] <- dplyr::bind_rows(lopo_focal)

    if (verbose) message(sprintf("[LOPO-fast] focal %d/%d (%s) done",
                                  fi, length(focal_levels), c_name))
  }

  full      <- dplyr::bind_rows(full_list)
  lopo_long <- dplyr::bind_rows(lopo_long_list)

  ## Per-(focal, drop_patient) aggregate (mean across genes)
  lopo_agg <- lopo_long |>
    dplyr::group_by(focal, drop_patient) |>
    dplyr::summarise(
      `Cell type %`               = mean(`Cell type %`,               na.rm = TRUE),
      `Spatial state %`           = mean(`Spatial state %`,           na.rm = TRUE),
      `Responder spatial state %` = mean(`Responder spatial state %`, na.rm = TRUE),
      `Spillover %`               = mean(`Spillover %`,               na.rm = TRUE),
      `Residual %`                = mean(`Residual %`,                na.rm = TRUE),
      .groups = "drop")

  full_agg <- full |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      `Cell type %`               = mean(`Cell type %`,               na.rm = TRUE),
      `Spatial state %`           = mean(`Spatial state %`,           na.rm = TRUE),
      `Responder spatial state %` = mean(`Responder spatial state %`, na.rm = TRUE),
      `Spillover %`               = mean(`Spillover %`,               na.rm = TRUE),
      `Residual %`                = mean(`Residual %`,                na.rm = TRUE),
      .groups = "drop")

  summary <- lopo_agg |>
    dplyr::group_by(focal) |>
    dplyr::summarise(
      n_lopo = dplyr::n(),
      mean   = mean(`Responder spatial state %`),
      sd     = stats::sd(`Responder spatial state %`),
      min    = min(`Responder spatial state %`),
      max    = max(`Responder spatial state %`),
      max_minus_min = max - min,
      .groups = "drop") |>
    dplyr::left_join(
      full_agg |> dplyr::select(focal, full = `Responder spatial state %`),
      by = "focal") |>
    dplyr::arrange(dplyr::desc(full))

  worst_offender <- lopo_agg |>
    dplyr::group_by(focal) |>
    dplyr::slice_min(`Responder spatial state %`, n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(focal, drop_patient,
                     resp_pct_drop = `Responder spatial state %`) |>
    dplyr::left_join(
      full_agg |> dplyr::select(focal, full = `Responder spatial state %`),
      by = "focal") |>
    dplyr::mutate(delta = full - resp_pct_drop) |>
    dplyr::arrange(dplyr::desc(delta))

  list(summary = summary,
       worst_offender = worst_offender,
       lopo_long = lopo_agg,
       full = full_agg,
       gene_lopo_long = lopo_long,
       gene_full = full)
}


# Convert mv-PQL output to canonical PACE ran_vals tibble (per-gene)
#
# Same shape as `extractRandomEffects()` returns -- one row per
# (gene, group, level, term) with estimate, std.error, scaled_estimate, pval,
# lower, upper. Drop-in for downstream `helpers/mashr_pipeline.R` etc.
mvpql_to_ran_vals <- function(fit, group_col = "celltype") {
  re <- fit$re_meta
  q  <- ncol(re$Z); g_n <- ncol(fit$U)
  gene_names <- colnames(fit$U)
  if (is.null(gene_names)) gene_names <- paste0("g", seq_len(g_n))

  out <- vector("list", g_n)
  for (gi in seq_len(g_n)) {
    est <- fit$U[, gi]
    se  <- fit$se_U[, gi]
    sca <- est / se
    out[[gi]] <- tibble::tibble(
      effect    = "ran_vals",
      component = "cond",
      group     = group_col,
      level     = re$group_of_col,
      term      = re$term_of_col,
      estimate  = est,
      std.error = se,
      gene      = gene_names[gi],
      lower     = est - 2 * se,
      upper     = est + 2 * se,
      scaled_estimate = sca,
      pval      = 2 * stats::pnorm(-abs(sca))
    )
  }
  do.call(rbind, out)
}


# Per-(focal, neighbour) pair-score per the manuscript Eq. (pair-score).
#
# For each (focal cell type c, neighbour cell type t) pair and each block in
# \{spatial, responder_spatial\}, computes
# \deqn{\mathrm{pair\!-\!score}_{c,t} = \overline{(\tilde{b}_{j,c,t}^{2}
#   + \tilde{s}_{j,c,t}^{2})}_j \cdot \widehat{\mathrm{Var}}_{i \in c}(N^{(c)}_{i,t})}
# where \eqn{\tilde{b}, \tilde{s}} are mashr-shrunken slope and posterior SD
# (from `apply_mashr_shrinkage`), and `Var(N_{i,t})` is the variance of the
# neighbour-t count over cells of focal celltype c, optionally computed under
# the within-patient ANOVA decomposition.
#
# Cluster LOPO: leave each patient out, recompute the variance factor (the
# shrunken posterior is held fixed -- refitting mashr per LOPO patient is
# too expensive). Isolates design-side cluster leverage in the pair score.
mvpql_pair_score <- function(shrunken_long = NULL, fit = NULL, df, vars,
                              patient_col, resp_term,
                              focal_levels = NULL,
                              variance_estimator = c("within_patient", "total"),
                              posterior_source = c("raw_pql", "mashr"),
                              include_lopo = TRUE,
                              verbose = TRUE) {
  variance_estimator <- match.arg(variance_estimator)
  posterior_source <- match.arg(posterior_source)
  if (!patient_col %in% colnames(df))
    stop(sprintf("patient_col '%s' not in df", patient_col))
  if (is.null(focal_levels)) focal_levels <- levels(df$celltype)
  if (is.null(focal_levels)) focal_levels <- sort(unique(as.character(df$celltype)))

  ## Posterior source dispatch:
  ## - "raw_pql": pull (u, se_U) from fit$U / fit$se_U keyed by "<focal>::<term>".
  ## - "mashr"  : pull (estimate_shrunk, sd_shrunk) from shrunken_long.
  if (posterior_source == "raw_pql") {
    if (is.null(fit)) stop("posterior_source='raw_pql' requires fit argument")
    rn <- rownames(fit$U)
    if (is.null(rn)) stop("fit$U has no rownames; cannot key by focal::term")
    gene_names_fit <- colnames(fit$U)
  } else {
    if (is.null(shrunken_long))
      stop("posterior_source='mashr' requires shrunken_long argument")
    shr <- shrunken_long |>
      dplyr::filter(!is.na(estimate_shrunk), !is.na(sd_shrunk))
  }

  patient_vec <- as.character(df[[patient_col]])
  patients <- unique(patient_vec)
  pat_idx <- match(patient_vec, patients)

  rows <- list()
  for (c_name in focal_levels) {
    cells_in_c <- which(as.character(df$celltype) == c_name)
    n_c <- length(cells_in_c)
    if (n_c < 5L) next
    pat_c <- pat_idx[cells_in_c]
    pats_in_c <- sort(unique(pat_c))
    np <- length(pats_in_c)

    N_c <- as.matrix(df[cells_in_c, vars, drop = FALSE])
    Sx_p_N  <- matrix(0, np, length(vars))
    Sxx_p_N <- matrix(0, np, length(vars))
    n_per_p <- integer(np)
    for (i in seq_along(pats_in_c)) {
      idx_p <- which(pat_c == pats_in_c[i])
      n_per_p[i] <- length(idx_p)
      if (!length(idx_p)) next
      Sx_p_N[i, ]  <- colSums(N_c[idx_p, , drop = FALSE])
      Sxx_p_N[i, ] <- colSums(N_c[idx_p, , drop = FALSE]^2)
    }
    Sx_N_full  <- colSums(N_c)
    Sxx_N_full <- colSums(N_c^2)

    qq <- which(n_per_p >= 2L)
    if (length(qq) > 0L) {
      ssw_per_p_q <- pmax(Sxx_p_N[qq, , drop = FALSE] -
                           Sx_p_N[qq, , drop = FALSE]^2 / n_per_p[qq], 0)
      ssw_full <- colSums(ssw_per_p_q)
      df_w_full <- sum(n_per_p[qq] - 1L)
    } else {
      ssw_full <- numeric(length(vars)); df_w_full <- 0L
    }

    var_factor_full <- if (variance_estimator == "within_patient") {
      if (df_w_full < 1L) numeric(length(vars)) else ssw_full / df_w_full
    } else {
      mean_N <- Sx_N_full / n_c
      pmax((Sxx_N_full - n_c * mean_N^2) / pmax(n_c - 1L, 1L), 0)
    }

    for (block in c("spatial", "responder_spatial")) {
      for (ti in seq_along(vars)) {
        t_name <- vars[ti]
        target_term <- if (block == "spatial") t_name
                        else paste0(resp_term, ":", t_name)
        if (posterior_source == "raw_pql") {
          row_key <- paste0(c_name, "::", target_term)
          ix <- match(row_key, rn)
          if (is.na(ix)) next  ## term not present for this focal (e.g. spatial block when fit has only Responder×TYPES random slopes)
          u_vec    <- fit$U[ix, ]
          seU_vec  <- fit$se_U[ix, ]
          ok <- is.finite(u_vec) & is.finite(seU_vec)
          if (!any(ok)) next
          num_per_gene <- u_vec[ok]^2 + seU_vec[ok]^2
          n_genes_pair_local <- sum(ok)
        } else {
          rows_jt <- shr |> dplyr::filter(focal == c_name, term == target_term)
          if (nrow(rows_jt) == 0L) next
          num_per_gene <- rows_jt$estimate_shrunk^2 + rows_jt$sd_shrunk^2
          n_genes_pair_local <- nrow(rows_jt)
        }
        num_mean <- mean(num_per_gene, na.rm = TRUE)

        if (isTRUE(include_lopo) && np > 1L) {
          lopo_scores <- numeric(np)
          for (pi_local in seq_len(np)) {
            n_p <- n_per_p[pi_local]
            n_keep <- n_c - n_p
            if (n_keep < 5L) { lopo_scores[pi_local] <- NA_real_; next }
            if (variance_estimator == "within_patient") {
              if (n_p >= 2L) {
                q_pos <- match(pi_local, qq)
                ssw_drop <- ssw_full[ti] - ssw_per_p_q[q_pos, ti]
                d_drop <- max(df_w_full - (n_p - 1L), 1L)
                vf <- ssw_drop / d_drop
              } else {
                vf <- if (df_w_full >= 1L) ssw_full[ti] / df_w_full else 0
              }
            } else {
              Sxx_drop <- Sxx_N_full[ti] - Sxx_p_N[pi_local, ti]
              Sx_drop  <- Sx_N_full[ti]  - Sx_p_N[pi_local, ti]
              mean_drop <- Sx_drop / n_keep
              vf <- pmax((Sxx_drop - n_keep * mean_drop^2) / pmax(n_keep - 1L, 1L), 0)
            }
            lopo_scores[pi_local] <- num_mean * vf
          }
          lopo_mean <- mean(lopo_scores, na.rm = TRUE)
          lopo_sd   <- stats::sd(lopo_scores, na.rm = TRUE)
          have_lopo <- any(is.finite(lopo_scores))
          lopo_min  <- if (have_lopo) min(lopo_scores, na.rm = TRUE) else NA_real_
          lopo_max  <- if (have_lopo) max(lopo_scores, na.rm = TRUE) else NA_real_
          worst_pi  <- which.min(lopo_scores)
          worst_pat <- if (length(worst_pi)) as.character(patients[pats_in_c[worst_pi]]) else NA_character_
          worst_delta <- (num_mean * var_factor_full[ti]) - lopo_min
        } else {
          lopo_mean <- lopo_sd <- lopo_min <- lopo_max <- worst_delta <- NA_real_
          worst_pat <- NA_character_
        }

        ## Cache scalars before the tibble{} call -- tibble evaluates columns
        ## lazily, so `var_factor_full = var_factor_full[ti]` would shadow the
        ## outer vector and the next line would index a 1-element column.
        vf_ti <- var_factor_full[ti]
        ps_ti <- num_mean * vf_ti
        rows[[length(rows) + 1L]] <- tibble::tibble(
          focal = c_name, neighbour = t_name, block = block,
          n_genes_pair = n_genes_pair_local,
          num_mean = num_mean,
          var_factor_full = vf_ti,
          pair_score_full = ps_ti,
          lopo_mean = lopo_mean, lopo_sd = lopo_sd,
          lopo_min  = lopo_min,  lopo_max = lopo_max,
          worst_patient = worst_pat, worst_delta = worst_delta)
      }
    }
    if (verbose) message(sprintf("[pair_score] focal '%s' done (%d cells)", c_name, n_c))
  }
  dplyr::bind_rows(rows)
}


# Manuscript-spec 6-component variance decomposition (Methods line 38).
#
# Components per Methods text:
#   1. Cell type identity     = (u_(c,Int))^2 + V_post
#   2. Condition main effect  = beta_R^2 * Var(R over focal c)
#   3. Spatial cell state     = sum_t (u_(c,t)^2 + V_post) * Var(N_t over focal c)
#   4. Condition x spatial    = sum_t (u_(c,R:t)^2 + V_post) * Var(N_t * R over focal c)
#   5. imageID block          = Var across focal-c cells of (Z_img u_img)
#   6. Residual dispersion    = log(1 + (1+alpha)/mu_bar)        (NB1)
#
# Spillover-near is *not* in the Total (it's a technical correction per
# Methods line 16, not a biological variance component).
#
# Per-pair % for (focal c, neighbour t, block in {spatial, responder_spatial}):
#   pair-pct = 100 * pair_numerator / Total_g, mean across genes.
# where pair_numerator(c,t) = (u^2 + V_post) * Var(N_t over focal c) for that
# specific block.
mvpql_variance_decomposition_v6 <- function(fit, df, vars, X_fixed,
                                              resp_term = NULL,
                                              focal_levels = NULL,
                                              patient_col = NULL,
                                              variance_estimator = c("total","within_patient"),
                                              numerator_shrinkage = c("none","marshall","signal_share"),
                                              pair_allocation = c("lmg","diagonal"),
                                              ct_means = NULL,
                                              weight_by_spec_sq = FALSE,
                                              aggregate = c("pooled","median","mean"),
                                              verbose = TRUE) {
  variance_estimator <- match.arg(variance_estimator)
  numerator_shrinkage <- match.arg(numerator_shrinkage)
  pair_allocation <- match.arg(pair_allocation)
  aggregate <- match.arg(aggregate)
  if (isTRUE(weight_by_spec_sq) && is.null(ct_means))
    stop("weight_by_spec_sq=TRUE requires ct_means (genes x celltypes matrix)")
  ## Weighted mean and weighted median (Hoaglin-Mosteller-Tukey 1983)
  wmean <- function(x, w) {
    keep <- is.finite(x) & is.finite(w) & w > 0
    if (!any(keep)) return(NA_real_)
    sum(x[keep] * w[keep]) / sum(w[keep])
  }
  wmedian <- function(x, w) {
    keep <- is.finite(x) & is.finite(w) & w > 0
    if (!any(keep)) return(NA_real_)
    o <- order(x[keep]); xo <- x[keep][o]; wo <- w[keep][o]
    cw <- cumsum(wo) / sum(wo)
    xo[which(cw >= 0.5)[1L]]
  }
  agg_unweighted <- function(x) {
    if (aggregate == "mean") mean(x, na.rm = TRUE) else stats::median(x, na.rm = TRUE)
  }
  agg_weighted <- function(x, w) {
    if (aggregate == "mean") wmean(x, w) else wmedian(x, w)
  }
  ## Pooled (sum/sum) aggregation works on the un-percentaged numerators and
  ## denominators directly: pct = 100 * sum_g pair_var(g) / sum_g Total(g).
  ## That's the canonical variance partition coefficient (Goldstein, Browne &
  ## Rasbash 2002 Understanding Statistics 1:223-231; Snijders & Bosker 2012
  ## Multilevel Analysis ch. 7), variance-conserving across genes.
  agg_pooled <- function(num_vec, tot_vec) {
    keep <- is.finite(num_vec) & is.finite(tot_vec) & tot_vec > 0
    if (!any(keep)) return(NA_real_)
    100 * sum(num_vec[keep]) / sum(tot_vec[keep])
  }
  ## Pair allocation rule:
  ##   diagonal : pair_var(c,t,g) = (u^2 + V_post) * Var(N_t over focal c)
  ##              (assumes neighbour counts are uncorrelated across pairs)
  ##   lmg      : pair_var(c,t,g) = row_sum(cov(contrib_mat))[t] + V_post*Var(N_t)
  ##              where contrib_mat[i,t] = N_c[i,t] * u_t. Captures cross-pair
  ##              covariance via Lindeman-Merenda-Gold / Shapley split (Gromping
  ##              2007 Am Stat 61:139-147; Snijders & Bosker 2012 Multilevel
  ##              Analysis ch.7). Total Var(eta_state) is conserved.
  ## Shrinkage transform: maps (u^2, V_post) -> shrunken numerator
  ##   none           : u^2 + V_post (raw fitted variance)
  ##   marshall       : pmax(u^2 - V_post, 0) (Marshall 2003 JRSS-C)
  ##   signal_share   : u^2 * u^2/(u^2 + V_post) = u^4/(u^2 + V_post)
  shrink_uv <- function(u2, vpost) {
    switch(numerator_shrinkage,
      none         = u2 + vpost,
      marshall     = pmax(u2 - vpost, 0),
      signal_share = u2 * u2 / pmax(u2 + vpost, 1e-12))
  }
  if (is.null(patient_col) || !patient_col %in% colnames(df))
    stop("patient_col must name a column of df")

  re <- fit$re_meta
  blk_idx_ct  <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "celltype")
  blk_idx_img <- which(vapply(re$blocks, `[[`, character(1), "group_col") == "imageID")
  if (length(blk_idx_ct) != 1L) stop("expected one celltype RE block")
  blk_ct <- re$blocks[[blk_idx_ct]]
  blk_img <- if (length(blk_idx_img) == 1L) re$blocks[[blk_idx_img]] else NULL
  groups <- blk_ct$group_levels
  if (is.null(focal_levels))
    focal_levels <- intersect(unique(as.character(df$celltype)), groups)
  ct_col_fn <- function(t_idx, g_idx) blk_ct$col_offset + (t_idx - 1L) * blk_ct$K_groups + g_idx
  term2t_ct <- setNames(seq_along(blk_ct$term_levels), blk_ct$term_levels)
  gene_names <- colnames(fit$U); G <- length(gene_names)
  fix_names  <- rownames(fit$B)
  alpha_g <- pmax(unname(fit$alpha), 0)

  ## Resp main-effect column in X_fixed
  resp_fixed_col <- if (!is.null(resp_term))
    match(paste0(resp_term, "PD"), fix_names, nomatch = match(resp_term, fix_names, nomatch = NA))
    else NA_integer_
  ## Try common name variants
  if (is.na(resp_fixed_col) && !is.null(resp_term)) {
    cand <- grep(paste0("^", resp_term, "(PD|mut|DKD)?$"), fix_names)
    if (length(cand) == 1L) resp_fixed_col <- cand
  }

  ## Helper: within-patient or total variance of a numeric vector across rows
  var_within_or_total <- function(x_vec, pat_vec, est) {
    if (length(x_vec) < 5L) return(0)
    if (est == "total")
      return(stats::var(x_vec, na.rm = TRUE))
    pats <- unique(pat_vec)
    ssw <- 0; df_w <- 0
    for (p in pats) {
      idx <- which(pat_vec == p)
      if (length(idx) < 2L) next
      m <- mean(x_vec[idx], na.rm = TRUE)
      ssw <- ssw + sum((x_vec[idx] - m)^2, na.rm = TRUE)
      df_w <- df_w + length(idx) - 1L
    }
    if (df_w < 1L) 0 else ssw / df_w
  }
  ## same for matrix-by-column (returns per-column var vector)
  var_within_mat <- function(M, pat_vec, est) {
    apply(M, 2, var_within_or_total, pat_vec = pat_vec, est = est)
  }

  pair_rows <- list(); focal_rows <- list()
  Z_full <- re$Z

  ## Pre-compute spec_focal per (gene, focal) if requested.
  ## spec_focal[g, c] = focal_mean[g, c] / (focal_mean[g, c] + max_other_mean[g, c])
  ## Genes that aren't expressed in focal c get spec ~ 0 (downweighted).
  ## Reference: matches MCSD spec exponent=2 weighting; legacy mvpql_variance
  ## _decomposition_multi agg_focal_5block_specw aggregate.
  spec_lookup <- if (isTRUE(weight_by_spec_sq)) {
    ## Align ct_means columns to groups
    if (!all(groups %in% colnames(ct_means)))
      stop("ct_means columns must include all celltype groups")
    ctm <- ct_means[, groups, drop = FALSE]
    ## Align rows to gene_names
    if (!all(gene_names %in% rownames(ctm)))
      stop("ct_means rows must include all gene names from fit$U")
    ctm <- ctm[gene_names, , drop = FALSE]
    out <- matrix(0, nrow = G, ncol = length(groups),
                  dimnames = list(gene_names, groups))
    for (ci in seq_along(groups)) {
      fmean    <- ctm[, ci]
      othermax <- apply(ctm[, -ci, drop = FALSE], 1, max, na.rm = TRUE)
      out[, ci] <- fmean / pmax(fmean + othermax, 1e-9)
    }
    out
  } else NULL

  for (c_name in focal_levels) {
    c_idx <- match(c_name, groups)
    cells_in_c <- which(as.character(df$celltype) == c_name)
    if (length(cells_in_c) < 5L) next
    n_c <- length(cells_in_c)
    pat_c <- as.character(df[[patient_col]][cells_in_c])
    R_c <- if (".resp_dummy" %in% colnames(df)) as.numeric(df$.resp_dummy[cells_in_c])
            else if (!is.null(resp_term)) as.numeric(df[[paste0(resp_term, "_dummy")]][cells_in_c])
            else rep(0, n_c)
    spec_g <- if (isTRUE(weight_by_spec_sq)) spec_lookup[, c_name] else NULL
    spec_w <- if (isTRUE(weight_by_spec_sq)) spec_g^2 else NULL

    ## (1) Cell type identity (level, not variance — constant within focal c)
    ct_u2  <- fit$U[ct_col_fn(term2t_ct[["(Intercept)"]], c_idx), ]^2
    ct_vp  <- fit$se_U[ct_col_fn(term2t_ct[["(Intercept)"]], c_idx), ]^2
    ct_offset_sq <- shrink_uv(ct_u2, ct_vp)

    ## (2) Condition main effect = beta_R^2 * Var(R over focal c).
    ## R is patient-level (constant within patient), so within-patient ANOVA
    ## var_R = 0 by construction. Use total variance for this block regardless
    ## of variance_estimator setting — the only meaningful Var(R) is across
    ## the patient × focal-cell mix.
    var_R <- var_within_or_total(R_c, pat_c, "total")
    V_condition <- if (!is.na(resp_fixed_col)) fit$B[resp_fixed_col, ]^2 * var_R else numeric(G)

    ## (3) Spatial cell state — per-pair contributions
    slope_term_idx <- vapply(vars, function(v) term2t_ct[[v]], integer(1))
    slope_blups <- fit$U[ct_col_fn(slope_term_idx, c_idx), , drop=FALSE]   # K_t × G
    slope_post  <- fit$se_U[ct_col_fn(slope_term_idx, c_idx), , drop=FALSE]^2
    N_c <- as.matrix(df[cells_in_c, vars, drop=FALSE])
    var_N  <- vapply(seq_along(vars),
                     function(ti) var_within_or_total(N_c[, ti], pat_c, variance_estimator),
                     numeric(1))
    if (pair_allocation == "diagonal") {
      V_spat_pair_g <- shrink_uv(slope_blups^2, slope_post) * var_N
    } else {
      ## LMG / Shapley split: per-gene per-pair = row_sum(cov(N_c %*% diag(u_t)))[t]
      ## plus posterior-variance contribution u_se^2 * Var(N_t) on the diagonal.
      V_spat_pair_g <- matrix(0, nrow = length(vars), ncol = G)
      ## For numerator_shrinkage, we apply a per-(t,g) scaling factor
      ## s_tg = shrink_uv(u_tg^2, V_post_tg) / max(u_tg^2, eps) so that the LMG
      ## row-sum scales appropriately for each pair-gene cell.
      for (gi in seq_len(G)) {
        u_g <- slope_blups[, gi]
        contrib <- sweep(N_c, 2, u_g, FUN = "*")    # n_c × T
        cv <- stats::cov(contrib)
        rs <- rowSums(cv)
        rs <- pmax(rs, 0)                             # clamp tiny negatives
        if (numerator_shrinkage != "none") {
          u2 <- u_g^2; vp <- slope_post[, gi]
          fac <- shrink_uv(u2, vp) / pmax(u2, 1e-12)
          rs <- rs * fac
        }
        V_spat_pair_g[, gi] <- rs + slope_post[, gi] * var_N
      }
    }
    V_state_baseline <- colSums(V_spat_pair_g)

    ## (4) Condition × spatial — per-pair contributions
    has_resp <- !is.null(resp_term) && any(grepl(paste0("^", resp_term, ":"), blk_ct$term_levels))
    if (has_resp) {
      resp_term_names <- paste0(resp_term, ":", vars)
      matched  <- vapply(resp_term_names, function(nm) {
        v <- term2t_ct[[nm]]; if (is.null(v)) NA_integer_ else as.integer(v)
      }, integer(1))
      keep_v   <- which(!is.na(matched))
      resp_idx <- matched[keep_v]
      rblups <- fit$U[ct_col_fn(resp_idx, c_idx), , drop=FALSE]
      rpost  <- fit$se_U[ct_col_fn(resp_idx, c_idx), , drop=FALSE]^2
      resp_vars <- vars[keep_v]                        # neighbour labels for the surviving rows
      N_r_c <- N_c[, resp_vars, drop=FALSE] * R_c      # surviving neighbours (was vars[seq_along] — misaligned)
      var_Nr <- vapply(seq_along(resp_idx),
                       function(ti) var_within_or_total(N_r_c[, ti], pat_c, variance_estimator),
                       numeric(1))
      if (pair_allocation == "diagonal") {
        V_resp_pair_g <- shrink_uv(rblups^2, rpost) * var_Nr
      } else {
        V_resp_pair_g <- matrix(0, nrow = length(resp_idx), ncol = G)
        for (gi in seq_len(G)) {
          u_g <- rblups[, gi]
          contrib_r <- sweep(N_r_c, 2, u_g, FUN = "*")
          cv_r <- stats::cov(contrib_r)
          rs_r <- pmax(rowSums(cv_r), 0)
          if (numerator_shrinkage != "none") {
            u2 <- u_g^2; vp <- rpost[, gi]
            fac <- shrink_uv(u2, vp) / pmax(u2, 1e-12)
            rs_r <- rs_r * fac
          }
          V_resp_pair_g[, gi] <- rs_r + rpost[, gi] * var_Nr
        }
      }
      V_state_responder <- colSums(V_resp_pair_g)
    } else {
      V_resp_pair_g <- matrix(0, length(slope_term_idx), G)
      resp_vars <- vars
      V_state_responder <- numeric(G)
    }

    ## (5) imageID block — Var across focal-c cells of Z_img × u_img per gene
    V_image <- if (!is.null(blk_img)) {
      img_cols <- (blk_img$col_offset + 1L):(blk_img$col_offset + blk_img$n_cols)
      Z_img_c <- Z_full[cells_in_c, img_cols, drop = FALSE]
      U_img <- fit$U[img_cols, , drop=FALSE]
      eta_img_mat <- as.matrix(Z_img_c %*% U_img)
      var_within_mat(eta_img_mat, pat_c, variance_estimator)
    } else numeric(G)

    ## (6) Residual dispersion
    mu_c <- as.matrix(fit$mu[cells_in_c, , drop=FALSE])
    mu_bar <- colMeans(mu_c, na.rm=TRUE)
    V_disp <- log(1 + (1 + alpha_g) / pmax(mu_bar, 1e-9))

    Total_g <- ct_offset_sq + V_condition + V_state_baseline +
               V_state_responder + V_image + V_disp

    ## Per-focal aggregate %
    safe_pct <- function(num) 100 * num / pmax(Total_g, 1e-12)
    aggfun_blk <- function(num_vec) {
      if (aggregate == "pooled") {
        agg_pooled(num_vec, Total_g)
      } else if (isTRUE(weight_by_spec_sq)) {
        agg_weighted(safe_pct(num_vec), spec_w)
      } else {
        agg_unweighted(safe_pct(num_vec))
      }
    }
    focal_rows[[length(focal_rows) + 1L]] <- tibble::tibble(
      focal = c_name, n_focal = n_c,
      `Cell type %`             = aggfun_blk(ct_offset_sq),
      `Condition main %`        = aggfun_blk(V_condition),
      `Spatial state %`         = aggfun_blk(V_state_baseline),
      `Resp x Spatial %`        = aggfun_blk(V_state_responder),
      `imageID %`               = aggfun_blk(V_image),
      `Residual %`              = aggfun_blk(V_disp))

    ## Per-pair % for spatial and condition × spatial blocks
    for (block in c("spatial","responder_spatial")) {
      M <- if (block == "spatial") V_spat_pair_g else V_resp_pair_g
      ## The responder block keeps only the neighbours whose interaction term
      ## matched, so it can have fewer rows than vars. Label rows from the
      ## matching vector, never from seq_along(vars).
      nbr <- if (block == "spatial") vars else resp_vars
      for (ti in seq_along(nbr)) {
        if (aggregate == "pooled") {
          pp <- agg_pooled(M[ti, ], Total_g)
        } else if (isTRUE(weight_by_spec_sq)) {
          pp <- agg_weighted(safe_pct(M[ti, ]), spec_w)
        } else {
          pp <- agg_unweighted(safe_pct(M[ti, ]))
        }
        pair_rows[[length(pair_rows) + 1L]] <- tibble::tibble(
          focal = c_name, neighbour = nbr[ti], block = block,
          pair_pct = pp,
          pair_pct_med = stats::median(safe_pct(M[ti, ]), na.rm=TRUE))
      }
    }
    if (verbose) message(sprintf("[v6] focal '%s' done (n=%d)", c_name, n_c))
  }
  list(focal_summary = dplyr::bind_rows(focal_rows),
       pair_pct      = dplyr::bind_rows(pair_rows))
}
