## Pipeline-step mashr shrinkage of per-celltype neighbour-slope BLUPs.
##
## Two entry points:
##   apply_mashr_shrinkage(results, focals, neighbours, resp_term = NULL)
##     -> long-form tibble (gene, focal, neighbour, estimate, std.error,
##                          estimate_shrunk, sd_shrunk, lfsr)
##     One mashr fit per neighbour slice, focals as conditions.
##
##   merge_shrunken_into_ranvals(ran_vals, shrunken_long, resp_term = NULL)
##     -> ran_vals + extra columns (estimate_shrunk, sd_shrunk, lfsr) joined
##     on (gene, level=focal, term=neighbour). Original estimate / std.error
##     untouched so downstream code can choose either.
##
## In melanoma (multi-sample) pass `resp_term = "ResponderPD"` and the
## function shrinks the ResponderPD:<neighbour> interaction BLUPs. In breast
## cancer (single-sample) leave `resp_term = NULL` and it shrinks the bare
## <neighbour> proximity-slope BLUPs.

suppressPackageStartupMessages({
})

# Null correlation between conditions, or NULL when it is not usable.
#
# estimate_null_correlation_simple() takes the empirical correlation of the z
# scores over effects with |z| < 2 in every condition. Sparse focal-neighbour
# pairs are zeroed by the kernel drop without their column being removed, so
# those conditions have a constant z column; cov2cor then divides by a zero
# variance and returns NaN for the WHOLE matrix. A constant column carries no
# estimable correlation, so the estimate is taken over the conditions that vary
# and the degenerate ones are left uncorrelated, which is the correct limit.
.null_correlation <- function(Bhat, Shat, target_term, verbose = TRUE) {
  Z    <- Bhat / Shat
  sd_z <- apply(Z, 2, function(v) stats::sd(v[is.finite(v)]))
  live <- which(is.finite(sd_z) & sd_z > 0)
  if (length(live) < 2L) {
    if (verbose) message("  [mashr] term '", target_term,
            "': fewer than two varying conditions, using independence")
    return(NULL)
  }
  V_live <- tryCatch(
    estimate_null_correlation_simple(
      mash_set_data(Bhat[, live, drop = FALSE], Shat[, live, drop = FALSE])),
    error = function(e) NULL)
  usable <- !is.null(V_live) && all(is.finite(V_live)) &&
    isTRUE(all.equal(V_live, t(V_live), tolerance = 1e-8)) &&
    min(eigen(V_live, symmetric = TRUE, only.values = TRUE)$values) > 1e-8
  if (!usable) {
    if (verbose) message("  [mashr] term '", target_term,
            "': null correlation not estimable, using independence")
    return(NULL)
  }
  V <- diag(ncol(Bhat))
  dimnames(V) <- list(colnames(Bhat), colnames(Bhat))
  V[live, live] <- V_live
  V
}


# Run mashr per-neighbour and return shrunken BLUPs in long form
#
# @param results PACE fit results from extractRandomEffects(*)
# @param focals  Character vector of focal cell-type levels (the columns
#                of each per-neighbour Bhat / Shat matrix).
# @param neighbours Character vector of neighbour cell types whose slopes
#                   are to be shrunk. The actual `term` looked up in
#                   ran_vals is `paste0(resp_term, ":", neighbour)` if
#                   resp_term is non-NULL, otherwise just `neighbour`.
# @param resp_term Responder term name for the multi-sample model
#                  (e.g., "ResponderPD"); NULL for single-sample.
# @param min_significant Minimum number of "significant" genes required
#                        before adding data-driven covariances; otherwise
#                        canonical-only.
# @param null_correlation Estimate the correlation between conditions under
#                        the null and pass it to mash as V. The focal cell
#                        types are conditions in the same per-gene solve, so
#                        they are correlated, and mash's default V = I
#                        mis-calibrates lfsr. FALSE restores independence.
# @return Long-form tibble.
apply_mashr_shrinkage <- function(results, focals, neighbours,
                                  resp_term = NULL,
                                  min_significant = 10,
                                  data_driven = TRUE,
                                  null_correlation = TRUE,
                                   ## How many R processes shrink the neighbour
                                   ## slices at once. 1 (default) keeps the
                                   ## serial order, which is what the fixtures
                                   ## were made with; see .pace_shrink_slices()
                                   ## for what changes above 1.
                                   shrink_threads = 1L,
                                   ## Progress reporting. FALSE silences both
                                   ## this function's messages and mashr's own,
                                   ## which it prints unconditionally otherwise.
                                   verbose = TRUE) {
  ## data_driven=FALSE -> canonical covariances only (skip the slow cov_pca/cov_ed
  ## extreme-deconvolution step). Much faster when there are many neighbour terms
  ## (e.g. 16-type fits); shrinkage is slightly less adaptive but driver rankings
  ## are nearly identical. Default TRUE preserves the canonical manuscript behaviour.
  rv_all <- get_ran_vals(results, well_fit_only = TRUE)
  build_term <- function(nb) {
    if (is.null(resp_term)) nb else paste0(resp_term, ":", nb)
  }
  se_cap <- 50  # any SE above this is a non-converged BLUP, drop
  se_floor <- 1e-4  # below this both Bhat and Shat are effectively zero

  ## Pass 1: build every neighbour slice and find the genes each one can use.
  slices <- list()
  for (nb in neighbours) {
    target_term <- build_term(nb)
    sub <- rv_all |>
      dplyr::filter(term == target_term, group == "celltype",
                    level %in% focals)
    if (nrow(sub) == 0) {
      if (verbose) message("  [mashr] no ran_vals for term '", target_term, "', skipping")
      next
    }

    Bhat <- sub |>
      dplyr::select(gene, focal = level, estimate) |>
      tidyr::pivot_wider(names_from = focal, values_from = estimate)
    Shat <- sub |>
      dplyr::select(gene, focal = level, std.error) |>
      tidyr::pivot_wider(names_from = focal, values_from = std.error)

    # Drop genes with NA / non-positive SE in any focal column
    cols_present <- intersect(focals, colnames(Bhat))
    if (length(cols_present) < 2) {
      if (verbose) message("  [mashr] term '", target_term,
              "' has <2 focal columns, skipping")
      next
    }
    genes <- Bhat$gene
    Bhat <- as.matrix(Bhat[, cols_present, drop = FALSE]); rownames(Bhat) <- genes
    Shat <- as.matrix(Shat[, cols_present, drop = FALSE]); rownames(Shat) <- genes
    ## Strict filter: require finite Bhat / Shat, positive Shat, and Shat
    ## not pathologically large (boundary-singular fits sometimes return
    ## absurd SEs that destabilise mash's likelihood matrix).
    keep <- rowSums(!is.finite(Bhat)) == 0 &
            rowSums(!is.finite(Shat)) == 0 &
            rowSums(Shat <= 0)        == 0 &
            rowSums(Shat > se_cap)    == 0 &
            ## drop rows where every focal's BLUP collapsed to ~zero with
            ## ~zero SE — these are degenerate boundary fits that crash
            ## mashr's likelihood-matrix check
            !(rowSums(abs(Bhat) < se_floor & Shat < se_floor) == ncol(Bhat))
    ## A slice with too few usable genes is skipped on its own; letting it into
    ## the common gene set below would empty every other slice too.
    if (sum(keep) < 5) {
      if (verbose) message("  [mashr] term '", target_term,
              "' has <5 well-fit genes after filter, skipping")
      next
    }
    slices[[nb]] <- list(target_term = target_term, Bhat = Bhat, Shat = Shat,
                         usable_genes = genes[keep])
  }
  if (length(slices) == 0) return(.empty_shrunken_long())

  ## Every slice is shrunk over the SAME genes. mash calibrates its null
  ## against the genes it is given, so slices filtered separately would carry
  ## lfsr on different bases, and one global lfsr cut across them would not be
  ## comparable. A gene unusable in any slice is dropped from all of them.
  common_genes <- Reduce(intersect, lapply(slices, `[[`, "usable_genes"))
  for (nb in names(slices)) {
    n_dropped <- length(setdiff(slices[[nb]]$usable_genes, common_genes))
    if (n_dropped > 0)
      if (verbose) message("  [mashr] term '", slices[[nb]]$target_term, "': ", n_dropped,
              " genes dropped to keep the gene set common across neighbours")
  }

  ## Pass 2: shrink each slice over the common genes. The slices are
  ## independent -- one mash fit each -- so `shrink_threads` can run them in
  ## separate R processes. mashr is R code, so processes are the only way to
  ## overlap them; the fit itself needs none of this and starts no workers.
  shrink_one <- function(nb) {
    target_term <- slices[[nb]]$target_term
    Bhat <- slices[[nb]]$Bhat[common_genes, , drop = FALSE]
    Shat <- slices[[nb]]$Shat[common_genes, , drop = FALSE]
    ## Floor any remaining tiny SEs to a small positive number to keep
    ## mashr's likelihood matrix well-conditioned.
    Shat <- pmax(Shat, se_floor)
    if (nrow(Bhat) < 5) {
      if (verbose) message("  [mashr] term '", target_term,
              "' has <5 genes common to every neighbour, skipping")
      return(NULL)
    }

    ## The focal cell types are conditions estimated in the SAME per-gene
    ## solve, so their errors are correlated. mash's default V = I treats them
    ## as independent, which mis-calibrates lfsr. Estimate the null correlation
    ## from the effects that look null in every condition (mashr's documented
    ## two-step workflow), and fall back to independence when that estimate is
    ## not usable, which happens on slices with too few null genes.
    V <- if (isTRUE(null_correlation)) .null_correlation(Bhat, Shat, target_term, verbose = verbose) else NULL
    data <- if (is.null(V)) mash_set_data(Bhat, Shat)
            else mash_set_data(Bhat, Shat, V = V)
    U_c  <- cov_canonical(data)
    ## Be defensive about data-driven covariances: cov_pca / cov_ed can
    ## produce non-PSD components on small / degenerate slices, which then
    ## break mash's likelihood-matrix check. Try canonical + ED; if it
    ## fails, fall back to canonical only.
    m <- tryCatch({
      m_c <- mash(data, U_c, verbose = verbose)
      top_idx <- get_significant_results(m_c)
      if (data_driven && length(top_idx) >= min_significant) {
        U_pca <- cov_pca(data, npc = min(5, length(top_idx) - 1),
                         subset = top_idx)
        U_ed  <- cov_ed(data, U_pca, subset = top_idx)
        mash(data, c(U_c, U_ed), verbose = verbose)
      } else {
        m_c
      }
    }, error = function(e) {
      if (verbose) message("  [mashr] data-driven cov failed for term '", target_term,
              "' (", conditionMessage(e), "); falling back to canonical only")
      tryCatch(mash(data, U_c, verbose = verbose),
               error = function(e2) {
                 if (verbose) message("  [mashr] canonical-only also failed (",
                         conditionMessage(e2), "); skipping")
                 NULL
               })
    })
    if (is.null(m)) return(NULL)

    pm   <- get_pm(m)
    psd  <- get_psd(m)
    lfsr <- get_lfsr(m)

    result <- tibble::tibble(
      gene = rownames(pm)[row(pm)],
      focal = colnames(pm)[col(pm)],
      neighbour = nb,
      term      = target_term,
      estimate  = as.numeric(Bhat),
      std.error = as.numeric(Shat),
      estimate_shrunk = as.numeric(pm),
      sd_shrunk       = as.numeric(psd),
      lfsr            = as.numeric(lfsr)
    )
    fsr <- expected_false_sign(as.numeric(lfsr))
    if (verbose) message(sprintf("  [mashr] %s: %d genes shrunk; gene-focal calls (lfsr<0.05) = %d, expected false sign = %.1f",
                    target_term, nrow(Bhat), fsr$n_calls, fsr$expected_false_sign))
    result
  }

  out <- .pace_shrink_slices(names(slices), shrink_one, shrink_threads, verbose = verbose)
  out <- out[!vapply(out, is.null, logical(1))]

  if (length(out) == 0) return(.empty_shrunken_long())
  dplyr::bind_rows(out)
}


# Run one function over the neighbour slices, in this process or in several.
#
# The slices are independent, and about 90% of the shrinkage runs inside them,
# so this is where the wall time is. Processes, not threads: mashr is R code.
# BiocParallel's MulticoreParam forks, which is cheap here (the slices are
# already built, and copy-on-write shares them); Windows has no fork, so it
# falls back to running them in order.
#
# RNG: mash() draws from the stream itself, and cov_pca() draws again for its
# irlba starting vectors, so a serial run consumes the stream slice by slice in
# an order no parallel run can reproduce -- with or without `data_driven`. The
# parallel path therefore takes its randomness from BiocParallel's RNGseed,
# which hands each task its own L'Ecuyer-CMRG substream: reproducible run to
# run, but NOT equal to the serial result. The serial default is what the
# fixtures were made with, and the measured difference is in the design note.
#
# Seeding each worker with set.seed() would do the same job and is what this
# did first, but package code must not call set.seed(): it overwrites whatever
# stream the caller had set up, and BiocCheck flags it. RNGseed is the
# supported way to get the same determinism without touching global state.
.pace_shrink_slices <- function(names_of_slices, shrink_one, shrink_threads,
                                verbose = TRUE) {
  threads <- .pace_thread_count(shrink_threads)
  if (threads == 1L || length(names_of_slices) < 2L) {
    return(stats::setNames(lapply(names_of_slices, shrink_one), names_of_slices))
  }
  if (.Platform$OS.type == "windows") {
    if (verbose) message("  [mashr] shrink_threads > 1 needs a forking platform; running the slices in order")
    return(stats::setNames(lapply(names_of_slices, shrink_one), names_of_slices))
  }
  ## Collect once here rather than once inside every worker: the children
  ## inherit this heap, and its garbage with it.
  invisible(gc(full = TRUE, verbose = FALSE))
  ## One task per slice, so a slow slice does not hold a worker's whole share:
  ## the slices differ several-fold in cost. RNGseed makes the substreams, and
  ## therefore the result, the same on every run.
  param <- BiocParallel::MulticoreParam(workers = min(threads, length(names_of_slices)),
                                        tasks = length(names_of_slices),
                                        RNGseed = 1L)
  stats::setNames(BiocParallel::bplapply(names_of_slices, shrink_one, BPPARAM = param),
                  names_of_slices)
}

.empty_shrunken_long <- function() {
  tibble::tibble(gene = character(), focal = character(),
                 neighbour = character(), term = character(),
                 estimate = numeric(), std.error = numeric(),
                 estimate_shrunk = numeric(), sd_shrunk = numeric(),
                 lfsr = numeric())
}

# Expected number of false sign calls among the effects called at lfsr < thresh
#
# lfsr is the posterior probability that a called effect has the wrong sign, so
# its sum over the calls is the posterior expected number of sign errors among
# them (Stephens 2017, Biostatistics 18:275), and its mean is the estimated
# false sign rate. Both are only as calibrated as the lfsr they are built from.
#
# @param lfsr Numeric vector of lfsr values (NA ignored).
# @param thresh The calling threshold.
# @return A list: `n_calls`, `expected_false_sign`, `false_sign_rate`.
expected_false_sign <- function(lfsr, thresh = 0.05) {
  called <- lfsr[!is.na(lfsr) & lfsr < thresh]
  n_calls <- length(called)
  expected <- sum(called)
  list(n_calls = n_calls,
       expected_false_sign = expected,
       false_sign_rate = if (n_calls > 0) expected / n_calls else NA_real_)
}

# Left-join shrunken values into ran_vals
#
# @param ran_vals  Output of get_ran_vals(results, ...).
# @param shrunken_long Output of apply_mashr_shrinkage(...).
# @return ran_vals plus columns estimate_shrunk, sd_shrunk, lfsr (NA where
#   the (gene, level, term) combination was not shrunk).
merge_shrunken_into_ranvals <- function(ran_vals, shrunken_long) {
  if (nrow(shrunken_long) == 0) {
    return(ran_vals |> dplyr::mutate(estimate_shrunk = NA_real_,
                                     sd_shrunk       = NA_real_,
                                     lfsr            = NA_real_))
  }
  ran_vals |>
    dplyr::left_join(shrunken_long |>
                       dplyr::select(gene, level = focal, term,
                                     estimate_shrunk, sd_shrunk, lfsr),
                     by = c("gene", "level", "term"))
}
