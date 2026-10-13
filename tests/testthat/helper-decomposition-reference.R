# The dense single-frame decomposition as implemented at f55d976 (R/engine-decomp.R),
# kept verbatim as the reference for the sparse-statistics rewrite.

reference_single_frame_decomp_obs <- function(Y, celltype, nCount, gene_focal_block) {
  genes <- colnames(Y); ct <- as.character(celltype); TYPES <- sort(unique(ct))
  Ylog <- log1p(Y * (1e4 / nCount))                      # n x G observed log1p CP10k
  global_mean <- as.numeric(Matrix::colMeans(Ylog)); names(global_mean) <- genes
  gf <- gene_focal_block
  # Condition cohorts split the spatial component into baseline + responder
  # (gene_focal_5block has V_state_responder); no-condition cohorts have V_state.
  has_resp <- "V_state_responder" %in% names(gf)

  out <- vector("list", length(TYPES)); k <- 0L
  for (fc_type in TYPES) {
    fc <- which(ct == fc_type); n_fc <- length(fc); if (n_fc < 5) next
    Yf <- as.matrix(Ylog[fc, , drop = FALSE])
    focal_mean <- colMeans(Yf)
    SS_within  <- colSums(sweep(Yf, 2, focal_mean, "-")^2)     # per gene, observed
    SS_lineage <- n_fc * (focal_mean - global_mean)^2          # per gene, observed
    den <- SS_lineage + SS_within
    ## fit within-proportions (locked gene_focal block) for this focal
    g4 <- gf[gf$focal == fc_type, , drop = FALSE]
    idx <- match(genes, g4$gene)
    vstate <- if (has_resp) g4$V_state_baseline[idx] else g4$V_state[idx]
    vresp  <- if (has_resp) g4$V_state_responder[idx] else 0
    vspill <- g4$V_spill[idx]; vdisp <- g4$V_disp[idx]
    vt <- vstate + vresp + vspill + vdisp
    p_sp   <- ifelse(vt > 0, vstate / vt, 0)
    p_resp <- ifelse(vt > 0, vresp  / vt, 0)
    p_bl   <- ifelse(vt > 0, vspill / vt, 0)
    p_rs   <- ifelse(vt > 0, vdisp  / vt, 0)
    k <- k + 1L
    row <- data.frame(
      focal = fc_type, gene = genes,
      `Cell type %` = 100 * SS_lineage       / den,
      `Spatial %`   = 100 * SS_within * p_sp / den,
      check.names = FALSE, stringsAsFactors = FALSE)
    if (has_resp) row[["Responder spatial %"]] <- 100 * SS_within * p_resp / den
    row[["Spillover %"]] <- 100 * SS_within * p_bl / den
    row[["Residual %"]]  <- 100 * SS_within * p_rs / den
    row$SS_lineage <- SS_lineage; row$SS_within <- SS_within
    row$denom <- den; row$n_focal <- n_fc
    out[[k]] <- row
  }
  do.call(rbind, out)
}

# A fit made by this version keeps only the per-cell-type statistics, so the
# per-cell reference below gets the n x G matrices rebuilt from the fit and the
# counts -- the same rebuild paceDecompose() uses for an older stripped fit.
with_rebuilt_matrices <- function(fit, spe) {
  if (!is.null(fit@fit$mu)) return(fit)      # the dense solver still stores them
  inputs <- PACE:::.pace_mu_inputs(fit, spe)
  parts <- PACE:::.pace_mu_block(fit, inputs)
  if (is.null(parts$mu_spill)) {
    fit@fit$mu <- parts$mu_bio
    fit@fit$technical_offset_mat <- matrix(0, nrow(parts$mu_bio), ncol(parts$mu_bio))
  } else {
    fit@fit$mu <- pmax(parts$mu_bio + parts$mu_spill, 1e-6)
    fit@fit$technical_offset_mat <- log1p(parts$mu_spill / parts$mu_bio)
  }
  fit@fit$bleed_offset_mat <- fit@fit$technical_offset_mat
  fit
}

# The f55d976 decomposition wrapper: the per-cell loop of
# mvpql_variance_decomposition_multi() on a fit that stores mu and the technical
# offset, plus the 4-block view.
reference_pace_decompose <- function(fit, df, Y, types, X_fixed, resp_term = NULL) {
  dec <- PACE:::mvpql_variance_decomposition_multi(
    fit = fit, df = df, Y = Y, vars = types,
    X_fixed = X_fixed, resp_term = resp_term, focal_levels = types)
  g5 <- dec$gene_focal_5block
  total4 <- with(g5, celltype_offset_sq + V_state_baseline + V_state_responder +
                     V_spill + V_disp)
  dec$gene_focal_4block <- tibble::tibble(
    focal = g5$focal,
    gene  = g5$gene,
    `Cell type %`     = 100 * g5$celltype_offset_sq / pmax(total4, 1e-12),
    `Spatial state %` = 100 * (g5$V_state_baseline + g5$V_state_responder) / pmax(total4, 1e-12),
    `Spillover %`     = 100 * g5$V_spill / pmax(total4, 1e-12),
    `Residual %`      = 100 * g5$V_disp  / pmax(total4, 1e-12),
    celltype_offset_sq = g5$celltype_offset_sq,
    V_state = g5$V_state_baseline + g5$V_state_responder,
    V_spill = g5$V_spill,
    V_disp  = g5$V_disp,
    Total   = total4,
    spec       = g5$spec,
    focal_mean = g5$focal_mean)
  dec
}

# How each entry is missing, so NaN, NA, +Inf and -Inf cannot stand in for one
# another in the comparison below.
.pace_value_kind <- function(x) {
  kind <- rep("finite", length(x))
  kind[is.nan(x)] <- "NaN"
  kind[is.na(x) & !is.nan(x)] <- "NA"
  kind[!is.na(x) & x == Inf] <- "Inf"
  kind[!is.na(x) & x == -Inf] <- "-Inf"
  kind
}

# Two objects are equal when their non-numeric parts are identical and every
# numeric entry satisfies |a - b| <= tol_abs + tol_rel * |b|, with the absolute
# floor scaled to the column (exact zeros on one side, 1e-17 residue on the other).
# tol_rel = 1e-10 leaves more than 200x headroom over the measured identity errors
# (<= 4.6e-13); tol_abs = 1e-14 x max(1, max |b|) covers quadratic-form residue
# where the per-cell variance is exactly zero.
expect_tables_close <- function(new, reference, tol_rel = 1e-10, tol_abs_scale = 1e-14, label = "") {
  if (is.data.frame(reference)) {
    expect_identical(names(new), names(reference), label = label)
    expect_identical(nrow(new), nrow(reference), label = label)
    for (column in names(reference)) {
      a <- new[[column]]
      b <- reference[[column]]
      if (is.numeric(b)) {
        expect_identical(names(a), names(b), label = paste(label, column, "names"))
        a <- as.numeric(a)
        b <- as.numeric(b)
        ## NaN, NA, +Inf and -Inf are told apart: is.finite() alone let a NaN
        ## stand in for an NA, and +Inf for -Inf.
        expect_identical(.pace_value_kind(a), .pace_value_kind(b),
                         label = paste(label, column, "missing pattern"))
        finite <- is.finite(b)
        tol_abs <- tol_abs_scale * max(1, abs(b[finite]))
        within <- abs(a[finite] - b[finite]) <= tol_abs + tol_rel * abs(b[finite])
        expect_true(all(within), label = paste(label, column))
      } else {
        expect_identical(a, b, label = paste(label, column))
      }
    }
    return(invisible(TRUE))
  }
  if (is.list(reference)) {
    expect_identical(names(new), names(reference), label = label)
    for (key in names(reference)) {
      expect_tables_close(new[[key]], reference[[key]], tol_rel, tol_abs_scale, paste(label, key))
    }
    return(invisible(TRUE))
  }
  if (is.numeric(reference)) {
    expect_tables_close(data.frame(value = as.numeric(new)), data.frame(value = as.numeric(reference)),
                        tol_rel, tol_abs_scale, label)
    return(invisible(TRUE))
  }
  expect_identical(new, reference, label = label)
}

# The f55d976 driver-score loop (R/pace-core.R), kept verbatim as the reference
# for the compiled driver scores. Only the function name changed.
reference_pace_top_drivers <- function(fit, shrunken_long, dec, types, mu_means, pairs = NULL,
                             resp_term = NULL, resp_dummy = NULL) {
  ## default: every ordered focal != neighbour pair (a chosen subset is only a
  ## reporting convenience, not part of the method).
  if (is.null(pairs)) {
    pairs <- list()
    for (fc in types) for (nc in types) if (fc != nc) pairs <- c(pairs, list(c(fc, nc)))
  }
  g5 <- dec$gene_focal_5block
  Z_re <- fit$re_meta$Z
  cells_by_ct <- lapply(types, function(c) which(Z_re[, paste0(c, "::(Intercept)")] != 0))
  names(cells_by_ct) <- types
  alpha_g <- pmax(fit$alpha, 0)
  gene_names_fit <- colnames(fit$U)

  out <- list()
  for (p in pairs) {
    fc <- p[1]
    nc <- p[2]
    pk <- paste(fc, nc, sep = "_")
    target_term <- if (is.null(resp_term)) nc else paste0(resp_term, ":", nc)
    s <- shrunken_long |>
      dplyr::filter(focal == fc, neighbour == nc, term == target_term) |>
      dplyr::distinct(gene, .keep_all = TRUE)
    fm <- g5 |>
      dplyr::filter(focal == fc) |>
      dplyr::select(gene, spec, focal_mean) |>
      dplyr::distinct(gene, .keep_all = TRUE)
    if (!nrow(s) || !nrow(fm)) next

    cells_c <- cells_by_ct[[fc]]
    col_N <- paste0(fc, "::", nc)
    if (!col_N %in% colnames(Z_re)) {
      out[[pk]] <- list(scores = s[0, ], status = "dropped (n_eff)",
                        expected_false_sign = 0, false_sign_rate = NA_real_)
      next
    }
    N_t <- as.numeric(Z_re[cells_c, col_N])
    var_N <- stats::var(N_t, na.rm = TRUE)
    mu_bar_per_gene <- mu_means[fc, ]
    names(alpha_g) <- gene_names_fit
    names(mu_bar_per_gene) <- gene_names_fit

    base <- s |>
      dplyr::inner_join(fm, by = "gene") |>
      dplyr::mutate(
        driver_score       = (estimate_shrunk^2) * (spec^2) * pmax(focal_mean, 0),
        driver_score_spec4 = (estimate_shrunk^2) * (spec^4) * pmax(focal_mean, 0),
        mu_bar  = mu_bar_per_gene[gene],
        alpha   = alpha_g[gene],
        V_resid = log(1 + (1 + pmax(alpha, 0)) / pmax(mu_bar, 1e-6)))
    if (is.null(resp_term)) {
      ## baseline neighbour effect only (no condition).
      res_all <- base |>
        dplyr::mutate(V_S = estimate_shrunk^2 * var_N,
                      V_total = V_S + V_resid,
                      R2_S = V_S / pmax(V_total, 1e-12))
      res <- res_all |>
        dplyr::filter(lfsr < 0.05) |>
        dplyr::arrange(dplyr::desc(driver_score)) |>
        dplyr::mutate(rank = dplyr::row_number()) |>
        dplyr::rename(b_clean = estimate_shrunk) |>
        dplyr::select(rank, gene, driver_score, driver_score_spec4, b_clean, spec, focal_mean,
                      R2_S, mu_bar, alpha, V_S, V_resid, V_total, lfsr, sd_shrunk)
    } else {
      ## condition x spatial: baseline slope V_S (from the raw BLUP u) plus the
      ## condition-interaction slope V_RxS (from the shrunken estimate).
      col_R <- match(paste0(fc, "::", resp_term), colnames(Z_re))
      R_c <- if (is.na(col_R)) {
        if (is.null(resp_dummy)) rep(0, length(cells_c)) else resp_dummy[cells_c]
      } else as.numeric(Z_re[cells_c, col_R])
      var_RN <- stats::var(R_c * N_t, na.rm = TRUE)
      row_u <- match(col_N, rownames(fit$U))
      u_vec <- if (is.na(row_u))
        setNames(rep(0, length(gene_names_fit)), gene_names_fit)
      else setNames(as.numeric(fit$U[row_u, ]), gene_names_fit)
      res_all <- base |>
        dplyr::mutate(u_raw   = u_vec[gene],
                      V_S     = u_raw^2 * var_N,
                      V_RxS   = estimate_shrunk^2 * var_RN,
                      V_total = V_S + V_RxS + V_resid,
                      R2_S    = V_S   / pmax(V_total, 1e-12),
                      R2_RxS  = V_RxS / pmax(V_total, 1e-12))
      res <- res_all |>
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
    fsr <- PACE:::expected_false_sign(res$lfsr)
    out[[pk]] <- list(scores = res, status = status,
                      expected_false_sign = fsr$expected_false_sign,
                      false_sign_rate = fsr$false_sign_rate)
  }
  out
}

# The f55d976 Pratt pair attribution (R/plots.R), kept verbatim as the reference
# for the compiled attribution. Only the function name changed.
reference_pair_variance_pratt <- function(mv, cond_prefix = NULL, focals = NULL,
                                     cohort_label = "cohort", block_label = NULL) {
  fit <- mv$fit
  Z <- fit$re_meta$Z
  gn <- mv$gene_set
  colnames(fit$U) <- gn
  TYPES <- if (!is.null(fit$re_meta$blocks))
             fit$re_meta$blocks[[1]]$group_levels
           else fit$re_meta$group_levels
  if (is.null(focals)) focals <- TYPES
  if (is.null(block_label))
    block_label <- if (is.null(cond_prefix)) "Spatial" else "RxS"

  pair_rows <- list()
  focal_rows <- list()

  for (fc in focals) {
    fc_int <- paste0(fc, "::(Intercept)")
    if (!(fc_int %in% colnames(Z))) next
    cells <- which(as.numeric(Z[, fc_int]) != 0)
    if (length(cells) < 50) next
    term_names <- if (is.null(cond_prefix))
                    paste0(fc, "::", TYPES)
                  else paste0(fc, "::", cond_prefix, ":", TYPES)
    keep <- term_names %in% colnames(Z) & term_names %in% rownames(fit$U)
    if (!any(keep)) next
    tn <- term_names[keep]
    tt <- TYPES[keep]
    Z_fc <- as.matrix(Z[cells, tn, drop = FALSE])
    Sigma_K <- stats::cov(Z_fc)
    U_c <- as.matrix(fit$U[tn, , drop = FALSE]); U_c[!is.finite(U_c)] <- 0
    SU <- Sigma_K %*% U_c
    V_pair_gene <- U_c * SU
    V_pair_t <- as.numeric(rowSums(V_pair_gene))
    V_pratt_total <- sum(V_pair_t)
    V_diag_t <- as.numeric(rowSums(U_c^2) * diag(Sigma_K))
    V_diag_total <- sum(V_diag_t)

    for (i in seq_along(tt)) {
      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        cohort = cohort_label, block = block_label,
        focal = fc, neighbour = tt[i],
        V_pair_pratt = V_pair_t[i],
        V_pair_diag  = V_diag_t[i],
        within_focal_share_pct = 100 * V_pair_t[i] / V_pratt_total,
        within_focal_diag_pct  = 100 * V_diag_t[i] / V_diag_total,
        sign = sign(V_pair_t[i]),
        stringsAsFactors = FALSE)
    }
    focal_rows[[length(focal_rows) + 1]] <- data.frame(
      cohort = cohort_label, block = block_label, focal = fc,
      n_cells = length(cells),
      V_block_pratt = V_pratt_total,
      V_block_diag  = V_diag_total,
      cross_cov_pct = 100 * (V_pratt_total - V_diag_total) / V_pratt_total,
      n_negative_pairs = sum(V_pair_t < 0),
      stringsAsFactors = FALSE)
  }
  list(pair_long = do.call(rbind, pair_rows),
       focal_summary = do.call(rbind, focal_rows))
}
