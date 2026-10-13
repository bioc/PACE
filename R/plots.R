# Manuscript-style decomposition, pairwise, and driver plots.
#
# The three engine functions below (pace_pair_variance_pratt,
# plot_stacked_bar_focal_no_resp, make_pair_composite) are ported verbatim from
# the PACE analysis code (scripts/helpers/pace_pair_variance_pratt.R and
# scripts/core/06a-plots-decomp.R). Only the data plumbing changed: the exported
# wrappers read from a PACEFit rather than a bare `mv` bundle, so the figures
# reproduce the manuscript style exactly.

## ---------------------------------------------------------------------------
## Reconstruct the minimal `mv` bundle the ported functions expect.
## ---------------------------------------------------------------------------
.pace_mv_shim <- function(object) {
  gs <- colnames(object@fit$U)
  if (is.null(gs)) gs <- as.character(seq_len(ncol(object@fit$U)))
  list(fit = object@fit,
       gene_set = gs,
       mcsd_canonical = object@topDrivers,
       decomposition = list(
         gene_focal_single_frame = object@varianceDecomposition$perGene))
}

## SS-weighted per-focal block table (Goldstein pooling), from the single frame:
## each block's percentage back to a sum of squares, pooled over genes, and back
## to a percentage. The pooling is pace::single_frame_focal_blocks; R names the
## blocks and fixes the focal order dplyr::group_by() would have produced.
.pace_focal_blocks <- function(sf) {
  has_resp <- "Responder spatial %" %in% names(sf)
  block_names <- c("Cell type", "Spatial cell state", "Responder spatial state",
                   "Spillover", "Residual")
  share_names <- c("Cell type %", "Spatial %", "Responder spatial %", "Spillover %", "Residual %")
  shares <- matrix(0, nrow(sf), length(share_names))
  for (b in seq_along(share_names)) {
    if (share_names[b] %in% names(sf)) shares[, b] <- sf[[share_names[b]]]
  }
  focal_levels <- dplyr::group_keys(dplyr::group_by(sf, .data$focal))$focal
  pooled <- pace_single_frame_focal_blocks_cpp(
    as.integer(match(sf$focal, focal_levels) - 1L), length(focal_levels), shares, sf$denom)
  fb <- tibble::tibble(focal = focal_levels)
  for (b in seq_along(block_names)) fb[[block_names[b]]] <- pooled[, b]
  if (!has_resp) fb[["Responder spatial state"]] <- NULL
  fb
}

## ---------------------------------------------------------------------------
## Pratt 1987 signed pair-level variance decomposition (verbatim engine port).
## ---------------------------------------------------------------------------
pace_pair_variance_pratt <- function(mv, cond_prefix = NULL, focals = NULL,
                                     cohort_label = "cohort", block_label = NULL) {
  fit <- mv$fit
  Z <- .pace_as_dgc(fit$re_meta$Z)
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
    cells <- pace_column_nonzero_rows_cpp(Z, match(fc_int, colnames(Z)))
    if (length(cells) < 50) next
    term_names <- if (is.null(cond_prefix))
                    paste0(fc, "::", TYPES)
                  else paste0(fc, "::", cond_prefix, ":", TYPES)
    keep <- term_names %in% colnames(Z) & term_names %in% rownames(fit$U)
    if (!any(keep)) next
    tn <- term_names[keep]
    tt <- TYPES[keep]
    ## Sigma = cov of the focal's kernel columns; the attribution is pace::pair_variance_pratt
    ## (which reads non-finite BLUPs as zero, as this function always has).
    Sigma_K <- pace_subset_covariance_cpp(Z, cells, match(tn, colnames(Z)), numeric(0))
    pratt <- pace_pair_variance_pratt_cpp(Sigma_K, as.matrix(fit$U[tn, , drop = FALSE]))

    for (i in seq_along(tt)) {
      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        cohort = cohort_label, block = block_label,
        focal = fc, neighbour = tt[i],
        V_pair_pratt = pratt$V_pair_pratt[i],
        V_pair_diag  = pratt$V_pair_diag[i],
        within_focal_share_pct = pratt$within_focal_share_pct[i],
        within_focal_diag_pct  = pratt$within_focal_diag_pct[i],
        sign = pratt$sign[i],
        stringsAsFactors = FALSE)
    }
    focal_rows[[length(focal_rows) + 1]] <- data.frame(
      cohort = cohort_label, block = block_label, focal = fc,
      n_cells = length(cells),
      V_block_pratt = pratt$V_block_pratt,
      V_block_diag  = pratt$V_block_diag,
      cross_cov_pct = pratt$cross_cov_pct,
      n_negative_pairs = pratt$n_negative_pairs,
      stringsAsFactors = FALSE)
  }
  list(pair_long = do.call(rbind, pair_rows),
       focal_summary = do.call(rbind, focal_rows))
}

## ---------------------------------------------------------------------------
## Per-focal stacked bar + zoom (verbatim engine port).
## ---------------------------------------------------------------------------
plot_stacked_bar_focal_no_resp <- function(
        blks_all,
        block_levels = c("Cell type", "Spatial cell state", "Spillover", "Residuals"),
        title = "") {
  keep_levels <- setdiff(block_levels, "Residuals")

  resid_by_focal <- blks_all %>%
    dplyr::filter(.data$block == "Residuals") %>%
    dplyr::mutate(focal = stringr::str_replace(.data$focal, "_", " ")) %>%
    dplyr::transmute(focal, resid_pct = .data$pct_total, explained = 1 - .data$pct_total)

  blks_plot <- blks_all %>%
    dplyr::filter(.data$block != "Residuals") %>%
    dplyr::mutate(focal = stringr::str_replace(.data$focal, "_", " "),
                  block = factor(.data$block, levels = keep_levels)) %>%
    tidyr::complete(focal, block, fill = list(pct_total = 0)) %>%
    dplyr::left_join(resid_by_focal, by = "focal")

  cols <- c("Cell type" = "#003049", "Spatial cell state" = "#780000",
            "Spillover" = "#E76F51")
  cols <- cols[keep_levels]

  focal_order <- blks_plot %>%
    dplyr::filter(block == "Spatial cell state") %>%
    dplyr::arrange(dplyr::desc(pct_total)) %>%
    dplyr::pull(focal)
  blks_plot$focal <- factor(blks_plot$focal, levels = focal_order)

  p_full <- ggplot(blks_plot, aes(x = focal, y = pct_total, fill = block)) +
    geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
    labs(x = "Focal cell type", y = "Proportion of total variance",
         fill = "Block", title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))

  zoom_data <- blks_plot %>%
    dplyr::filter(block %in% c("Spatial cell state", "Spillover"))

  p_zoom <- ggplot(zoom_data, aes(x = focal, y = pct_total, fill = block)) +
    geom_col(width = 0.7, position = position_dodge(width = 0.75)) +
    # Vertical labels: dodged bars of similar height put their labels side by
    # side at the same y, where horizontal text from neighbouring bars overlaps.
    # Rotated, each label stays within its own bar's width.
    geom_text(aes(label = sprintf("%.2f%%", pct_total * 100)),
              position = position_dodge(width = 0.75), angle = 90,
              hjust = -0.15, vjust = 0.5, size = 2.5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.28))) +
    scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
    labs(x = "Focal cell type", y = NULL, fill = "Block",
         title = "Spatial cell state & Spillover") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none")

  p_full <- p_full + theme(legend.position = "bottom")
  patchwork::wrap_plots(p_full, p_zoom, widths = c(1, 1))
}

## ---------------------------------------------------------------------------
## Per-focal stacked bar + zoom WITH the responder spatial block (verbatim port
## of plot_stacked_bar_focal, for condition cohorts).
## ---------------------------------------------------------------------------
plot_stacked_bar_focal <- function(
        blks_all,
        block_levels = c("Cell type", "Spatial cell state", "Responder spatial state",
                         "Responder status", "Spillover", "Residuals"),
        title = "",
        zoom_blocks = c("Spatial cell state", "Responder spatial state", "Spillover")) {
  keep_levels <- setdiff(block_levels, "Residuals")

  resid_by_focal <- blks_all %>%
    dplyr::filter(.data$block == "Residuals") %>%
    dplyr::mutate(focal = stringr::str_replace(.data$focal, "_", " ")) %>%
    dplyr::transmute(focal, resid_pct = .data$pct_total, explained = 1 - .data$pct_total)

  blks_plot <- blks_all %>%
    dplyr::filter(.data$block != "Residuals") %>%
    dplyr::mutate(focal = stringr::str_replace(.data$focal, "_", " "),
                  block = factor(.data$block, levels = keep_levels)) %>%
    tidyr::complete(focal, block, fill = list(pct_total = 0)) %>%
    dplyr::left_join(resid_by_focal, by = "focal")

  cols <- c("Cell type" = "#003049", "Spatial cell state" = "#780000",
            "Responder spatial state" = "#c1121f", "Responder status" = "#669bbc",
            "Spillover" = "#E5A100")
  cols <- cols[keep_levels]

  focal_order <- blks_plot %>%
    dplyr::filter(block == "Spatial cell state") %>%
    dplyr::arrange(dplyr::desc(pct_total)) %>%
    dplyr::pull(focal)
  blks_plot$focal <- factor(blks_plot$focal, levels = focal_order)

  p_full <- ggplot(blks_plot, aes(x = focal, y = pct_total, fill = block)) +
    geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
    labs(x = "Focal cell type", y = "Proportion of total variance",
         fill = "Block", title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))

  zoom_blocks <- intersect(zoom_blocks, keep_levels)
  zoom_data <- blks_plot %>% dplyr::filter(block %in% zoom_blocks)

  p_zoom <- ggplot(zoom_data, aes(x = focal, y = pct_total, fill = block)) +
    geom_col(width = 0.7, position = position_dodge(width = 0.75)) +
    # Vertical labels: dodged bars of similar height put their labels side by
    # side at the same y, where horizontal text from neighbouring bars overlaps.
    # Rotated, each label stays within its own bar's width.
    geom_text(aes(label = sprintf("%.2f%%", pct_total * 100)),
              position = position_dodge(width = 0.75), angle = 90,
              hjust = -0.15, vjust = 0.5, size = 2.5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.28))) +
    scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
    labs(x = "Focal cell type", y = NULL, fill = "Block",
         title = paste(zoom_blocks, collapse = " & ")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none")

  p_full <- p_full + theme(legend.position = "bottom")
  patchwork::wrap_plots(p_full, p_zoom, widths = c(1, 1))
}

## ---------------------------------------------------------------------------
## Per-pair driver-score bar + per-gene single-frame decomposition (verbatim).
## ---------------------------------------------------------------------------
make_pair_composite <- function(focus_focal, focus_neighbour, mv, sf = NULL, n_top = 5,
                                panels = c("both", "score", "gene")) {
  panels      <- match.arg(panels)
  focus_key   <- paste(focus_focal, focus_neighbour, sep = "_")
  top_drivers <- utils::head(mv$mcsd_canonical[[focus_key]]$scores, n_top)

  p_score <- ggplot2::ggplot(top_drivers, ggplot2::aes(driver_score, stats::reorder(gene, driver_score))) +
    ggplot2::geom_col(fill = "#F5A623", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", driver_score)), hjust = -0.15, size = 3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
    ggplot2::labs(x = "Driver score (spatial cell state)", y = NULL,
                  title = sprintf("%s <- %s", focus_focal, focus_neighbour)) +
    ggplot2::theme_classic(base_size = 11)
  if (panels == "score") return(p_score)

  if (is.null(sf)) stop("`sf` (single-frame decomposition) is required for the gene panel.")
  has_resp    <- "Responder spatial %" %in% names(sf)
  gene_blocks <- c("Cell type", "Spatial cell state",
                   if (has_resp) "Responder spatial state", "Spillover")
  block_cols  <- c("Cell type" = "#1F3B57", "Spatial cell state" = "#8B1A1A",
                   "Responder spatial state" = "#E07B39", "Spillover" = "#F6C9C9")

  gd <- dplyr::filter(sf, .data$focal == focus_focal, .data$gene %in% top_drivers$gene)
  gene_var <- data.frame(gene = gd$gene,
                         `Cell type`          = gd[["Cell type %"]],
                         `Spatial cell state` = gd[["Spatial %"]],
                         `Spillover`          = gd[["Spillover %"]],
                         check.names = FALSE)
  if (has_resp) gene_var[["Responder spatial state"]] <- gd[["Responder spatial %"]]
  gene_var <- gene_var |>
    tidyr::pivot_longer(-gene, names_to = "Block", values_to = "pct") |>
    dplyr::mutate(Block = factor(Block, levels = gene_blocks),
                  gene  = factor(gene, levels = top_drivers$gene))

  p_gene <- ggplot2::ggplot(gene_var, ggplot2::aes(gene, pct / 100, fill = Block)) +
    ggplot2::geom_col(width = 0.78, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_fill_manual(values = block_cols[gene_blocks]) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Gene", y = "Proportion of total variance",
                  title = "Per-gene decomposition (single-frame)") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (panels == "gene") return(p_gene)

  patchwork::wrap_plots(p_score, p_gene)
}

## ---------------------------------------------------------------------------
## Exported, PACEFit-facing wrappers
## ---------------------------------------------------------------------------

#' Per-focal variance decomposition plot
#'
#' Reproduces the manuscript per-focal decomposition figure: a full stacked bar
#' of the cell-type identity, spatial cell state, and contamination (spillover)
#' components, with a zoomed grouped bar of the two small spatial components.
#'
#' @param object A [PACEFit].
#' @param title Plot title.
#' @return A `patchwork` / `ggplot` object.
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' plotDecomposition(fit)
#' @export
plotDecomposition <- function(object, title = "Per-focal decomposition") {
  stopifnot(methods::is(object, "PACEFit"))
  sf <- object@varianceDecomposition$perGene
  fb <- .pace_focal_blocks(sf)
  has_resp <- "Responder spatial state" %in% names(fb)

  if (has_resp) {
    blks_all <- fb |>
      dplyr::transmute(.data$focal,
                       `Cell type`               = .data$`Cell type` / 100,
                       `Spatial cell state`      = .data$`Spatial cell state` / 100,
                       `Responder spatial state` = .data$`Responder spatial state` / 100,
                       `Spillover`               = .data$Spillover / 100,
                       `Residuals`               = .data$Residual / 100) |>
      tidyr::pivot_longer(-.data$focal, names_to = "block", values_to = "pct_total")
    plot_stacked_bar_focal(
      blks_all,
      block_levels = c("Cell type", "Spatial cell state", "Responder spatial state",
                       "Spillover", "Residuals"),
      title = title)
  } else {
    blks_all <- fb |>
      dplyr::transmute(.data$focal,
                       `Cell type`          = .data$`Cell type` / 100,
                       `Spatial cell state` = .data$`Spatial cell state` / 100,
                       `Spillover`          = .data$Spillover / 100,
                       `Residuals`          = .data$Residual / 100) |>
      tidyr::pivot_longer(-.data$focal, names_to = "block", values_to = "pct_total")
    plot_stacked_bar_focal_no_resp(
      blks_all,
      block_levels = c("Cell type", "Spatial cell state", "Spillover", "Residuals"),
      title = title)
  }
}

#' Pairwise spatial variance table
#'
#' Returns the tidy per-(focal, neighbour) spatial variance table that
#' [plotPairHeatmap()] draws: for each focal cell type, the percentage of its
#' total expression variance contributed by spatial interaction with each
#' neighbour, obtained by splitting the focal's spatial share across neighbours
#' by a normalised (off-diagonal-renormalised) Pratt attribution. Off-diagonal
#' pairs only.
#'
#' @param object A [PACEFit].
#' @param block Either `"spatial"` (default) or `"responder"` (condition cohorts).
#' @param ... Unused.
#' @return A data frame with columns `focal`, `neighbour`, and `val` (spatial %
#'   of the focal type's total variance).
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' head(pairVariance(fit))
#' @rdname pairVariance
#' @export
setMethod("pairVariance", "PACEFit",
  function(object, block = c("spatial", "responder"), ...) {
    block <- match.arg(block)
    mv <- .pace_mv_shim(object)
    sf <- object@varianceDecomposition$perGene
    fb <- .pace_focal_blocks(sf)

    if (block == "responder") {
      if (!"Responder spatial state" %in% names(fb))
        stop("this fit has no responder block; use block = 'spatial'.", call. = FALSE)
      cond_prefix <- object@params$resp_term
      foc         <- stats::setNames(fb$`Responder spatial state`, as.character(fb$focal))
      block_label <- "RxS"
    } else {
      cond_prefix <- NULL
      foc         <- stats::setNames(fb$`Spatial cell state`, as.character(fb$focal))
      block_label <- "Spatial"
    }

    pace_pair_variance_pratt(mv, cond_prefix = cond_prefix,
                             cohort_label = "cohort",
                             block_label = block_label)$pair_long |>
      dplyr::filter(as.character(.data$focal) != as.character(.data$neighbour)) |>
      dplyr::group_by(.data$focal) |>
      dplyr::mutate(val = foc[as.character(.data$focal)] *
                          .data$V_pair_pratt / sum(.data$V_pair_pratt)) |>
      dplyr::ungroup() |>
      dplyr::transmute(focal = as.character(.data$focal),
                       neighbour = as.character(.data$neighbour), .data$val)
  })

#' Pairwise spatial-percent heatmap
#'
#' Reproduces the manuscript pairwise heatmap: each cell is the percentage of the
#' focal type's total variance contributed by spatial interaction with that
#' neighbour, computed as the focal's spatial share split across neighbours by the
#' off-diagonal-renormalised Pratt attribution.
#'
#' @param object A [PACEFit].
#' @param block Which spatial block to map: `"spatial"` (default, the baseline
#'   spatial cell state) or `"responder"` (the responder-by-proximity block, for
#'   condition cohorts).
#' @param title Plot title; a sensible default is chosen per `block`.
#' @return A `ggplot` object.
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' plotPairHeatmap(fit)
#' @export
plotPairHeatmap <- function(object, block = c("spatial", "responder"), title = NULL) {
  stopifnot(methods::is(object, "PACEFit"))
  block <- match.arg(block)
  sf    <- object@varianceDecomposition$perGene
  TYPES <- sort(unique(as.character(sf$focal)))

  if (block == "responder") {
    fill_name <- "responder %\nof total"
    if (is.null(title)) title <- "Pairwise responder spatial % (normalised Pratt)"
  } else {
    fill_name <- "spatial %\nof total"
    if (is.null(title)) title <- "Pairwise spatial % (variance, normalised Pratt)"
  }

  pratt_long <- pairVariance(object, block = block)

  heat <- expand.grid(focal = TYPES, neighbour = TYPES, stringsAsFactors = FALSE) |>
    dplyr::filter(.data$focal != .data$neighbour) |>
    dplyr::left_join(pratt_long, by = c("focal", "neighbour")) |>
    dplyr::mutate(val = ifelse(is.na(.data$val), 0, .data$val),
                  focal = factor(.data$focal, levels = TYPES),
                  neighbour = factor(.data$neighbour, levels = rev(TYPES)))

  ggplot2::ggplot(heat, ggplot2::aes(.data$focal, .data$neighbour, fill = .data$val)) +
    ggplot2::geom_tile(colour = "grey92", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(abs(.data$val) < 0.005, "0.00",
                                                   sprintf("%.2f", .data$val))),
                       size = 2.4, colour = "white") +
    ggplot2::scale_fill_viridis_c(option = "inferno", na.value = "white",
                                  name = fill_name) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Focal", y = "Neighbour", title = title) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Per-pair driver composite
#'
#' Reproduces the manuscript per-pair driver figure: the top genes ranked by
#' their driver score for a focal-neighbour pair, alongside their per-gene
#' single-frame variance decomposition. The score is carried in the
#' `driver_score` column of `topDrivers(object)`.
#'
#' @param object A [PACEFit].
#' @param focal,neighbour The focal and neighbour cell types.
#' @param n_top Number of top genes to show.
#' @param panels One of `"both"`, `"score"`, or `"gene"`.
#' @return A `patchwork` / `ggplot` object.
#' @examples
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' plotDrivers(fit, "Macrophage", "Tumour")
#' @export
plotDrivers <- function(object, focal, neighbour, n_top = 5,
                        panels = c("both", "score", "gene")) {
  stopifnot(methods::is(object, "PACEFit"))
  mv <- .pace_mv_shim(object)
  sf <- object@varianceDecomposition$perGene
  make_pair_composite(focal, neighbour, mv, sf = sf, n_top = n_top,
                      panels = match.arg(panels))
}

#' Expression versus number of neighbours (per-bin boxplot)
#'
#' Reproduces the manuscript proximity figure (verbatim style port of
#' `density_bin_boxplot` in the analysis helpers): for each focal cell the number
#' of `neighbour`-type cells within `radius` micrometres is counted and binned,
#' and each gene's raw counts are drawn as a boxplot per bin with the per-bin mean
#' overlaid as a point. Because these genes are zero-inflated the box median sits
#' at zero in most bins, so the mean point carries the trend.
#'
#' @param object A [PACEFit] (supplies the cell-type column used at fit time).
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted.
#' @param genes Character vector of genes (shown as facets).
#' @param focal,neighbour Focal and neighbour cell types.
#' @param radius Neighbour search radius in micrometres (default 30).
#' @param breaks Neighbour-count bin breaks (default `c(0, 2, 4, 6, 8, Inf)`).
#' @param condition Optional `colData` column splitting the boxes by group, so
#'   the two arms of a condition cohort can be compared within each bin. `TRUE`
#'   uses the fit's own `condition_col`; `NULL` (default) pools all cells.
#' @param condition_colours Colours for the condition levels.
#' @param box_colour Box and mean-point colour, used when `condition` is `NULL`.
#' @param assay_name Counts assay name (default `"counts"`).
#' @return A `ggplot` object.
#' @examples
#' spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
#' fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
#' g <- head(subset(neighbourSlopes(fit),
#'                   focal == "Macrophage" & neighbour == "Tumour")$gene, 2)
#' plotProximity(fit, spe, g, "Macrophage", "Tumour")
#' @importFrom dbscan frNN
#' @export
plotProximity <- function(object, spe, genes, focal, neighbour,
                          radius = 30, breaks = c(0, 2, 4, 6, 8, Inf),
                          condition = NULL, condition_colours = c("#3B6FB6", "#C0392B"),
                          box_colour = "#4F8B5E", assay_name = "counts") {
  stopifnot(methods::is(object, "PACEFit"))
  if (isTRUE(condition)) condition <- object@params$condition_col
  if (!is.null(condition) && !condition %in%
      names(SummarizedExperiment::colData(spe)))
    stop("`condition` column '", condition, "' is not in colData(spe).",
         call. = FALSE)
  ct     <- as.character(SummarizedExperiment::colData(spe)[[object@params$celltype_col]])
  coords <- as.matrix(SpatialExperiment::spatialCoords(spe))
  focal_idx <- which(ct == focal)
  nbr_idx   <- which(ct == neighbour)

  ## number of neighbour-type cells within `radius` um of each focal cell
  fr <- dbscan::frNN(x = coords[nbr_idx, , drop = FALSE], eps = radius,
                     query = coords[focal_idx, , drop = FALSE])
  n_neighbours <- lengths(fr$id)

  counts   <- SummarizedExperiment::assay(spe, assay_name)
  expr_mat <- as.matrix(counts[genes, focal_idx, drop = FALSE])   # genes x cells (RAW)
  plot_data <- data.frame(
    n_neighbours = rep(n_neighbours, times = length(genes)),
    gene = factor(rep(genes, each = length(focal_idx)), levels = genes),
    expr = as.numeric(t(expr_mat)))
  if (!is.null(condition)) {
    arm <- as.character(SummarizedExperiment::colData(spe)[[condition]])[focal_idx]
    plot_data$condition <- factor(rep(arm, times = length(genes)))
  }
  plot_data$bin <- cut(plot_data$n_neighbours, breaks = breaks, include.lowest = TRUE)
  plot_data <- plot_data[!is.na(plot_data$bin), ]
  ## relabel the open top bin with the observed maximum count
  lev <- levels(plot_data$bin)
  lev[length(lev)] <- sprintf("(%g,%d]", breaks[length(breaks) - 1L],
                              max(plot_data$n_neighbours))
  levels(plot_data$bin) <- lev

  if (is.null(condition)) {
    p <- ggplot(plot_data, aes(.data$bin, .data$expr)) +
      geom_boxplot(colour = box_colour, fill = "white", linewidth = 0.6,
                   outlier.colour = box_colour, outlier.alpha = 0.25,
                   outlier.size = 1.2) +
      stat_summary(fun = mean, geom = "point", colour = box_colour, size = 3.2)
  } else {
    ## One box per (bin, arm); the arms are the comparison, so they are dodged
    ## side by side rather than faceted apart.
    cols <- condition_colours
    if (is.null(names(cols)))
      names(cols) <- levels(plot_data$condition)[seq_along(cols)]
    p <- ggplot(plot_data, aes(.data$bin, .data$expr, colour = .data$condition)) +
      geom_boxplot(fill = "white", linewidth = 0.6, outlier.alpha = 0.25,
                   outlier.size = 1.2,
                   position = position_dodge(width = 0.8)) +
      stat_summary(fun = mean, geom = "point", size = 3.2,
                   position = position_dodge(width = 0.8)) +
      scale_colour_manual(values = cols, name = condition)
  }
  p <- p +
    labs(x = paste("Number of", neighbour, "neighbours"),
         y = paste(focal, "RAW counts")) +
    theme_minimal(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

  if (length(genes) == 1L)
    p + labs(title = paste(genes, "expression"),
             y = paste(focal, genes, "counts"))
  else
    p + facet_wrap(~ gene, scales = "free_y")
}

## ---------------------------------------------------------------------------
## Condition-stratified response figures.
##
## Ported from pace_mv/scripts/helpers/pace_response_figures.R, the code behind
## the manuscript's melanoma figures, so the package draws them rather than
## leaving each user to reassemble them.
## ---------------------------------------------------------------------------

## PACE shrunken slope of `gene` for focal <- neighbour, per condition arm, read
## from the fit rather than refitted. The reference arm is the focal::neighbour
## BLUP; the other arm adds the focal::<resp_term>:neighbour interaction.
.pace_arm_slopes <- function(object, gene, focal, neighbour) {
  U <- object@fit$U
  if (!gene %in% colnames(U)) stop("gene not in the fit: ", gene, call. = FALSE)
  u <- stats::setNames(U[, gene], rownames(U))
  base <- paste0(focal, "::", neighbour)
  if (!base %in% names(u))
    stop("no slope for ", focal, " <- ", neighbour, call. = FALSE)
  s_ref <- unname(u[[base]])
  s_alt <- s_ref
  rt <- object@params$resp_term
  if (!is.null(rt)) {
    int <- paste0(focal, "::", rt, ":", neighbour)
    if (int %in% names(u)) s_alt <- s_ref + unname(u[[int]])
  }
  c(ref = s_ref, alt = s_alt)
}

## Focal-cell kernel density of `neighbour` and log1p CP10K of `gene`, plus the
## arm label. The density is the Gaussian kernel the model itself used, which
## the fit already carries in its working frame.
.pace_response_data <- function(object, spe, gene, focal, neighbour) {
  df <- object@context$df
  if (!neighbour %in% names(df))
    stop("no kernel column for neighbour '", neighbour, "' in the fit.",
         call. = FALSE)
  cond <- object@params$condition_col
  if (is.null(cond))
    stop("this fit has no condition; use plotProximity() instead.", call. = FALSE)
  Y <- SummarizedExperiment::assay(spe, object@params$assay_name)
  if (!gene %in% rownames(Y)) stop("gene not in `spe`: ", gene, call. = FALSE)
  ## The fit's own kernel column is centred within (image, celltype), which puts
  ## half the density below zero. Recompute the RAW kernel with the fit's own
  ## settings: a slope is unchanged by centring, so it still overlays correctly.
  types <- object@cellTypes
  K <- pace_neighbour_kernel(
    coords = as.matrix(df[, c("x", "y")]), celltype = as.character(df$celltype),
    types = types, h_bio = object@params$h_bio, h_tech = object@params$h_tech,
    eps = 3 * object@params$h_bio, image = as.character(df$imageID),
    per_image = isTRUE(object@params$kernel_per_image))$K_bio

  keep <- as.character(df$celltype) == focal
  cp10k <- as.numeric(Y[gene, keep]) / pmax(df$nCount[keep], 1) * 1e4
  data.frame(density = as.numeric(K[keep, neighbour]),
             expr    = log1p(cp10k),
             arm     = factor(as.character(df[[cond]])[keep]),
             x       = df$x[keep], y = df$y[keep],
             image   = as.character(df$imageID)[keep])
}

#' Response curve: binned means with the PACE slope overlaid
#'
#' Population view of one gene in `focal` cells against `neighbour` density,
#' split by the fit's condition. The solid line and ribbon are binned means and
#' standard errors of the observed data; the dashed line is the PACE slope read
#' straight from the fitted model, anchored at each arm's data centroid. It is
#' the model's estimate drawn over the data, not a smoother refitted to it.
#'
#' @param object A [PACEFit] fitted with a `condition_col`.
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted.
#' @param gene,focal,neighbour Gene, focal cell type and neighbour cell type.
#' @param n_bins Number of density quantile bins (default 10).
#' @param min_bin_n Bins with fewer focal cells than this are dropped.
#' @param colours Colours for the two arms.
#' @return A `ggplot` object.
#' @examples
#' spe <- readRDS(system.file("extdata", "mel_cosmx_subset.rds", package = "PACE"))
#' \donttest{
#' fit <- paceFit(spe, celltype_col = "cellType", condition_col = "Responder",
#'                image_col = "image", kernel_per_image = TRUE,
#'                image_re = "intercept", verbose = FALSE)
#' plotResponseCurve(fit, spe, "SPP1", "Macrophage", "Tumour")
#' }
#' @export
plotResponseCurve <- function(object, spe, gene, focal, neighbour,
                              n_bins = 10, min_bin_n = 10,
                              colours = c("#3B6FB6", "#C0392B")) {
  stopifnot(methods::is(object, "PACEFit"))
  d <- .pace_response_data(object, spe, gene, focal, neighbour)
  d <- d[is.finite(d$density) & is.finite(d$expr), ]

  ## Which arm carries the interaction comes from `resp_term`, recorded when the
  ## model was fitted, NOT from the level order here. `.pace_response_data()`
  ## rebuilds the arm with factor(as.character(...)), which sorts alphabetically
  ## and throws away the order `model.matrix()` actually used, so the two agree
  ## only by coincidence. They do coincide on the shipped melanoma subset, whose
  ## reference "nonPD" is also the alphabetically first of its two levels, which
  ## is why this went unseen; on a cohort whose arms are "R" and "NR" the model's
  ## reference is R while the alphabetically first is NR, and both dashed slopes
  ## were then drawn on, and labelled with, the wrong arm.
  cond      <- object@params$condition_col
  resp_term <- object@params$resp_term
  if (is.null(resp_term) || !startsWith(resp_term, cond)) {
    stop("this fit records resp_term '", resp_term, "', which does not begin with ",
         "condition_col '", cond, "', so the arm carrying the interaction cannot ",
         "be identified.", call. = FALSE)
  }
  ## substring() rather than sub(), so a condition column whose name contains a
  ## regular-expression metacharacter is not silently mismatched.
  case_arm <- substring(resp_term, nchar(cond) + 1L)
  if (!case_arm %in% levels(d$arm)) {
    stop("resp_term '", resp_term, "' names the level '", case_arm, "', which is ",
         "not among the ", cond, " levels present here (",
         paste(levels(d$arm), collapse = ", "), ").", call. = FALSE)
  }
  ref_arm <- setdiff(levels(d$arm), case_arm)[1]
  slopes  <- .pace_arm_slopes(object, gene, focal, neighbour)
  slope_by_arm <- stats::setNames(c(slopes[["ref"]], slopes[["alt"]]),
                                  c(ref_arm, case_arm))
  d$arm <- factor(as.character(d$arm), levels = c(case_arm, ref_arm))
  arms  <- levels(d$arm)

  brks <- unique(stats::quantile(d$density, seq(0, 1, length.out = n_bins + 1),
                                 na.rm = TRUE))
  d$bin <- cut(d$density, breaks = brks, include.lowest = TRUE)
  binned <- d |>
    dplyr::group_by(.data$arm, .data$bin) |>
    dplyr::summarise(x = mean(.data$density), mean_e = mean(.data$expr),
                     se = stats::sd(.data$expr) / sqrt(dplyr::n()),
                     n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(.data$n >= min_bin_n)

  x_range <- c(0, max(d$density))
  centro <- d |>
    dplyr::group_by(.data$arm) |>
    dplyr::summarise(dx = mean(.data$density), ey = mean(.data$expr),
                     .groups = "drop")
  pace_lines <- do.call(rbind, lapply(seq_len(nrow(centro)), function(i) {
    a <- as.character(centro$arm[i])
    data.frame(arm = centro$arm[i], x = x_range,
               y = centro$ey[i] + slope_by_arm[[a]] * (x_range - centro$dx[i]))
  }))

  cols <- colours
  if (is.null(names(cols))) names(cols) <- c(ref_arm, case_arm)

  ggplot() +
    geom_ribbon(data = binned,
                aes(.data$x, ymin = .data$mean_e - .data$se,
                    ymax = .data$mean_e + .data$se, fill = .data$arm),
                alpha = 0.18) +
    geom_line(data = binned, aes(.data$x, .data$mean_e, colour = .data$arm),
              linewidth = 0.9) +
    geom_point(data = binned, aes(.data$x, .data$mean_e, colour = .data$arm),
               size = 1.8) +
    geom_line(data = pace_lines, aes(.data$x, .data$y, colour = .data$arm),
              linetype = "dashed", linewidth = 1.1) +
    scale_colour_manual(values = cols, name = object@params$condition_col) +
    scale_fill_manual(values = cols, guide = "none") +
    coord_cartesian(xlim = x_range, ylim = c(0, max(binned$mean_e + binned$se) * 1.05)) +
    labs(title = sprintf("%s in %s vs %s-neighbour density", gene, focal, neighbour),
         subtitle = sprintf("Solid = binned means +/- SE; dashed = PACE slope (%s)",
                            paste(sprintf("%s %+.3f", arms,
                                          slope_by_arm[arms]),
                                  collapse = ", ")),
         x = sprintf("%s-neighbour density (Gaussian kernel)", neighbour),
         y = sprintf("%s expression (CP10K, log1p)", gene)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9))
}

#' Tissue maps and scatters for one exemplar image per arm
#'
#' The spatial view behind [plotResponseCurve()]. For one image per condition
#' arm: a smooth kernel-density map of the neighbour type with the focal cells
#' overlaid and coloured by expression (top), and expression against neighbour
#' density for the same cells with a fitted line (bottom). Colour scales are
#' shared across the panels so the arms are directly comparable.
#'
#' @param object A [PACEFit] fitted with a `condition_col`.
#' @param spe The [SpatialExperiment::SpatialExperiment] that was fitted.
#' @param gene,focal,neighbour Gene, focal cell type and neighbour cell type.
#' @param images Optional named character vector of image identifiers, one per
#'   arm. The default picks, for each arm, the image with the most focal cells.
#' @param kde_h,ngrid Kernel bandwidth and grid size for the density map.
#' @param scalebar_len Scale-bar length in micrometres.
#' @return A `patchwork` object.
#' @examples
#' spe <- readRDS(system.file("extdata", "mel_cosmx_subset.rds", package = "PACE"))
#' \donttest{
#' fit <- paceFit(spe, celltype_col = "cellType", condition_col = "Responder",
#'                image_col = "image", kernel_per_image = TRUE,
#'                image_re = "intercept", verbose = FALSE)
#' plotResponseMap(fit, spe, "SPP1", "Macrophage", "Tumour")
#' }
#' @export
plotResponseMap <- function(object, spe, gene, focal, neighbour, images = NULL,
                            kde_h = 60, ngrid = 200, scalebar_len = 100) {
  stopifnot(methods::is(object, "PACEFit"))
  df   <- object@context$df
  cond <- object@params$condition_col
  if (is.null(cond))
    stop("this fit has no condition; use plotProximity() instead.", call. = FALSE)
  fd <- .pace_response_data(object, spe, gene, focal, neighbour)

  ## One exemplar image per arm: by default the one with the most focal cells,
  ## since a map of a handful of cells shows nothing.
  ref_arm  <- levels(fd$arm)[1]
  case_arm <- setdiff(levels(fd$arm), ref_arm)[1]
  if (is.null(images)) {
    images <- vapply(split(fd, fd$arm), function(a) {
      tb <- sort(table(a$image), decreasing = TRUE); names(tb)[1]
    }, character(1))
  }
  ## Case arm on the left, as in the manuscript figure.
  ord <- c(case_arm, ref_arm)
  ord <- ord[ord %in% names(images)]
  if (length(ord) == length(images)) images <- images[ord]
  arm_cols <- stats::setNames(c("#C0392B", "#3B6FB6"), c(case_arm, ref_arm))

  ## Shared scales across arms, or the two panels cannot be compared. Ranged on
  ## the sections actually drawn, not the whole cohort, which would leave most
  ## of each panel empty.
  shown       <- fd[fd$image %in% images, , drop = FALSE]
  expr_limits <- range(shown$expr, na.rm = TRUE)
  x_limits    <- range(shown$density, na.rm = TRUE)

  kdes <- lapply(images, function(img) {
    cells <- df[as.character(df$imageID) == img, , drop = FALSE]
    nb <- cells[as.character(cells$celltype) == neighbour, , drop = FALSE]
    xr <- range(cells$x) + c(-20, 20); yr <- range(cells$y) + c(-20, 20)
    k  <- MASS::kde2d(nb$x, nb$y, n = ngrid, h = c(kde_h, kde_h), lims = c(xr, yr))
    g  <- expand.grid(x = k$x, y = k$y)
    g$z <- as.vector(k$z) * nrow(nb)
    g
  })
  density_limits <- range(unlist(lapply(kdes, `[[`, "z")))

  heat_cols <- c("#000000", "#4D0000", "#8B0000", "#E41A1C", "#FF7F00",
                 "#FFD92F", "#FFFFB2")
  emax <- expr_limits[2]
  heat_values <- scales::rescale(c(0, 2, 3.2, 4.2, 5.2, 6.2, emax),
                                 from = c(0, emax))

  maps <- scatters <- vector("list", length(images))
  for (i in seq_along(images)) {
    arm <- names(images)[i]
    fc  <- fd[fd$image == images[i], , drop = FALSE]
    kde <- kdes[[i]]
    sb_x0 <- min(kde$x) + 0.06 * diff(range(kde$x))
    sb_y0 <- min(kde$y) + 0.06 * diff(range(kde$y))

    maps[[i]] <- ggplot() +
      geom_raster(data = kde, aes(.data$x, .data$y, fill = .data$z),
                  interpolate = TRUE) +
      geom_point(data = fc, aes(.data$x, .data$y, colour = .data$expr),
                 size = 2.4, alpha = 0.97, shape = 16) +
      scale_fill_gradientn(
        colours = c("white", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B"),
        limits = density_limits, name = paste(neighbour, "density"),
        guide = guide_colourbar(order = 2, direction = "horizontal",
                                title.position = "top", barwidth = 8,
                                barheight = 0.5)) +
      scale_colour_gradientn(
        colours = heat_cols, values = heat_values, limits = expr_limits,
        name = sprintf("%s %s (log CP10K)", focal, gene),
        guide = guide_colourbar(order = 1, direction = "horizontal",
                                title.position = "top", barwidth = 8,
                                barheight = 0.5)) +
      annotate("segment", x = sb_x0, xend = sb_x0 + scalebar_len,
               y = sb_y0, yend = sb_y0, colour = "black", linewidth = 1.4) +
      annotate("text", x = sb_x0 + scalebar_len / 2,
               y = sb_y0 + 0.035 * diff(range(kde$y)),
               label = paste(scalebar_len, "um"), size = 2.8, vjust = 0) +
      coord_equal(expand = FALSE) +
      labs(title = arm) +
      theme_void(base_size = 11) +
      theme(legend.position = "bottom", legend.box = "horizontal",
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            plot.title = element_text(size = 11, hjust = 0.5,
                                      margin = margin(b = 4)))

    f  <- stats::lm(expr ~ density, data = fc)
    co <- summary(f)$coefficients
    ann <- sprintf("slope = %+.3f\nn = %d", co[2, 1], nrow(fc))
    scatters[[i]] <- ggplot(fc, aes(.data$density, .data$expr)) +
      geom_point(colour = arm_cols[[arm]], alpha = 0.5, size = 1.1) +
      geom_smooth(method = "lm", formula = y ~ x, colour = "black",
                  fill = "grey60", alpha = 0.3) +
      annotate("label", x = x_limits[1], y = expr_limits[2], hjust = 0, vjust = 1,
               label = ann, size = 2.9,
               fill = scales::alpha("white", 0.7)) +
      coord_cartesian(xlim = x_limits, ylim = expr_limits) +
      labs(title = arm, x = sprintf("%s density (Gaussian kernel)", neighbour),
           y = sprintf("%s %s (log CP10K)", focal, gene)) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(size = 10))
  }

  ## wrap_plots rather than Reduce(`+`): patchwork is imported, not attached,
  ## so `+` on ggplots would not dispatch to it here.
  patchwork::wrap_plots(
    patchwork::wrap_plots(maps, nrow = 1),
    patchwork::wrap_plots(scatters, nrow = 1),
    nrow = 2, heights = c(1.4, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
}
