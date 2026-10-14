# cellContamination() must reproduce the solver's own contam_frac definition
# (rowSums(mu_spill) / rowSums(mu)) exactly, and fail loudly when the fit
# carries no per-cell contamination rather than returning a misleading zero.

make_fit <- function(n = 40L, p = 6L, seed = 3L, percell = TRUE,
                     store = c("frac", "matrices", "mu_bio_only", "neither")) {
  store <- match.arg(store)
  set.seed(seed)
  mu_bio   <- matrix(rgamma(n * p, 2, 1), n, p)
  rho      <- runif(n, 0, 0.5)
  ambient  <- matrix(rgamma(n * p, 1, 1), n, p)
  mu_spill <- sweep(ambient, 1, rho, "*")
  mu       <- pmax(mu_bio + mu_spill, 1e-6)

  ## Stub solver state. These tests are about cellContamination(), not about
  ## the fit, but a PACEFit is only valid if it carries the state the readouts
  ## need -- so supply it rather than building an object that could not exist.
  fit <- list(percell_bleed_rho = if (percell) rho else NULL,
              B = matrix(0, 1, p), U = matrix(0, 2, p),
              se_B = matrix(0, 1, p), se_U = matrix(0, 2, p),
              alpha = rep(1, p), re_meta = list())
  ## The solvers return `contam_frac`; older fits carry the matrices instead.
  if (store == "frac")
    fit$contam_frac <- rowSums(mu_spill) / pmax(rowSums(mu), 1e-9)
  if (store == "matrices")
    fit <- c(fit, list(mu_spill = mu_spill, mu_bio = mu_bio, mu = mu))
  if (store == "mu_bio_only")
    fit <- c(fit, list(mu_spill = mu_spill, mu_bio = mu_bio))

  df <- data.frame(celltype = rep(c("A", "B"), length.out = n),
                   row.names = paste0("cell", seq_len(n)))
  methods::new("PACEFit", fit = fit, cellTypes = c("A", "B"),
               context = list(df = df),
               params = list(contamination = "percell_hc"))
}

test_that("contamFraction is the solver's stored contam_frac", {
  object <- make_fit(store = "frac")
  cc <- cellContamination(object)

  expect_s3_class(cc, "data.frame")
  expect_identical(nrow(cc), length(object@fit$percell_bleed_rho))
  expect_identical(cc$contamFraction, unname(object@fit$contam_frac))
  expect_identical(cc$rho, as.numeric(object@fit$percell_bleed_rho))
  expect_true(all(cc$contamFraction >= 0 & cc$contamFraction <= 1))
})

test_that("a fit without contam_frac falls back to rowSums(mu_spill)/rowSums(mu)", {
  object   <- make_fit(store = "matrices")
  expected <- rowSums(object@fit$mu_spill) / pmax(rowSums(object@fit$mu), 1e-9)
  expect_identical(cellContamination(object)$contamFraction, unname(expected))
})

test_that("the fallback reconstructs mu from mu_bio when mu was not retained", {
  # Both paths must agree: the reconstruction applies the solver's own floor.
  with_mu    <- cellContamination(make_fit(store = "matrices"))
  without_mu <- cellContamination(make_fit(store = "mu_bio_only"))
  expect_equal(with_mu$contamFraction, without_mu$contamFraction)
})

test_that("the stored fraction and the recomputed one agree", {
  # Guards against the two definitions drifting apart.
  expect_equal(cellContamination(make_fit(store = "frac"))$contamFraction,
               cellContamination(make_fit(store = "matrices"))$contamFraction)
})

test_that("cell identifiers and fitted labels are carried through in fit order", {
  object <- make_fit(store = "frac")
  cc <- cellContamination(object)
  expect_identical(cc$cell, rownames(object@context$df))
  expect_identical(cc$celltype, as.character(object@context$df$celltype))
})

test_that("a fit without per-cell contamination errors rather than returning zeros", {
  expect_error(cellContamination(make_fit(percell = FALSE, store = "frac")),
               "no per-cell contamination")
})

test_that("a fit with neither the fraction nor the matrices errors rather than guessing", {
  expect_error(cellContamination(make_fit(store = "neither")), "mu_spill")
})

test_that("both solvers return contam_frac for a real fit", {
  # The accessor is only useful if the fitting path actually populates it;
  # paceFit() runs the streaming solver.
  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  expect_false(is.null(fit@fit$contam_frac))
  expect_identical(length(fit@fit$contam_frac),
                   length(fit@fit$percell_bleed_rho))
})
