test_that("the packaged fit is valid and its getters return the slots", {
  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  expect_true(validObject(fit))
  expect_identical(cellTypes(fit), fit@cellTypes)
  expect_identical(paceParams(fit), fit@params)
  expect_identical(fitState(fit), fit@fit)
  expect_gt(length(cellTypes(fit)), 1)
  expect_true(all(c("B", "U", "alpha") %in% names(fitState(fit))))
})

test_that("validity refuses a fit whose solver state has been replaced", {
  # the reviewer's example: this passed silently before, then broke downstream
  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  broken <- fit
  broken@fit <- list()
  expect_error(validObject(broken), "missing the solver state")
})

test_that("validity refuses a malformed cell type vector", {
  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  empty <- fit; empty@cellTypes <- character(0)
  expect_error(validObject(empty), "at least one cell type")
  dup <- fit; dup@cellTypes <- rep(cellTypes(fit)[1], 2)
  expect_error(validObject(dup), "unique")
})

test_that("the setters validate what they are given", {
  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  expect_error(neighbourSlopes(fit) <- data.frame(x = 1), "missing columns")
  expect_error(varianceDecomposition(fit) <- list(nope = 1), "must contain .perGene")
  expect_error(topDrivers(fit) <- list(a = list()), "needs a .scores. table")
})

test_that("the setters accept what the pipeline produces", {
  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  round_trip <- fit
  neighbourSlopes(round_trip) <- neighbourSlopes(fit)
  varianceDecomposition(round_trip) <- varianceDecomposition(fit)
  topDrivers(round_trip) <- topDrivers(fit)
  expect_identical(neighbourSlopes(round_trip), neighbourSlopes(fit))
  expect_identical(topDrivers(round_trip), topDrivers(fit))
  expect_true(validObject(round_trip))
})

test_that("an empty reporting layer is valid, since its stage has not run", {
  fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
  bare <- fit
  bare@neighbourSlopes <- data.frame()
  bare@varianceDecomposition <- list()
  bare@topDrivers <- list()
  expect_true(validObject(bare))
})
