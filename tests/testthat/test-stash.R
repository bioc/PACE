fit <- readRDS(system.file("extdata", "pace_fit_example.rds", package = "PACE"))
spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))

test_that("stashing copies the readouts and leaves the solver state behind", {
  out <- stashPACE(spe, fit)
  stashed <- S4Vectors::metadata(out)$PACE
  expect_identical(stashed$neighbourSlopes, neighbourSlopes(fit))
  expect_identical(stashed$topDrivers, topDrivers(fit))
  expect_identical(stashed$cellTypes, cellTypes(fit))
  # the point of stashing readouts rather than the fit
  expect_false("fit" %in% names(stashed))
  expect_lt(utils::object.size(stashed), utils::object.size(fit))
})

test_that("stashing records which fit produced the readouts", {
  stashed <- S4Vectors::metadata(stashPACE(spe, fit))$PACE
  expect_identical(stashed$params, paceParams(fit))
  expect_identical(stashed$packageVersion, as.character(utils::packageVersion("PACE")))
  expect_s3_class(stashed$stashedAt, "POSIXct")
})

test_that("stashing refuses a SpatialExperiment the fit was not built on", {
  expect_error(stashPACE(spe[, seq_len(10)], fit), "10 cells but the fit has")
})

test_that("stashing will not silently clobber", {
  out <- stashPACE(spe, fit)
  expect_error(stashPACE(out, fit), "already exists")
  expect_silent(again <- stashPACE(out, fit, overwrite = TRUE))
  expect_identical(S4Vectors::metadata(again)$PACE$cellTypes, cellTypes(fit))
})

test_that("stashing takes a name, so several fits can sit side by side", {
  out <- stashPACE(stashPACE(spe, fit), fit, name = "second")
  expect_setequal(names(S4Vectors::metadata(out)), c("PACE", "second"))
})

test_that("stashing refuses anything that is not a fit", {
  expect_error(stashPACE(spe, list()), "must be a PACEFit")
})
