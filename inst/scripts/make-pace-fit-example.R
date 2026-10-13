## Rebuilds inst/extdata/pace_fit_example.rds, the small fitted object used by
## the man-page examples and by tests/testthat/test-pace-mu.R.
##
## Settings match the object it replaces: no condition, no image random effect,
## a single image, NB1 dispersion and per-cell contamination. The n x G matrices
## are stripped before saving because `mu` is rebuilt on demand, not stored.
suppressPackageStartupMessages({
  library(PACE); library(SpatialExperiment)
})

spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
fit <- paceFit(spe,
               celltype_col     = "cellType",
               kernel_per_image = FALSE,
               image_re         = "none",
               dispersion       = "nb1",
               verbose          = FALSE)

## drop the large per-cell matrices; `mu` is rebuilt by the accessors
for (slot_name in c("mu", "mu_bio", "mu_spill",
                    "technical_offset_mat", "bleed_offset_mat", "stats")) {
  fit@fit[[slot_name]] <- NULL
}

saveRDS(fit, "inst/extdata/pace_fit_example.rds", compress = "xz")
message("wrote inst/extdata/pace_fit_example.rds  (",
        round(file.size("inst/extdata/pace_fit_example.rds") / 1024), " KB)")
