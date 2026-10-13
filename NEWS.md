# PACE 0.99.13

- `verbose = FALSE` is now silent. The shrinkage step reported its progress
  regardless, because it never received the flag, and `mashr` printed its own
  on top of that.
- Installation instructions use `BiocManager::install()`.
- The vignette evaluates every chunk.
- Minimum R version is 4.6.0.

# PACE 0.99.12

- `plotDecomposition()` draws the zoom panel's percentage labels upright, so
  that neighbouring bars of similar height no longer overlap each other's text.

# PACE 0.99.11

- Added a canonical multivariate adaptive shrinkage implementation to the
  compiled core, alongside the existing `mashr` path. It reproduces `mashr`'s
  calls exactly on both reference cohorts.

# PACE 0.99.10

- The driver score columns are named for what the method calls them:
  `driver_score` and `driver_score_spec4`, previously `MCSD` and `MCSD4`.
  `plotDrivers()` takes `panels = "score"` in place of `"mcsd"`.

# PACE 0.99.9

- `topDrivers()` documents which model term each column holds. On a fit with a
  condition the table carries two different slopes, and they are not a raw and
  shrunken pair of the same quantity.
- The condition vignette shows the responder decomposition and driver figures,
  and prints the driver table.

# PACE 0.99.6

- The per-pair gene score is called the driver score throughout the
  documentation and plot labels, matching the method description.

# PACE 0.99.4

- Fixed portability defects found by the Bioconductor build matrix: bit
  identity on arm64, a crash on Windows from `static thread_local` destructors,
  an unsupported `system2(env=)` on Windows, and missing `-pthread` on compile
  flags for Wasm.

# PACE 0.99.0

- Initial Bioconductor submission.
