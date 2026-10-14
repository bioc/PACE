# The core carries its own trigamma because a Python caller has no Rf_trigamma.
# A fit must not depend on which language called it, so the two have to agree.
test_that("the core's trigamma agrees with R's", {
  x <- c(0.001, 0.01, 0.1, 0.5, 1, 1.5, 2, 3, 5, 6, 10,
         19.9, 20, 20.1,            # the recurrence / asymptotic boundary
         50, 100, 1000, 1e5, 1e8)
  expect_equal(PACE:::pace_trigamma_cpp(x), trigamma(x), tolerance = 1e-14)
})

test_that("trigamma handles the awkward arguments the way R does", {
  expect_true(is.nan(PACE:::pace_trigamma_cpp(NaN)))
  # non-positive integers are poles
  expect_equal(PACE:::pace_trigamma_cpp(c(0, -1, -2)), rep(Inf, 3))
  # and the reflection formula holds off the poles
  expect_equal(PACE:::pace_trigamma_cpp(c(-0.5, -1.5, -2.5)),
               trigamma(c(-0.5, -1.5, -2.5)), tolerance = 1e-13)
})
