// pybind11 bindings for the PACE core.
//
// The core is R-free by construction, so the same translation units back both
// the R package and this module. Arrays cross the boundary as numpy arrays in
// Fortran order, which is what the core expects and what R already uses, so
// neither side copies.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "mash.hpp"

namespace py = pybind11;

namespace {

using ColumnMajor = py::array_t<double, py::array::f_style | py::array::forcecast>;

void raise_if_failed(const pace::Status& status, const char* what) {
  if (status.is_ok()) return;
  throw std::runtime_error(std::string(what) + ": " + status.message);
}

pace::Span<const double> span_of(const ColumnMajor& array) {
  return pace::Span<const double>(array.data(), static_cast<std::int64_t>(array.size()));
}

pace::Span<double> out_span_of(ColumnMajor& array) {
  return pace::Span<double>(array.mutable_data(), static_cast<std::int64_t>(array.size()));
}

// Canonical mash for one slice, end to end. `bhat` and `shat` are
// n_effects x n_conditions; `v` may be empty for the identity, which is what
// PACE publishes with.
py::dict mash(const ColumnMajor& bhat, const ColumnMajor& shat, const ColumnMajor& v,
              double null_weight, double grid_mult, std::int64_t max_iterations,
              double tolerance, const std::vector<std::string>& condition_names) {
  if (bhat.ndim() != 2 || shat.ndim() != 2) throw std::runtime_error("bhat and shat must be 2-D");
  if (bhat.shape(0) != shat.shape(0) || bhat.shape(1) != shat.shape(1)) {
    throw std::runtime_error("bhat and shat must have the same shape");
  }
  const std::int64_t n_effects = bhat.shape(0);
  const std::int64_t n_conditions = bhat.shape(1);

  const std::vector<pace::CovarianceComponent> covariances =
      pace::canonical_covariances(n_conditions, condition_names);

  std::vector<double> grid;
  raise_if_failed(pace::autoselect_grid(span_of(bhat), span_of(shat), n_effects, n_conditions,
                                        grid_mult, &grid),
                  "grid");

  const std::int64_t n_components = pace::n_mixture_components(
      static_cast<std::int64_t>(covariances.size()), static_cast<std::int64_t>(grid.size()), true);

  ColumnMajor loglik({n_effects, n_components});
  raise_if_failed(pace::mixture_loglik(span_of(bhat), span_of(shat), span_of(v), n_effects,
                                       n_conditions, covariances, grid, true, 1,
                                       out_span_of(loglik)),
                  "likelihood");

  std::vector<double> weights;
  std::int64_t iterations = 0;
  raise_if_failed(pace::optimize_mixture_weights(span_of(loglik), n_effects, n_components,
                                                 null_weight, max_iterations, tolerance, &weights,
                                                 &iterations),
                  "mixture weights");

  ColumnMajor posterior_mean({n_effects, n_conditions});
  ColumnMajor posterior_sd({n_effects, n_conditions});
  ColumnMajor lfsr({n_effects, n_conditions});
  raise_if_failed(pace::posterior_summaries(span_of(bhat), span_of(shat), span_of(v), n_effects,
                                            n_conditions, covariances, grid, true, weights,
                                            span_of(loglik), out_span_of(posterior_mean),
                                            out_span_of(posterior_sd), out_span_of(lfsr)),
                  "posteriors");

  py::dict out;
  out["posterior_mean"] = posterior_mean;
  out["posterior_sd"] = posterior_sd;
  out["lfsr"] = lfsr;
  out["pi"] = weights;
  out["grid"] = grid;
  out["iterations"] = iterations;
  return out;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "PACE compiled core";
  m.def("mash", &mash, py::arg("bhat"), py::arg("shat"),
        py::arg("v") = ColumnMajor({0, 0}), py::arg("null_weight") = 10.0,
        py::arg("grid_mult") = 1.4142135623730951, py::arg("max_iterations") = 200000,
        py::arg("tolerance") = 1e-13,
        py::arg("condition_names") = std::vector<std::string>(),
        "Canonical multivariate adaptive shrinkage for one slice.");
}
