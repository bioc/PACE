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

#include "decomposition.hpp"
#include "neighbourhood.hpp"
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

// Driver scores for one focal-neighbour pair. Pure arithmetic over per-gene
// vectors, so the same call backs `topDrivers()` in R. `u_raw` is only read on
// a condition fit, where V_S comes from the unshrunken baseline BLUP and V_RxS
// from the shrunken interaction; pass it empty otherwise.
py::dict driver_scores(const ColumnMajor& estimate_shrunk, const ColumnMajor& spec,
                       const ColumnMajor& focal_mean, const ColumnMajor& mu_bar,
                       const ColumnMajor& alpha, const ColumnMajor& u_raw, double var_n,
                       double var_rn, bool has_responder) {
  const std::int64_t n = estimate_shrunk.size();
  for (const auto* a : {&spec, &focal_mean, &mu_bar, &alpha}) {
    if (a->size() != n) throw std::runtime_error("every per-gene input must have the same length");
  }
  if (has_responder && u_raw.size() != n) {
    throw std::runtime_error("a condition fit needs u_raw, one value per gene");
  }
  ColumnMajor score(n), score4(n), v_resid(n), v_s(n), v_total(n), r2_s(n);
  ColumnMajor v_rxs(has_responder ? n : 0), r2_rxs(has_responder ? n : 0);
  raise_if_failed(pace::driver_scores(span_of(estimate_shrunk), span_of(spec), span_of(focal_mean),
                                      span_of(mu_bar), span_of(alpha), span_of(u_raw), var_n,
                                      var_rn, has_responder, n, out_span_of(score),
                                      out_span_of(score4), out_span_of(v_resid), out_span_of(v_s),
                                      out_span_of(v_rxs), out_span_of(v_total), out_span_of(r2_s),
                                      out_span_of(r2_rxs)),
                  "driver scores");
  py::dict out;
  out["driver_score"] = score;
  out["driver_score_spec4"] = score4;
  out["V_resid"] = v_resid;
  out["V_S"] = v_s;
  out["V_RxS"] = v_rxs;
  out["V_total"] = v_total;
  out["R2_S"] = r2_s;
  out["R2_RxS"] = r2_rxs;
  return out;
}

// Pratt attribution of a focal type's spatial variance across its neighbours.
// `sigma` is the n_pairs x n_pairs covariance of the neighbour kernels and `u`
// the n_pairs x n_genes slopes. This is what the pair heatmap reads.
py::dict pair_variance_pratt(const ColumnMajor& sigma, const ColumnMajor& u) {
  if (sigma.ndim() != 2 || sigma.shape(0) != sigma.shape(1)) {
    throw std::runtime_error("sigma must be a square n_pairs x n_pairs matrix");
  }
  if (u.ndim() != 2) throw std::runtime_error("u must be n_pairs x n_genes");
  const std::int64_t n_pairs = sigma.shape(0);
  if (u.shape(0) != n_pairs) throw std::runtime_error("u must have one row per pair");
  const std::int64_t n_genes = u.shape(1);

  ColumnMajor v_pair(n_pairs), v_diag(n_pairs), share_pct(n_pairs), diag_pct(n_pairs), sign(n_pairs);
  double v_total = 0.0, v_diag_total = 0.0, cross_cov_pct = 0.0;
  int n_negative = 0;
  raise_if_failed(pace::pair_variance_pratt(span_of(sigma), span_of(u), n_pairs, n_genes,
                                            out_span_of(v_pair), out_span_of(v_diag),
                                            out_span_of(share_pct), out_span_of(diag_pct),
                                            out_span_of(sign), &v_total, &v_diag_total,
                                            &cross_cov_pct, &n_negative),
                  "pair variance");
  py::dict out;
  out["v_pair"] = v_pair;
  out["v_diag"] = v_diag;
  out["share_pct"] = share_pct;
  out["diag_pct"] = diag_pct;
  out["sign"] = sign;
  out["v_total"] = v_total;
  out["v_diag_total"] = v_diag_total;
  out["cross_cov_pct"] = cross_cov_pct;
  out["n_negative"] = n_negative;
  return out;
}

using IntColumn = py::array_t<int, py::array::c_style | py::array::forcecast>;

// Gaussian biological and exponential technical neighbour kernels, the first
// stage of a fit. For each cell, the kernel-weighted count of every other cell
// type within the radius.
//
// `cell_type` is 0-based type codes, one per cell. `group` assigns cells to
// images; pass `per_group = true` to build the kernels within each image
// separately, which is what a cohort of separate samples needs. `eps` is the
// radius, defaulting to 3 * h_bio as the method specifies.
py::dict neighbour_kernels(const ColumnMajor& x, const ColumnMajor& y,
                           const IntColumn& cell_type, int n_types, const IntColumn& group,
                           bool per_group, double h_bio, double h_tech, double eps,
                           int n_threads) {
  const std::int64_t n_cells = x.size();
  if (y.size() != n_cells || cell_type.size() != n_cells) {
    throw std::runtime_error("x, y and cell_type must have one entry per cell");
  }
  if (group.size() != n_cells) throw std::runtime_error("group must have one entry per cell");
  if (eps <= 0.0) eps = 3.0 * h_bio;

  ColumnMajor k_bio({n_cells, static_cast<std::int64_t>(n_types)});
  ColumnMajor k_tech({n_cells, static_cast<std::int64_t>(n_types)});
  const pace::Status status = pace::neighbour_kernels(
      span_of(x), span_of(y),
      pace::Span<const int>(cell_type.data(), n_cells), n_types,
      pace::Span<const int>(group.data(), n_cells), per_group,
      h_bio, h_tech, eps, out_span_of(k_bio), out_span_of(k_tech),
      n_threads, pace::InterruptCheck());
  raise_if_failed(status, "neighbour kernels");

  py::dict out;
  out["k_bio"] = k_bio;
  out["k_tech"] = k_tech;
  out["eps"] = eps;
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
  m.def("driver_scores", &driver_scores, py::arg("estimate_shrunk"), py::arg("spec"),
        py::arg("focal_mean"), py::arg("mu_bar"), py::arg("alpha"),
        py::arg("u_raw") = ColumnMajor(0), py::arg("var_n") = 0.0, py::arg("var_rn") = 0.0,
        py::arg("has_responder") = false,
        "Driver scores for one focal-neighbour pair.");
  m.def("pair_variance_pratt", &pair_variance_pratt, py::arg("sigma"), py::arg("u"),
        "Pratt attribution of spatial variance across a focal type's neighbours.");
  m.def("neighbour_kernels", &neighbour_kernels, py::arg("x"), py::arg("y"),
        py::arg("cell_type"), py::arg("n_types"), py::arg("group"),
        py::arg("per_group") = false, py::arg("h_bio") = 30.0, py::arg("h_tech") = 5.0,
        py::arg("eps") = 0.0, py::arg("n_threads") = 1,
        "Gaussian biological and exponential technical neighbour kernels.");
}
