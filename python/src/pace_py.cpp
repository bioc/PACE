// pybind11 bindings for the PACE core.
//
// The core is R-free by construction, so the same translation units back both
// the R package and this module. Arrays cross the boundary as numpy arrays in
// Fortran order, which is what the core expects and what R already uses, so
// neither side copies.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "decomposition.hpp"
#include "neighbourhood.hpp"
#include "random_design.hpp"
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

// Rcpp's NumericMatrix arrives zero-filled and parts of the core rely on that,
// writing only the entries they compute. pybind11 does NOT zero a new array, so
// zero it here -- otherwise untouched entries are whatever was on the heap.
pace::Span<double> out_span_of(ColumnMajor& array) {
  double* data = array.mutable_data();
  std::fill(data, data + array.size(), 0.0);
  return pace::Span<double>(data, static_cast<std::int64_t>(array.size()));
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

// Sparse cross-type ambient weight matrix W, the contamination field:
//
//   W[i, j] = exp(-d_ij / h_tech) / edge_fraction[i]
//
// over cells i, j of the same image, different type, within 3 * h_tech.
// Returned in compressed-column form, so
//
//   scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)
//
// rebuilds it without this module depending on scipy.
//
// The edge correction integrates over a quadrature of `n_angles` directions,
// which R fixes at 1000; it is built here rather than asked of the caller.
py::dict ambient_field(const ColumnMajor& x, const ColumnMajor& y, const IntColumn& cell_type,
                       const IntColumn& image, int n_images, double h_tech, bool edge_correct,
                       int n_angles, int n_threads) {
  const std::int64_t n = x.size();
  if (y.size() != n || cell_type.size() != n || image.size() != n) {
    throw std::runtime_error("x, y, cell_type and image must have one entry per cell");
  }
  if (n_angles < 1) throw std::runtime_error("n_angles must be positive");

  // theta_k = 2 pi k / n_angles, matching the R quadrature exactly
  std::vector<double> cos_theta(static_cast<std::size_t>(n_angles));
  std::vector<double> sin_theta(static_cast<std::size_t>(n_angles));
  const double two_pi = 6.283185307179586476925286766559;
  for (int k = 0; k < n_angles; ++k) {
    const double theta = two_pi * static_cast<double>(k + 1) / static_cast<double>(n_angles);
    cos_theta[static_cast<std::size_t>(k)] = std::cos(theta);
    sin_theta[static_cast<std::size_t>(k)] = std::sin(theta);
  }

  ColumnMajor edge_fraction(n);
  pace::AmbientFieldBuilder builder;
  raise_if_failed(builder.prepare(span_of(x), span_of(y),
                                  pace::Span<const int>(cell_type.data(), n),
                                  pace::Span<const int>(image.data(), n), n_images, h_tech,
                                  edge_correct,
                                  pace::Span<const double>(cos_theta.data(), n_angles),
                                  pace::Span<const double>(sin_theta.data(), n_angles),
                                  out_span_of(edge_fraction), n_threads, pace::InterruptCheck()),
                  "ambient field");

  const std::int64_t nnz = builder.nnz();
  py::array_t<int> indices(nnz);
  std::fill(indices.mutable_data(), indices.mutable_data() + nnz, 0);
  ColumnMajor data(nnz);
  py::array_t<std::int64_t> indptr(n + 1);
  const std::vector<std::int64_t>& pointer = builder.column_pointer();
  auto indptr_out = indptr.mutable_unchecked<1>();
  for (std::int64_t j = 0; j <= n; ++j) indptr_out(j) = pointer[static_cast<std::size_t>(j)];

  raise_if_failed(builder.fill(pace::Span<int>(indices.mutable_data(), nnz), out_span_of(data),
                               pace::Span<const double>(edge_fraction.data(), n), n_threads,
                               pace::InterruptCheck()),
                  "ambient field");

  py::dict out;
  out["data"] = data;
  out["indices"] = indices;
  out["indptr"] = indptr;
  out["shape"] = py::make_tuple(n, n);
  out["edge_fraction"] = edge_fraction;
  out["nnz"] = nnz;
  return out;
}

// Zero a freshly allocated integer array, for the same reason out_span_of()
// zeroes a double one.
pace::Span<int> out_int_span(py::array_t<int>& array) {
  int* data = array.mutable_data();
  std::fill(data, data + array.size(), 0);
  return pace::Span<int>(data, static_cast<std::int64_t>(array.size()));
}

// One block of the random-effect design: the cells-by-terms matrix spread
// across groups, so column (group, term) carries that group's cells only.
// Returned in compressed-column form, with the grouping the core computed on
// the same pass -- `cells_by_group` holds cell indices ordered by group and
// `group_start` indexes into it, so group g owns
// cells_by_group[group_start[g]:group_start[g + 1]].
py::dict random_design_block(const ColumnMajor& x_terms, const IntColumn& cell_group,
                             int n_groups) {
  if (x_terms.ndim() != 2) throw std::runtime_error("x_terms must be n_cells x n_terms");
  const std::int64_t n = x_terms.shape(0);
  const std::int64_t n_terms = x_terms.shape(1);
  if (cell_group.size() != n) throw std::runtime_error("cell_group must have one entry per cell");
  if (n_groups < 1) throw std::runtime_error("n_groups must be positive");
  const std::int64_t n_columns = n_terms * static_cast<std::int64_t>(n_groups);
  const std::int64_t nnz = n * n_terms;

  py::array_t<int> column_pointer(n_columns + 1);
  py::array_t<int> row_index(nnz);
  ColumnMajor values(nnz);
  py::array_t<int> cells_by_group(n);
  py::array_t<int> group_start(static_cast<std::int64_t>(n_groups) + 1);

  raise_if_failed(pace::random_design_block(
                      span_of(x_terms), pace::Span<const int>(cell_group.data(), n), n, n_terms,
                      n_groups, out_int_span(column_pointer), out_int_span(row_index),
                      out_span_of(values), out_int_span(cells_by_group),
                      out_int_span(group_start)),
                  "random design block");

  py::dict out;
  out["data"] = values;
  out["indices"] = row_index;
  out["indptr"] = column_pointer;
  out["shape"] = py::make_tuple(n, n_columns);
  out["cells_by_group"] = cells_by_group;
  out["group_start"] = group_start;
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
  m.def("ambient_field", &ambient_field, py::arg("x"), py::arg("y"), py::arg("cell_type"),
        py::arg("image"), py::arg("n_images"), py::arg("h_tech") = 5.0,
        py::arg("edge_correct") = true, py::arg("n_angles") = 1000, py::arg("n_threads") = 1,
        "Sparse cross-type ambient weight matrix, in compressed-column form.");
  m.def("random_design_block", &random_design_block, py::arg("x_terms"),
        py::arg("cell_group"), py::arg("n_groups"),
        "One block of the random-effect design, in compressed-column form.");
}
