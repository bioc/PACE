// bindings.cpp -- thin Rcpp layer over the R-free core in src/core/.
//
// This file only (1) converts R objects to core spans without copying, (2)
// allocates R outputs on the main thread, (3) supplies the interrupt check, and
// (4) turns a core Status into an R error or interrupt. No computation lives here.
#include <Rcpp.h>
#include <Rmath.h>

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "core/mash.hpp"
#include "core/count_stats.hpp"
#include "core/core_types.hpp"
#include "core/decomposition.hpp"
#include "core/dispersion.hpp"
#include "core/gene_solve.hpp"
#include "core/hyperparameters.hpp"
#include "core/irls_chunk.hpp"
#include "core/irls_driver.hpp"
#include "core/neighbourhood.hpp"
#include "core/linear_predictor.hpp"
#include "core/preprocess.hpp"
#include "core/sparse_product.hpp"
#include "core/statistics.hpp"

namespace {

// R_CheckUserInterrupt() longjmps; running it under R_ToplevelExec() turns a
// pending interrupt into a FALSE return instead.
void check_interrupt_callback(void*) { R_CheckUserInterrupt(); }

bool user_interrupted() { return R_ToplevelExec(check_interrupt_callback, nullptr) == FALSE; }

// Raise the core status in R: an interrupt stays an interrupt, anything else is an error.
void raise_if_failed(const pace::Status& status, const char* where) {
  if (status.is_ok()) return;
  if (status.code == pace::StatusCode::interrupted) throw Rcpp::internal::InterruptedException();
  Rcpp::stop(std::string(where) + ": " + status.message);
}

pace::Span<const double> column_span(const Rcpp::NumericMatrix& matrix, int column) {
  const std::int64_t n_rows = matrix.nrow();
  return pace::Span<const double>(REAL(matrix) + column * n_rows, n_rows);
}

pace::Span<const int> int_span(const Rcpp::IntegerVector& vector) {
  return pace::Span<const int>(INTEGER(vector), Rf_xlength(vector));
}

pace::Span<const double> double_span(const Rcpp::NumericVector& vector) {
  return pace::Span<const double>(REAL(vector), Rf_xlength(vector));
}

void check_coordinates(const Rcpp::NumericMatrix& coords) {
  if (coords.ncol() != 2) Rcpp::stop("coords must have exactly two columns");
}

// The dgCMatrix slots, held so the core's view stays valid for the call.
struct CscHolder {
  Rcpp::IntegerVector dims;
  Rcpp::IntegerVector column_pointer;
  Rcpp::IntegerVector row_index;
  Rcpp::NumericVector values;
  pace::CscView view;

  explicit CscHolder(const Rcpp::S4& matrix)
      : dims(matrix.slot("Dim")),
        column_pointer(matrix.slot("p")),
        row_index(matrix.slot("i")),
        values(matrix.slot("x")) {
    view.column_pointer = int_span(column_pointer);
    view.row_index = int_span(row_index);
    view.values = double_span(values);
    view.n_rows = dims[0];
    view.n_cols = dims[1];
  }
};

pace::Span<double> out_span(Rcpp::NumericVector& vector) {
  return pace::Span<double>(REAL(vector), Rf_xlength(vector));
}

pace::Span<int> out_span(Rcpp::IntegerVector& vector) {
  return pace::Span<int>(INTEGER(vector), Rf_xlength(vector));
}

pace::Span<double> out_span(Rcpp::NumericMatrix& matrix) {
  return pace::Span<double>(REAL(matrix), Rf_xlength(matrix));
}

pace::Span<const double> const_span(const Rcpp::NumericMatrix& matrix) {
  return pace::Span<const double>(REAL(matrix), Rf_xlength(matrix));
}

pace::Span<const int> const_span(const Rcpp::IntegerMatrix& matrix) {
  return pace::Span<const int>(INTEGER(matrix), Rf_xlength(matrix));
}

// R's NA for the slots where stats::var() has fewer than two values.
void set_na_where_too_few(Rcpp::NumericVector& variance, const Rcpp::IntegerVector& n_used) {
  for (R_xlen_t k = 0; k < variance.size(); ++k) {
    if (n_used[k] < 2) variance[k] = NA_REAL;
  }
}

}  // namespace

// Neighbour kernels (see pace::neighbour_kernels). coords: n x 2 double;
// neighbour_type, group: 0-based integer codes, -1 for none. Returns a list with
// K_bio and K_tech, n x n_types double matrices.
// [[Rcpp::export]]
Rcpp::List pace_neighbour_kernels_cpp(const Rcpp::NumericMatrix& coords,
                                      const Rcpp::IntegerVector& neighbour_type, int n_types,
                                      const Rcpp::IntegerVector& group, bool per_group,
                                      double h_bio, double h_tech, double eps, int n_threads) {
  check_coordinates(coords);
  const std::int64_t n = coords.nrow();
  Rcpp::NumericMatrix k_bio(n, n_types);
  Rcpp::NumericMatrix k_tech(n, n_types);
  const pace::Status status = pace::neighbour_kernels(
      column_span(coords, 0), column_span(coords, 1), int_span(neighbour_type), n_types,
      int_span(group), per_group, h_bio, h_tech, eps,
      pace::Span<double>(REAL(k_bio), n * n_types), pace::Span<double>(REAL(k_tech), n * n_types),
      n_threads, user_interrupted);
  raise_if_failed(status, "neighbour kernels");
  return Rcpp::List::create(Rcpp::Named("K_bio") = k_bio, Rcpp::Named("K_tech") = k_tech);
}

// Neighbour count per cell under the kernel search rules (test helper).
// [[Rcpp::export]]
Rcpp::IntegerVector pace_neighbour_counts_cpp(const Rcpp::NumericMatrix& coords,
                                              const Rcpp::IntegerVector& group, bool per_group,
                                              double eps, int n_threads) {
  check_coordinates(coords);
  Rcpp::IntegerVector counts(coords.nrow());
  const pace::Status status = pace::neighbour_counts(
      column_span(coords, 0), column_span(coords, 1), int_span(group), per_group, eps,
      pace::Span<int>(INTEGER(counts), Rf_xlength(counts)), n_threads, user_interrupted);
  raise_if_failed(status, "neighbour counts");
  return counts;
}

// Neighbour lists (test helper): list(offsets = n + 1 offsets stored as doubles
// (exact up to 2^53, beyond the int range), neighbours = 1-based cell indices sorted
// within each cell, distances).
// [[Rcpp::export]]
Rcpp::List pace_neighbour_lists_cpp(const Rcpp::NumericMatrix& coords, const Rcpp::IntegerVector& group,
                                    bool per_group, double eps, int n_threads) {
  check_coordinates(coords);
  std::vector<std::int64_t> offsets;
  std::vector<int> neighbours;
  std::vector<double> distances;
  const pace::Status status = pace::neighbour_lists(
      column_span(coords, 0), column_span(coords, 1), int_span(group), per_group, eps,
      offsets, neighbours, distances, n_threads, user_interrupted);
  raise_if_failed(status, "neighbour lists");
  Rcpp::NumericVector offsets_out(offsets.begin(), offsets.end());
  Rcpp::IntegerVector neighbours_out(neighbours.size());
  for (std::size_t k = 0; k < neighbours.size(); ++k) neighbours_out[k] = neighbours[k] + 1;
  Rcpp::NumericVector distances_out(distances.begin(), distances.end());
  return Rcpp::List::create(Rcpp::Named("offsets") = offsets_out,
                            Rcpp::Named("neighbours") = neighbours_out,
                            Rcpp::Named("distances") = distances_out);
}

// Edge fractions for all rows of coords against one rectangle (see pace::area_fractions).
// [[Rcpp::export]]
Rcpp::NumericVector pace_area_fraction_cpp(const Rcpp::NumericMatrix& coords, double r,
                                           double x_min, double x_max, double y_min, double y_max,
                                           const Rcpp::NumericVector& cos_theta,
                                           const Rcpp::NumericVector& sin_theta, int n_threads) {
  check_coordinates(coords);
  const std::int64_t n = coords.nrow();
  Rcpp::NumericVector fraction(n);
  if (n == 0) return fraction;
  Rcpp::IntegerVector rows = Rcpp::seq(0, static_cast<int>(n) - 1);
  const pace::Status status = pace::area_fractions(
      column_span(coords, 0), column_span(coords, 1), int_span(rows), r, x_min, x_max, y_min, y_max,
      double_span(cos_theta), double_span(sin_theta), pace::Span<double>(REAL(fraction), n),
      n_threads, user_interrupted);
  raise_if_failed(status, "area fraction");
  return fraction;
}

// Ambient weight matrix W (see pace::AmbientFieldBuilder). type_code, image: 0-based
// codes (image -1: no image). Returns the dgCMatrix slots i, p, x and the edge fractions.
// [[Rcpp::export]]
Rcpp::List pace_ambient_field_cpp(const Rcpp::NumericMatrix& coords,
                                  const Rcpp::IntegerVector& type_code,
                                  const Rcpp::IntegerVector& image, int n_images, double h_tech,
                                  bool edge_correct, const Rcpp::NumericVector& cos_theta,
                                  const Rcpp::NumericVector& sin_theta, int n_threads) {
  check_coordinates(coords);
  const std::int64_t n = coords.nrow();
  Rcpp::NumericVector edge_fraction(n);
  pace::AmbientFieldBuilder builder;
  pace::Status status = builder.prepare(
      column_span(coords, 0), column_span(coords, 1), int_span(type_code), int_span(image), n_images,
      h_tech, edge_correct, double_span(cos_theta), double_span(sin_theta),
      pace::Span<double>(REAL(edge_fraction), n), n_threads, user_interrupted);
  raise_if_failed(status, "ambient field");

  const std::int64_t nnz = builder.nnz();
  Rcpp::IntegerVector row_index(nnz);
  Rcpp::NumericVector values(nnz);
  Rcpp::IntegerVector column_pointer(n + 1);
  const std::vector<std::int64_t>& pointer = builder.column_pointer();
  for (std::int64_t j = 0; j <= n; ++j) column_pointer[j] = static_cast<int>(pointer[j]);
  status = builder.fill(pace::Span<int>(INTEGER(row_index), nnz), pace::Span<double>(REAL(values), nnz),
                        pace::Span<const double>(REAL(edge_fraction), n), n_threads, user_interrupted);
  raise_if_failed(status, "ambient field");
  return Rcpp::List::create(Rcpp::Named("i") = row_index, Rcpp::Named("p") = column_pointer,
                            Rcpp::Named("x") = values, Rcpp::Named("edge_fraction") = edge_fraction);
}

// Same-type neighbour fraction per cell (see pace::same_type_fraction); NA where
// the cell's image is too small or missing.
// [[Rcpp::export]]
Rcpp::NumericVector pace_same_type_fraction_cpp(const Rcpp::NumericMatrix& coords,
                                                const Rcpp::IntegerVector& type_code,
                                                const Rcpp::IntegerVector& image,
                                                double radius, int min_image_cells, int n_threads) {
  check_coordinates(coords);
  const std::int64_t n = coords.nrow();
  Rcpp::NumericVector fraction(n);
  const pace::Status status = pace::same_type_fraction(
      column_span(coords, 0), column_span(coords, 1), int_span(type_code), int_span(image),
      radius, min_image_cells, pace::Span<double>(REAL(fraction), n), n_threads, user_interrupted);
  raise_if_failed(status, "same-type fraction");
  for (std::int64_t i = 0; i < n; ++i) {
    if (ISNAN(fraction[i])) fraction[i] = NA_REAL;
  }
  return fraction;
}

// Per-group column means of a cells x genes dgCMatrix (see pace::group_column_means).
// Returns an n_groups x n_genes double matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix pace_group_column_means_cpp(const Rcpp::S4& counts, const Rcpp::IntegerVector& group,
                                                int n_groups, bool detection, int n_threads) {
  const Rcpp::IntegerVector dims = counts.slot("Dim");
  const Rcpp::IntegerVector column_pointer = counts.slot("p");
  const Rcpp::IntegerVector row_index = counts.slot("i");
  const Rcpp::NumericVector values = counts.slot("x");
  pace::CscView view;
  view.column_pointer = int_span(column_pointer);
  view.row_index = int_span(row_index);
  view.values = double_span(values);
  view.n_rows = dims[0];
  view.n_cols = dims[1];
  Rcpp::NumericMatrix means(n_groups, dims[1]);
  const pace::Status status = pace::group_column_means(
      view, int_span(group), n_groups, detection,
      pace::Span<double>(REAL(means), static_cast<std::int64_t>(n_groups) * dims[1]), n_threads,
      user_interrupted);
  raise_if_failed(status, "group column means");
  return means;
}

// Per-group covariance of the columns of a dense n x p matrix (see
// pace::group_covariances). Returns a p x p x n_groups array (slice [, , g] is
// group g's covariance matrix), matching the core's contiguous per-group layout.
// [[Rcpp::export]]
Rcpp::NumericVector pace_group_covariances_cpp(const Rcpp::NumericMatrix& values,
                                               const Rcpp::IntegerVector& group, int n_groups,
                                               int n_threads) {
  const std::int64_t n = values.nrow();
  const std::int64_t p = values.ncol();
  Rcpp::NumericVector covariance(static_cast<R_xlen_t>(n_groups) * p * p);
  const pace::Status status = pace::group_covariances(
      pace::Span<const double>(REAL(values), n * p), n, p, int_span(group), n_groups,
      pace::Span<double>(REAL(covariance), Rf_xlength(covariance)), n_threads, user_interrupted);
  raise_if_failed(status, "group covariances");
  covariance.attr("dim") = Rcpp::IntegerVector::create(static_cast<int>(p), static_cast<int>(p), n_groups);
  return covariance;
}

// Single-frame statistics of a cells x genes dgCMatrix (see
// pace::single_frame_statistics). Returns focal_mean and within_ss
// (n_groups x n_genes) and global_mean (n_genes).
// [[Rcpp::export]]
Rcpp::List pace_single_frame_statistics_cpp(const Rcpp::S4& counts, const Rcpp::NumericVector& n_count,
                                            const Rcpp::IntegerVector& group, int n_groups,
                                            int n_threads) {
  const Rcpp::IntegerVector dims = counts.slot("Dim");
  const Rcpp::IntegerVector column_pointer = counts.slot("p");
  const Rcpp::IntegerVector row_index = counts.slot("i");
  const Rcpp::NumericVector values = counts.slot("x");
  pace::CscView view;
  view.column_pointer = int_span(column_pointer);
  view.row_index = int_span(row_index);
  view.values = double_span(values);
  view.n_rows = dims[0];
  view.n_cols = dims[1];
  Rcpp::NumericMatrix focal_mean(n_groups, dims[1]);
  Rcpp::NumericMatrix within_ss(n_groups, dims[1]);
  Rcpp::NumericVector global_mean(dims[1]);
  const std::int64_t n_slots = static_cast<std::int64_t>(n_groups) * dims[1];
  const pace::Status status = pace::single_frame_statistics(
      view, double_span(n_count), int_span(group), n_groups,
      pace::Span<double>(REAL(focal_mean), n_slots), pace::Span<double>(REAL(within_ss), n_slots),
      pace::Span<double>(REAL(global_mean), dims[1]), n_threads, user_interrupted);
  raise_if_failed(status, "single-frame statistics");
  return Rcpp::List::create(Rcpp::Named("focal_mean") = focal_mean,
                            Rcpp::Named("within_ss") = within_ss,
                            Rcpp::Named("global_mean") = global_mean);
}

namespace {

// A gene range of a cells x genes CSC matrix; `first_gene` is 1-based.
pace::GeneBlock gene_block(const CscHolder& holder, int first_gene) {
  pace::GeneBlock block;
  block.matrix = holder.view;
  block.first_gene = first_gene - 1;
  return block;
}

// The ambient field as the core now takes it. These wrappers are the CACHED
// path -- R has already built W %*% Y and hands the whole thing over -- so the
// source just carries that block. The streamed path is configured in
// pace_irls_driver_cpp(), which is the only caller that has W and Y separately.
pace::AmbientSource cached_ambient(const CscHolder& holder, int first_gene) {
  pace::AmbientSource source;
  source.streamed = false;
  source.cached = gene_block(holder, first_gene);
  return source;
}

}  // namespace

// Per-group mean and variance of the columns of a dense matrix (see
// pace::dense_group_moments). `values` may be a 0 x 0 matrix, which means an
// all-zero n x p matrix. Returns mean, variance (R's NA where fewer than two
// values entered it) and any_nonzero.
// [[Rcpp::export]]
Rcpp::List pace_dense_group_moments_cpp(const Rcpp::NumericMatrix& values, int n, int p,
                                        const Rcpp::IntegerVector& group, int n_groups,
                                        bool want_mean, bool want_variance, int n_threads) {
  const bool all_zero = values.size() == 0;
  if (!all_zero && (values.nrow() != n || values.ncol() != p)) {
    Rcpp::stop("values must be empty or n x p");
  }
  Rcpp::NumericMatrix mean(want_mean ? n_groups : 0, want_mean ? p : 0);
  Rcpp::NumericMatrix variance(want_variance ? n_groups : 0, want_variance ? p : 0);
  Rcpp::IntegerVector n_used(want_variance ? static_cast<R_xlen_t>(n_groups) * p : 0);
  bool any_nonzero = false;
  const pace::Status status = pace::dense_group_moments(
      all_zero ? pace::Span<const double>() : const_span(values), n, p, int_span(group), n_groups,
      out_span(mean), out_span(variance), out_span(n_used), &any_nonzero, n_threads, user_interrupted);
  raise_if_failed(status, "group moments");
  if (want_variance) set_na_where_too_few(variance, n_used);
  return Rcpp::List::create(Rcpp::Named("mean") = mean, Rcpp::Named("variance") = variance,
                            Rcpp::Named("any_nonzero") = any_nonzero);
}

// One gene chunk of the streaming solver's final pass (see pace::final_pass_statistics).
// `return_matrices` also returns this chunk's mu and technical offset.
// [[Rcpp::export]]
Rcpp::List pace_final_pass_statistics_cpp(const Rcpp::NumericMatrix& eta,
                                          const Rcpp::NumericVector& offset,
                                          const Rcpp::S4& ambient, int first_gene,
                                          const Rcpp::NumericVector& rho,
                                          const Rcpp::IntegerVector& group, int n_groups,
                                          bool want_groups, bool return_matrices) {
  CscHolder ambient_holder(ambient);
  const std::int64_t n = eta.nrow();
  const std::int64_t n_genes = eta.ncol();
  Rcpp::NumericMatrix mu_group_sum(want_groups ? n_groups : 0, want_groups ? n_genes : 0);
  Rcpp::NumericMatrix toff_variance(want_groups ? n_groups : 0, want_groups ? n_genes : 0);
  Rcpp::NumericVector mu_column_sum(n_genes);
  Rcpp::NumericVector spill_row_sum(n);
  Rcpp::NumericVector total_row_sum(n);
  Rcpp::NumericMatrix mu(return_matrices ? n : 0, return_matrices ? n_genes : 0);
  Rcpp::NumericMatrix toff(return_matrices ? n : 0, return_matrices ? n_genes : 0);
  bool any_nonzero = false;
  const pace::Status status = pace::final_pass_statistics(
      const_span(eta), n, n_genes, double_span(offset), gene_block(ambient_holder, first_gene),
      double_span(rho), int_span(group), n_groups, out_span(mu_group_sum), out_span(toff_variance),
      out_span(mu_column_sum), out_span(spill_row_sum), out_span(total_row_sum), out_span(mu),
      out_span(toff), &any_nonzero, user_interrupted);
  raise_if_failed(status, "final pass statistics");
  // stats::var() of fewer than two values is NA, not NaN.
  std::vector<int> group_size(want_groups ? n_groups : 0, 0);
  if (want_groups) {
    for (R_xlen_t i = 0; i < group.size(); ++i) {
      const int g = group[i];
      if (g >= 0 && g < n_groups) group_size[g] += 1;
    }
    for (std::int64_t j = 0; j < n_genes; ++j) {
      for (int g = 0; g < n_groups; ++g) {
        if (group_size[g] < 2) toff_variance[g + j * n_groups] = NA_REAL;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("mu_group_sum") = mu_group_sum,
                            Rcpp::Named("toff_variance") = toff_variance,
                            Rcpp::Named("mu_column_sum") = mu_column_sum,
                            Rcpp::Named("spill_row_sum") = spill_row_sum,
                            Rcpp::Named("total_row_sum") = total_row_sum,
                            Rcpp::Named("any_nonzero") = any_nonzero,
                            Rcpp::Named("mu") = mu, Rcpp::Named("toff") = toff);
}

// The ambient product A = W Y over a contiguous range of the right operand's
// columns (see pace::sparse_product_csc). Operands and result are dgCMatrix;
// `first_column` is 1-based, and the result's columns are numbered from one.
//
// R's `%*%` computes the same product, but it is a SECOND implementation, in a
// translation unit this package does not control and cannot set flags on. Cache
// mode used to take its product from there while streamed mode used this code,
// so the two modes agreed only as long as two independent implementations
// happened to round identically. They did not: the modes differed by about
// 1.5e-8 in the coefficients on every aarch64 build platform while agreeing on
// both x86-64 ones. The likeliest reason is fused multiply-add -- aarch64 has
// one in its baseline ISA and x86-64 does not, and this package's -ffp-contract
// setting cannot reach Matrix -- but that was never confirmed on a failing
// runner, and a different accumulation inside Matrix would look the same.
//
// Which it is does not matter here. Routing BOTH modes through this function
// leaves them differing only in chunking, and a column range is exactly the
// same columns of the whole product because a sparse product is
// column-independent. That is an invariant of this file rather than a
// coincidence between two libraries.
// [[Rcpp::export]]
Rcpp::S4 pace_sparse_product_cpp(const Rcpp::S4& left, const Rcpp::S4& right, int first_column,
                                 int n_columns, int n_threads) {
  CscHolder left_holder(left);
  CscHolder right_holder(right);
  std::vector<int> column_pointer, row_index;
  std::vector<double> values;
  const pace::Status status = pace::sparse_product_csc(
      left_holder.view, right_holder.view, static_cast<std::int64_t>(first_column) - 1, n_columns,
      column_pointer, row_index, values, n_threads, user_interrupted);
  raise_if_failed(status, "sparse product");
  Rcpp::S4 out("dgCMatrix");
  out.slot("Dim") =
      Rcpp::IntegerVector::create(static_cast<int>(left_holder.view.n_rows), n_columns);
  out.slot("p") = Rcpp::IntegerVector(column_pointer.begin(), column_pointer.end());
  out.slot("i") = Rcpp::IntegerVector(row_index.begin(), row_index.end());
  out.slot("x") = Rcpp::NumericVector(values.begin(), values.end());
  // Carry dimnames the way `%*%` does -- the left operand's rows, and the
  // SELECTED columns of the right operand. Dropping them made this an almost
  // drop-in replacement, which is worse than either alternative: the readout
  // rebuild returns its matrix to the caller, and its gene names went missing
  // while every value stayed identical.
  Rcpp::RObject left_names = Rcpp::S4(left).slot("Dimnames");
  Rcpp::RObject right_names = Rcpp::S4(right).slot("Dimnames");
  Rcpp::RObject out_rows = R_NilValue, out_cols = R_NilValue;
  if (left_names.inherits("list") || TYPEOF(left_names) == VECSXP) {
    const Rcpp::List names(left_names);
    if (names.size() > 0) out_rows = names[0];
  }
  if (right_names.inherits("list") || TYPEOF(right_names) == VECSXP) {
    const Rcpp::List names(right_names);
    if (names.size() > 1 && !Rf_isNull(names[1])) {
      const Rcpp::CharacterVector all(names[1]);
      Rcpp::CharacterVector picked(n_columns);
      for (int j = 0; j < n_columns; ++j) picked[j] = all[first_column - 1 + j];
      out_cols = picked;
    }
  }
  out.slot("Dimnames") = Rcpp::List::create(out_rows, out_cols);
  return out;
}

// Per-(focal, gene) variance decomposition (see pace::variance_decomposition).
// Index arguments are 1-based; a 0 x 0 matrix means the block is absent. Every
// returned vector has n_focals * n_genes entries, focal-major.
// [[Rcpp::export]]
Rcpp::List pace_variance_decomposition_cpp(
    const Rcpp::NumericMatrix& ct_means, const Rcpp::IntegerVector& group_size,
    const Rcpp::IntegerVector& focal_group, const Rcpp::IntegerVector& n_focal,
    const Rcpp::NumericMatrix& u, const Rcpp::NumericMatrix& se_u,
    const Rcpp::IntegerMatrix& slope_rows, const Rcpp::NumericVector& kernel_cov,
    const Rcpp::IntegerMatrix& responder_rows, const Rcpp::IntegerVector& responder_keep,
    const Rcpp::NumericVector& responder_cov, const Rcpp::IntegerVector& intercept_rows,
    const Rcpp::NumericMatrix& toff_var, const Rcpp::NumericVector& spill_cov,
    const Rcpp::NumericMatrix& beta_spill, const Rcpp::NumericMatrix& mu_mean,
    const Rcpp::NumericVector& alpha, bool nb1, int n_threads) {
  pace::DecompositionInput input;
  input.n_genes = ct_means.ncol();
  input.n_groups = group_size.size();
  input.n_focals = focal_group.size();
  input.ct_means = const_span(ct_means);
  input.group_size = int_span(group_size);
  input.focal_group = int_span(focal_group);
  input.n_focal = int_span(n_focal);
  input.u = const_span(u);
  input.se_u = const_span(se_u);
  input.q = u.nrow();
  input.n_kernel = slope_rows.nrow();
  input.slope_rows = const_span(slope_rows);
  input.kernel_cov = double_span(kernel_cov);
  input.n_responder = responder_rows.nrow();
  input.responder_rows = const_span(responder_rows);
  input.responder_keep = int_span(responder_keep);
  input.responder_cov = double_span(responder_cov);
  input.intercept_rows = int_span(intercept_rows);
  input.toff_var = const_span(toff_var);
  input.n_spill = beta_spill.nrow();
  input.spill_cov = double_span(spill_cov);
  input.beta_spill = const_span(beta_spill);
  input.mu_mean = const_span(mu_mean);
  input.alpha = double_span(alpha);
  input.nb1 = nb1;

  const R_xlen_t rows = static_cast<R_xlen_t>(input.n_focals) * input.n_genes;
  Rcpp::NumericVector v_state_baseline(rows), v_state_responder(rows), v_spill(rows), v_disp(rows);
  Rcpp::NumericVector celltype_offset_sq(rows), focal_mean(rows), max_other_mean(rows), spec(rows);
  Rcpp::NumericVector focal_other_ratio(rows), pct_celltype(rows), pct_state(rows);
  Rcpp::NumericVector pct_responder(rows), pct_spill(rows), pct_residual(rows);
  Rcpp::IntegerVector is_contaminated(rows), keep(rows);
  pace::DecompositionOutput output;
  output.v_state_baseline = out_span(v_state_baseline);
  output.v_state_responder = out_span(v_state_responder);
  output.v_spill = out_span(v_spill);
  output.v_disp = out_span(v_disp);
  output.celltype_offset_sq = out_span(celltype_offset_sq);
  output.focal_mean = out_span(focal_mean);
  output.max_other_mean = out_span(max_other_mean);
  output.spec = out_span(spec);
  output.focal_other_ratio = out_span(focal_other_ratio);
  output.is_contaminated = out_span(is_contaminated);
  output.pct_celltype = out_span(pct_celltype);
  output.pct_state = out_span(pct_state);
  output.pct_responder = out_span(pct_responder);
  output.pct_spill = out_span(pct_spill);
  output.pct_residual = out_span(pct_residual);
  output.keep = out_span(keep);
  const pace::Status status = pace::variance_decomposition(input, output, n_threads, user_interrupted);
  raise_if_failed(status, "variance decomposition");
  return Rcpp::List::create(
      Rcpp::Named("V_state_baseline") = v_state_baseline,
      Rcpp::Named("V_state_responder") = v_state_responder, Rcpp::Named("V_spill") = v_spill,
      Rcpp::Named("V_disp") = v_disp, Rcpp::Named("celltype_offset_sq") = celltype_offset_sq,
      Rcpp::Named("focal_mean") = focal_mean, Rcpp::Named("max_other_mean") = max_other_mean,
      Rcpp::Named("spec") = spec, Rcpp::Named("focal_other_ratio") = focal_other_ratio,
      Rcpp::Named("is_contaminated") = is_contaminated, Rcpp::Named("pct_celltype") = pct_celltype,
      Rcpp::Named("pct_state") = pct_state, Rcpp::Named("pct_responder") = pct_responder,
      Rcpp::Named("pct_spill") = pct_spill, Rcpp::Named("pct_residual") = pct_residual,
      Rcpp::Named("keep") = keep);
}

// The aggregate tables of the decomposition (see pace::decomposition_aggregates).
// `focal_code` is 0-based, in the row order the caller wants the tables in.
// [[Rcpp::export]]
Rcpp::List pace_decomposition_aggregates_cpp(const Rcpp::IntegerVector& focal_code, int n_focals,
                                             const Rcpp::NumericMatrix& pct,
                                             const Rcpp::NumericMatrix& components,
                                             const Rcpp::NumericVector& spec,
                                             const Rcpp::IntegerVector& is_contaminated,
                                             bool weight_by_spec_sq) {
  const std::int64_t n_rows = focal_code.size();
  Rcpp::NumericMatrix mean_table(n_focals, 5);
  Rcpp::NumericMatrix specw_table(weight_by_spec_sq ? n_focals : 0, weight_by_spec_sq ? 5 : 0);
  Rcpp::NumericMatrix pooled_table(n_focals, 5);
  Rcpp::NumericVector total(n_rows);
  Rcpp::NumericVector total_ss(n_focals);
  Rcpp::IntegerVector n_genes(n_focals), n_specific(n_focals);
  const pace::Status status = pace::decomposition_aggregates(
      int_span(focal_code), n_focals, n_rows, const_span(pct), const_span(components),
      double_span(spec), int_span(is_contaminated), out_span(mean_table), out_span(specw_table),
      out_span(pooled_table), out_span(total), out_span(total_ss), out_span(n_genes),
      out_span(n_specific));
  raise_if_failed(status, "decomposition aggregates");
  return Rcpp::List::create(Rcpp::Named("mean") = mean_table, Rcpp::Named("specw") = specw_table,
                            Rcpp::Named("pooled") = pooled_table, Rcpp::Named("Total") = total,
                            Rcpp::Named("total_SS") = total_ss, Rcpp::Named("n_genes") = n_genes,
                            Rcpp::Named("n_specific_genes") = n_specific);
}

// The four-block view of the decomposition (see pace::four_block_shares).
// [[Rcpp::export]]
Rcpp::List pace_four_block_shares_cpp(const Rcpp::NumericVector& celltype_offset_sq,
                                      const Rcpp::NumericVector& v_state_baseline,
                                      const Rcpp::NumericVector& v_state_responder,
                                      const Rcpp::NumericVector& v_spill,
                                      const Rcpp::NumericVector& v_disp) {
  const std::int64_t n_rows = celltype_offset_sq.size();
  Rcpp::NumericVector total(n_rows), v_state(n_rows), pct_celltype(n_rows), pct_state(n_rows);
  Rcpp::NumericVector pct_spill(n_rows), pct_residual(n_rows);
  const pace::Status status = pace::four_block_shares(
      double_span(celltype_offset_sq), double_span(v_state_baseline), double_span(v_state_responder),
      double_span(v_spill), double_span(v_disp), n_rows, out_span(total), out_span(v_state),
      out_span(pct_celltype), out_span(pct_state), out_span(pct_spill), out_span(pct_residual));
  raise_if_failed(status, "four-block shares");
  return Rcpp::List::create(Rcpp::Named("Total") = total, Rcpp::Named("V_state") = v_state,
                            Rcpp::Named("pct_celltype") = pct_celltype,
                            Rcpp::Named("pct_state") = pct_state,
                            Rcpp::Named("pct_spill") = pct_spill,
                            Rcpp::Named("pct_residual") = pct_residual);
}

// The observed single-frame shares (see pace::single_frame_shares). The four V
// matrices are n_groups x n_genes, aligned to the same genes as the statistics.
// Returned vectors are n_groups * n_genes, group-major.
// [[Rcpp::export]]
Rcpp::List pace_single_frame_shares_cpp(const Rcpp::NumericMatrix& focal_mean,
                                        const Rcpp::NumericMatrix& within_ss,
                                        const Rcpp::NumericVector& global_mean,
                                        const Rcpp::IntegerVector& group_size,
                                        const Rcpp::NumericMatrix& v_state,
                                        const Rcpp::NumericMatrix& v_responder,
                                        const Rcpp::NumericMatrix& v_spill,
                                        const Rcpp::NumericMatrix& v_disp, int n_threads) {
  const int n_groups = group_size.size();
  const std::int64_t n_genes = focal_mean.ncol();
  const R_xlen_t rows = static_cast<R_xlen_t>(n_groups) * n_genes;
  const bool has_responder = v_responder.size() != 0;
  Rcpp::NumericVector pct_celltype(rows), pct_spatial(rows), pct_spill(rows), pct_residual(rows);
  Rcpp::NumericVector pct_responder(has_responder ? rows : 0);
  Rcpp::NumericVector ss_lineage(rows), ss_within(rows), denom(rows);
  Rcpp::IntegerVector keep(n_groups);
  const pace::Status status = pace::single_frame_shares(
      const_span(focal_mean), const_span(within_ss), double_span(global_mean), int_span(group_size),
      n_groups, n_genes, const_span(v_state), const_span(v_responder), const_span(v_spill),
      const_span(v_disp), out_span(pct_celltype), out_span(pct_spatial), out_span(pct_responder),
      out_span(pct_spill), out_span(pct_residual), out_span(ss_lineage), out_span(ss_within),
      out_span(denom), out_span(keep), n_threads, user_interrupted);
  raise_if_failed(status, "single-frame shares");
  return Rcpp::List::create(Rcpp::Named("pct_celltype") = pct_celltype,
                            Rcpp::Named("pct_spatial") = pct_spatial,
                            Rcpp::Named("pct_responder") = pct_responder,
                            Rcpp::Named("pct_spill") = pct_spill,
                            Rcpp::Named("pct_residual") = pct_residual,
                            Rcpp::Named("SS_lineage") = ss_lineage,
                            Rcpp::Named("SS_within") = ss_within, Rcpp::Named("denom") = denom,
                            Rcpp::Named("keep") = keep);
}

// Canonical mash for one neighbour slice, end to end: covariances, grid,
// likelihood, EM weights and posteriors (see src/core/mash.*). Reproduces
// mashr's call sets exactly on both canonical cohorts, with no mashr dependency.
//
// `bhat` and `shat` are n_effects x n_conditions. Pass an empty `v` for the
// identity, which is what `null_correlation = FALSE` means.
// [[Rcpp::export]]
Rcpp::List pace_mash_cpp(const Rcpp::NumericMatrix& bhat, const Rcpp::NumericMatrix& shat,
                         const Rcpp::NumericMatrix& v, double null_weight, double grid_mult,
                         std::int64_t max_iterations, double tolerance) {
  const std::int64_t n_effects = bhat.nrow();
  const std::int64_t n_conditions = bhat.ncol();
  if (shat.nrow() != bhat.nrow() || shat.ncol() != bhat.ncol()) {
    Rcpp::stop("bhat and shat must have the same dimensions");
  }

  std::vector<std::string> condition_names;
  const Rcpp::List dimnames = bhat.attr("dimnames");
  if (dimnames.size() == 2 && !Rf_isNull(dimnames[1])) {
    const Rcpp::CharacterVector columns = dimnames[1];
    for (R_xlen_t i = 0; i < columns.size(); ++i) condition_names.push_back(Rcpp::as<std::string>(columns[i]));
  }

  const std::vector<pace::CovarianceComponent> covariances =
      pace::canonical_covariances(n_conditions, condition_names);

  std::vector<double> grid;
  raise_if_failed(pace::autoselect_grid(double_span(bhat), double_span(shat), n_effects,
                                        n_conditions, grid_mult, &grid),
                  "mash grid");

  const std::int64_t n_components = pace::n_mixture_components(
      static_cast<std::int64_t>(covariances.size()), static_cast<std::int64_t>(grid.size()), true);
  Rcpp::NumericMatrix loglik(n_effects, n_components);
  raise_if_failed(pace::mixture_loglik(double_span(bhat), double_span(shat), double_span(v),
                                       n_effects, n_conditions, covariances, grid, true, 1,
                                       out_span(loglik)),
                  "mash likelihood");

  std::vector<double> weights;
  std::int64_t iterations = 0;
  raise_if_failed(pace::optimize_mixture_weights(double_span(loglik), n_effects, n_components,
                                                 null_weight, max_iterations, tolerance, &weights,
                                                 &iterations),
                  "mash mixture weights");

  Rcpp::NumericMatrix posterior_mean(n_effects, n_conditions);
  Rcpp::NumericMatrix posterior_sd(n_effects, n_conditions);
  Rcpp::NumericMatrix lfsr(n_effects, n_conditions);
  raise_if_failed(pace::posterior_summaries(double_span(bhat), double_span(shat), double_span(v),
                                            n_effects, n_conditions, covariances, grid, true,
                                            weights, double_span(loglik), out_span(posterior_mean),
                                            out_span(posterior_sd), out_span(lfsr)),
                  "mash posteriors");

  posterior_mean.attr("dimnames") = bhat.attr("dimnames");
  posterior_sd.attr("dimnames") = bhat.attr("dimnames");
  lfsr.attr("dimnames") = bhat.attr("dimnames");
  return Rcpp::List::create(Rcpp::Named("PosteriorMean") = posterior_mean,
                            Rcpp::Named("PosteriorSD") = posterior_sd,
                            Rcpp::Named("lfsr") = lfsr,
                            Rcpp::Named("pi") = Rcpp::wrap(weights),
                            Rcpp::Named("grid") = Rcpp::wrap(grid),
                            Rcpp::Named("iterations") = iterations);
}

// Driver scores for one (focal, neighbour) pair (see pace::driver_scores).
// [[Rcpp::export]]
Rcpp::List pace_driver_scores_cpp(const Rcpp::NumericVector& estimate_shrunk,
                                  const Rcpp::NumericVector& spec,
                                  const Rcpp::NumericVector& focal_mean,
                                  const Rcpp::NumericVector& mu_bar,
                                  const Rcpp::NumericVector& alpha,
                                  const Rcpp::NumericVector& u_raw, double var_n, double var_rn,
                                  bool has_responder) {
  const std::int64_t n = estimate_shrunk.size();
  Rcpp::NumericVector mcsd(n), mcsd4(n), v_resid(n), v_s(n), v_total(n), r2_s(n);
  Rcpp::NumericVector v_rxs(has_responder ? n : 0), r2_rxs(has_responder ? n : 0);
  const pace::Status status = pace::driver_scores(
      double_span(estimate_shrunk), double_span(spec), double_span(focal_mean), double_span(mu_bar),
      double_span(alpha), double_span(u_raw), var_n, var_rn, has_responder, n, out_span(mcsd),
      out_span(mcsd4), out_span(v_resid), out_span(v_s), out_span(v_rxs), out_span(v_total),
      out_span(r2_s), out_span(r2_rxs));
  raise_if_failed(status, "driver scores");
  return Rcpp::List::create(Rcpp::Named("MCSD") = mcsd, Rcpp::Named("MCSD4") = mcsd4,
                            Rcpp::Named("V_resid") = v_resid, Rcpp::Named("V_S") = v_s,
                            Rcpp::Named("V_RxS") = v_rxs, Rcpp::Named("V_total") = v_total,
                            Rcpp::Named("R2_S") = r2_s, Rcpp::Named("R2_RxS") = r2_rxs);
}

// stats::cov() of selected columns of a sparse matrix over selected rows
// (see pace::subset_covariance). `rows` and `cols` are 1-based.
// [[Rcpp::export]]
Rcpp::NumericMatrix pace_subset_covariance_cpp(const Rcpp::S4& matrix, const Rcpp::IntegerVector& rows,
                                               const Rcpp::IntegerVector& cols,
                                               const Rcpp::NumericVector& scale) {
  CscHolder holder(matrix);
  Rcpp::NumericMatrix covariance(cols.size(), cols.size());
  const pace::Status status = pace::subset_covariance(holder.view, int_span(rows), int_span(cols),
                                                      double_span(scale), out_span(covariance));
  raise_if_failed(status, "subset covariance");
  return covariance;
}

// The values of one (1-based) column of a sparse matrix at selected (1-based) rows.
// [[Rcpp::export]]
Rcpp::NumericVector pace_subset_column_cpp(const Rcpp::S4& matrix, const Rcpp::IntegerVector& rows,
                                           int column) {
  CscHolder holder(matrix);
  Rcpp::NumericVector values(rows.size());
  const pace::Status status = pace::subset_column(holder.view, int_span(rows), column, out_span(values));
  raise_if_failed(status, "subset column");
  return values;
}

// The 1-based rows at which one column of a sparse matrix is not zero.
// [[Rcpp::export]]
Rcpp::IntegerVector pace_column_nonzero_rows_cpp(const Rcpp::S4& matrix, int column) {
  CscHolder holder(matrix);
  std::int64_t n_found = 0;
  pace::Status status = pace::column_nonzero_rows(holder.view, column, pace::Span<int>(), &n_found);
  raise_if_failed(status, "column non-zero rows");
  Rcpp::IntegerVector rows(n_found);
  status = pace::column_nonzero_rows(holder.view, column, out_span(rows), &n_found);
  raise_if_failed(status, "column non-zero rows");
  return rows;
}

// Pratt's pair attribution for one focal cell type (see pace::pair_variance_pratt).
// [[Rcpp::export]]
Rcpp::List pace_pair_variance_pratt_cpp(const Rcpp::NumericMatrix& sigma,
                                        const Rcpp::NumericMatrix& u) {
  const std::int64_t n_pairs = sigma.nrow();
  const std::int64_t n_genes = u.ncol();
  Rcpp::NumericVector v_pair(n_pairs), v_diag(n_pairs), share_pct(n_pairs), diag_pct(n_pairs);
  Rcpp::NumericVector sign(n_pairs);
  double v_total = 0;
  double v_diag_total = 0;
  double cross_cov_pct = 0;
  int n_negative = 0;
  const pace::Status status = pace::pair_variance_pratt(
      const_span(sigma), const_span(u), n_pairs, n_genes, out_span(v_pair), out_span(v_diag),
      out_span(share_pct), out_span(diag_pct), out_span(sign), &v_total, &v_diag_total,
      &cross_cov_pct, &n_negative);
  raise_if_failed(status, "pair variance");
  return Rcpp::List::create(Rcpp::Named("V_pair_pratt") = v_pair, Rcpp::Named("V_pair_diag") = v_diag,
                            Rcpp::Named("within_focal_share_pct") = share_pct,
                            Rcpp::Named("within_focal_diag_pct") = diag_pct,
                            Rcpp::Named("sign") = sign, Rcpp::Named("V_block_pratt") = v_total,
                            Rcpp::Named("V_block_diag") = v_diag_total,
                            Rcpp::Named("cross_cov_pct") = cross_cov_pct,
                            Rcpp::Named("n_negative_pairs") =
                                n_negative < 0 ? NA_INTEGER : n_negative);
}

// SS-weighted per-focal blocks of the single frame (see pace::single_frame_focal_blocks).
// [[Rcpp::export]]
Rcpp::NumericMatrix pace_single_frame_focal_blocks_cpp(const Rcpp::IntegerVector& focal_code,
                                                       int n_focals, const Rcpp::NumericMatrix& pct,
                                                       const Rcpp::NumericVector& denom) {
  Rcpp::NumericMatrix table(n_focals, pct.ncol());
  const pace::Status status = pace::single_frame_focal_blocks(
      int_span(focal_code), n_focals, focal_code.size(), const_span(pct), pct.ncol(),
      double_span(denom), out_span(table));
  raise_if_failed(status, "single-frame focal blocks");
  return table;
}

// True when every value is finite (see pace::all_finite): the same check as
// all(is.finite(x)) without allocating the logical vector.
// [[Rcpp::export]]
bool pace_all_finite_cpp(const Rcpp::NumericVector& values) {
  return pace::all_finite(double_span(values));
}

// ---------------------------------------------------------------------------
// The streaming solver's per-chunk arithmetic (see core/irls_chunk.hpp) and its
// variance-component updates (see core/hyperparameters.hpp). `first_gene` is the
// 1-based first column of the chunk in the counts / ambient matrices.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List pace_working_response_cpp(const Rcpp::NumericMatrix& eta, const Rcpp::S4& counts,
                                     const Rcpp::S4& ambient, int first_gene, int n_genes,
                                     const Rcpp::NumericVector& offset,
                                     const Rcpp::NumericVector& rho,
                                     const Rcpp::NumericVector& alpha,
                                     const Rcpp::NumericVector& sample_weight, bool nb2,
                                     bool gaussian, bool seed_iteration, int n_cells,
                                     int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  Rcpp::NumericMatrix z(n_cells, n_genes);
  Rcpp::NumericMatrix w(n_cells, n_genes);
  Rcpp::NumericVector colsum_w(n_genes);
  const pace::Status status = pace::working_response(
      seed_iteration ? pace::Span<const double>() : const_span(eta), gene_block(count_holder, first_gene),
      gene_block(ambient_holder, first_gene), double_span(offset), double_span(rho),
      double_span(alpha), double_span(sample_weight), nb2, gaussian, seed_iteration, n_cells,
      n_genes, out_span(z), out_span(w), out_span(colsum_w), n_threads, user_interrupted);
  raise_if_failed(status, "working response");
  return Rcpp::List::create(Rcpp::Named("z") = z, Rcpp::Named("w") = w,
                            Rcpp::Named("colsum_w") = colsum_w);
}

// [[Rcpp::export]]
Rcpp::List pace_rho_accumulate_cpp(const Rcpp::NumericMatrix& eta,
                                   const Rcpp::NumericMatrix& prev_eta, const Rcpp::S4& counts,
                                   const Rcpp::S4& ambient, int first_gene,
                                   const Rcpp::NumericVector& offset,
                                   const Rcpp::NumericVector& rho,
                                   const Rcpp::NumericVector& alpha,
                                   const Rcpp::NumericMatrix& mask,
                                   const Rcpp::IntegerVector& mask_index,
                                   const Rcpp::NumericVector& num_in,
                                   const Rcpp::NumericVector& den_in, bool nb2,
                                   bool seed_iteration, bool seed_previous, int n_cells,
                                   int n_genes_in, bool want_tail_counts, int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const std::int64_t n = seed_iteration ? n_cells : eta.nrow();
  const std::int64_t n_genes = seed_iteration ? n_genes_in : eta.ncol();
  Rcpp::NumericVector num(Rcpp::clone(num_in));
  Rcpp::NumericVector den(Rcpp::clone(den_in));
  double rel_delta_max = 0;
  double rel_delta_sum = 0;
  std::int64_t n_finite = 0;
  std::int64_t n_nonfinite = 0;
  Rcpp::NumericVector tail_counts(want_tail_counts ? 4 : 0);
  const pace::Status status = pace::rho_accumulate(
      const_span(eta), const_span(prev_eta), gene_block(count_holder, first_gene),
      gene_block(ambient_holder, first_gene), double_span(offset), double_span(rho),
      double_span(alpha), const_span(mask), int_span(mask_index), mask.nrow(), nb2, seed_iteration,
      seed_previous, n, n_genes, out_span(num), out_span(den), &rel_delta_max, &rel_delta_sum, &n_finite,
      &n_nonfinite, out_span(tail_counts), n_threads, user_interrupted);
  raise_if_failed(status, "rho accumulation");
  return Rcpp::List::create(Rcpp::Named("num") = num, Rcpp::Named("den") = den,
                            Rcpp::Named("rel_delta_max") = rel_delta_max,
                            Rcpp::Named("rel_delta_sum") = rel_delta_sum,
                            Rcpp::Named("n_finite") = static_cast<double>(n_finite),
                            Rcpp::Named("n_nonfinite") = static_cast<double>(n_nonfinite),
                            Rcpp::Named("tail_counts") = tail_counts);
}

// [[Rcpp::export]]
Rcpp::List pace_fitted_mean_column_cpp(const Rcpp::NumericVector& eta_column,
                                       const Rcpp::S4& counts, const Rcpp::S4& ambient,
                                       int first_gene, int gene,
                                       const Rcpp::NumericVector& offset,
                                       const Rcpp::NumericVector& rho) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const std::int64_t n = eta_column.size();
  Rcpp::NumericVector mu(n);
  Rcpp::NumericVector y(n);
  const pace::Status status = pace::fitted_mean_column(
      double_span(eta_column), gene_block(count_holder, first_gene),
      gene_block(ambient_holder, first_gene), double_span(offset), double_span(rho), n, gene - 1,
      out_span(mu), out_span(y));
  raise_if_failed(status, "fitted mean column");
  return Rcpp::List::create(Rcpp::Named("mu") = mu, Rcpp::Named("y") = y);
}

// [[Rcpp::export]]
Rcpp::List pace_rho_shrink_cpp(const Rcpp::NumericVector& num, const Rcpp::NumericVector& den) {
  const std::int64_t n = num.size();
  Rcpp::NumericVector rho(n);
  std::int64_t n_nonfinite = 0;
  double rho_bar = 0;
  double den_quantile = 0;
  const pace::Status status = pace::rho_shrink(double_span(num), double_span(den), n,
                                               out_span(rho), &n_nonfinite, &rho_bar, &den_quantile);
  raise_if_failed(status, "rho shrink");
  return Rcpp::List::create(Rcpp::Named("rho") = rho,
                            Rcpp::Named("n_nonfinite") = static_cast<double>(n_nonfinite),
                            Rcpp::Named("rho_bar") = rho_bar,
                            Rcpp::Named("den_quantile") = den_quantile);
}

// [[Rcpp::export]]
Rcpp::NumericVector pace_tau_em_update_cpp(const Rcpp::NumericMatrix& u,
                                           const Rcpp::NumericMatrix& re_var) {
  Rcpp::NumericVector tau(u.nrow());
  const pace::Status status = pace::tau_em_update(const_span(u), const_span(re_var), u.nrow(),
                                                  u.ncol(), out_span(tau));
  raise_if_failed(status, "tau EM update");
  return tau;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pace_tau_hierarchical_cpp(const Rcpp::NumericMatrix& tau,
                                              const Rcpp::IntegerVector& n_group,
                                              double lambda_factor) {
  Rcpp::NumericMatrix out(tau.nrow(), tau.ncol());
  const pace::Status status = pace::tau_hierarchical(const_span(tau), int_span(n_group), tau.nrow(),
                                                     tau.ncol(), lambda_factor, out_span(out));
  raise_if_failed(status, "hierarchical tau");
  return out;
}

// [[Rcpp::export]]
Rcpp::List pace_tau_eb_summaries_cpp(const Rcpp::NumericMatrix& s2, double reml_factor) {
  const std::int64_t q = s2.nrow();
  Rcpp::NumericVector panel(q), panel_median(q), log_variance(q);
  Rcpp::IntegerVector n_log_finite(q);
  const pace::Status status = pace::tau_eb_summaries(const_span(s2), q, s2.ncol(), reml_factor,
                                                     out_span(panel), out_span(panel_median),
                                                     out_span(log_variance), out_span(n_log_finite));
  raise_if_failed(status, "tau EB summaries");
  return Rcpp::List::create(Rcpp::Named("panel") = panel,
                            Rcpp::Named("panel_median") = panel_median,
                            Rcpp::Named("log_variance") = log_variance,
                            Rcpp::Named("n_log_finite") = n_log_finite);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pace_tau_eb_apply_cpp(const Rcpp::NumericMatrix& s2, double reml_factor,
                                          const Rcpp::NumericVector& panel,
                                          const Rcpp::NumericVector& panel_median,
                                          const Rcpp::NumericVector& d0, double tau_floor) {
  Rcpp::NumericMatrix out(s2.nrow(), s2.ncol());
  const pace::Status status = pace::tau_eb_apply(const_span(s2), s2.nrow(), s2.ncol(), reml_factor,
                                                 double_span(panel), double_span(panel_median),
                                                 double_span(d0), tau_floor, out_span(out));
  raise_if_failed(status, "tau EB shrinkage");
  return out;
}

// [[Rcpp::export]]
Rcpp::List pace_tau_half_cauchy_cpp(const Rcpp::NumericMatrix& s2, int n_em_iter, double tau_floor,
                                    const Rcpp::NumericVector& lambda2_prev,
                                    const Rcpp::NumericMatrix& a_prev) {
  Rcpp::NumericMatrix out(s2.nrow(), s2.ncol());
  Rcpp::NumericMatrix a_out(s2.nrow(), s2.ncol());
  Rcpp::NumericVector lambda2_out(s2.nrow());
  Rcpp::NumericVector panel(s2.nrow());
  const pace::Status status = pace::tau_half_cauchy(
      const_span(s2), s2.nrow(), s2.ncol(), n_em_iter, tau_floor, double_span(lambda2_prev),
      const_span(a_prev), out_span(out), out_span(lambda2_out), out_span(a_out), out_span(panel));
  raise_if_failed(status, "half-Cauchy tau");
  return Rcpp::List::create(Rcpp::Named("tau") = out, Rcpp::Named("lambda_sq") = lambda2_out,
                            Rcpp::Named("a") = a_out, Rcpp::Named("panel") = panel);
}

// [[Rcpp::export]]
Rcpp::List pace_tau_clamp_cpp(const Rcpp::NumericMatrix& tau, double tau_max) {
  Rcpp::NumericMatrix out(Rcpp::clone(tau));
  std::int64_t n_binding = 0;
  const pace::Status status = pace::tau_clamp(out_span(out), Rf_xlength(out), tau_max, &n_binding);
  raise_if_failed(status, "tau clamp");
  return Rcpp::List::create(Rcpp::Named("tau") = out,
                            Rcpp::Named("n_binding") = static_cast<double>(n_binding));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pace_data_informed_weights_cpp(const Rcpp::NumericMatrix& detection_rate,
                                                   const Rcpp::IntegerVector& focal_of_row,
                                                   const Rcpp::NumericVector& scale_of_row) {
  const std::int64_t q = focal_of_row.size();
  const std::int64_t n_genes = detection_rate.ncol();
  Rcpp::NumericMatrix weights(q, n_genes);
  const pace::Status status = pace::data_informed_weights(
      const_span(detection_rate), detection_rate.nrow(), n_genes, int_span(focal_of_row),
      double_span(scale_of_row), q, out_span(weights));
  raise_if_failed(status, "data-informed weights");
  return weights;
}

// The random-effect design, converted out of R once. The per-chunk solve
// binding used to rebuild this on EVERY call, and `group_cells` alone is one int
// per cell per block -- 9.8 MB at 1.2M cells, rebuilt 39 chunks x 13 iterations.
// Holding the R objects here keeps the spans the core sees alive for as long as
// the design does.
struct SolveDesign {
  std::vector<Rcpp::NumericMatrix> terms;
  std::vector<Rcpp::IntegerVector> cell_group;
  std::vector<std::vector<int>> group_offsets;
  std::vector<std::vector<int>> group_cells;
  std::vector<pace::SolveBlock> blocks;
};

SolveDesign build_solve_design(const Rcpp::List& blocks, const Rcpp::List& terms_list,
                               const Rcpp::List& cells_by_group_list,
                               const Rcpp::List& cell_group_list) {
  const int n_blocks = blocks.size();
  SolveDesign design;
  design.terms.resize(n_blocks);
  design.cell_group.resize(n_blocks);
  design.group_offsets.resize(n_blocks);
  design.group_cells.resize(n_blocks);
  design.blocks.resize(n_blocks);
  for (int b = 0; b < n_blocks; ++b) {
    const Rcpp::List block = blocks[b];
    design.blocks[b].col_offset = Rcpp::as<int>(block["col_offset"]);
    design.blocks[b].n_terms = Rcpp::as<int>(block["K_terms"]);
    design.blocks[b].n_groups = Rcpp::as<int>(block["K_groups"]);
    design.terms[b] = Rcpp::as<Rcpp::NumericMatrix>(terms_list[b]);
    design.cell_group[b] = Rcpp::as<Rcpp::IntegerVector>(cell_group_list[b]);
    const Rcpp::List cells = cells_by_group_list[b];
    design.group_offsets[b].reserve(design.blocks[b].n_groups + 1);
    design.group_offsets[b].push_back(0);
    for (int g = 0; g < design.blocks[b].n_groups; ++g) {
      const Rcpp::IntegerVector rows = cells[g];
      for (R_xlen_t i = 0; i < rows.size(); ++i) design.group_cells[b].push_back(rows[i] - 1);
      design.group_offsets[b].push_back(static_cast<int>(design.group_cells[b].size()));
    }
  }
  // The spans are set after every push_back is done, so no reallocation can
  // leave one dangling.
  for (int b = 0; b < n_blocks; ++b) {
    design.blocks[b].terms = const_span(design.terms[b]);
    design.blocks[b].cell_group = int_span(design.cell_group[b]);
    design.blocks[b].group_offsets = pace::Span<const int>(
        design.group_offsets[b].data(), static_cast<std::int64_t>(design.group_offsets[b].size()));
    design.blocks[b].group_cells = pace::Span<const int>(
        design.group_cells[b].data(), static_cast<std::int64_t>(design.group_cells[b].size()));
  }
  return design;
}

// The per-gene penalised WLS solve of one chunk (see core/gene_solve.hpp).
// `blocks` carries col_offset / K_terms / K_groups; `cells_by_grp_list` and
// `cell_grp_list` are 1-based, as R built them.
// [[Rcpp::export]]
Rcpp::List pace_solve_genes_chunk_cpp(const Rcpp::NumericMatrix& x_fixed,
                                      const Rcpp::NumericMatrix& w, const Rcpp::NumericMatrix& z,
                                      const Rcpp::NumericMatrix& lam_diag, int q_total,
                                      const Rcpp::List& blocks, const Rcpp::List& terms_list,
                                      const Rcpp::List& cells_by_group_list,
                                      const Rcpp::List& cell_group_list, bool single_precision,
                                      int n_threads) {
  const std::int64_t n = x_fixed.nrow();
  const std::int64_t p = x_fixed.ncol();
  const std::int64_t n_genes = w.ncol();
  const SolveDesign design =
      build_solve_design(blocks, terms_list, cells_by_group_list, cell_group_list);
  const std::vector<pace::SolveBlock>& core_blocks = design.blocks;

  Rcpp::NumericMatrix beta(p, n_genes);
  Rcpp::NumericMatrix u(q_total, n_genes);
  Rcpp::NumericMatrix ainv_diag(p + q_total, n_genes);
  const pace::Status status = pace::solve_genes_chunk(
      const_span(x_fixed), n, p, core_blocks, const_span(w), const_span(z), const_span(lam_diag),
      q_total, n_genes, single_precision, n_threads, user_interrupted, out_span(beta), out_span(u),
      out_span(ainv_diag));
  raise_if_failed(status, "gene solve");
  return Rcpp::List::create(Rcpp::Named("B") = beta, Rcpp::Named("U") = u,
                            Rcpp::Named("Ainv_diag") = ainv_diag);
}

// ---------------------------------------------------------------------------
// The dispersion MLE and the prior degrees of freedom (see core/dispersion.hpp).
// The core owns the optimisers and the likelihood; the two elementary special
// functions are R's own, passed in, so the numbers are identical to the
// optimize()/dnbinom() and uniroot()/trigamma() calls they replace while the
// core keeps no R dependency of its own.
// ---------------------------------------------------------------------------

namespace {

double r_log_nbinom(double x, double size, double mu) { return Rf_dnbinom_mu(x, size, mu, 1); }

double r_trigamma(double x) { return Rf_trigamma(x); }

}  // namespace

// [[Rcpp::export]]
Rcpp::List pace_dispersion_chunk_cpp(const Rcpp::NumericMatrix& eta, const Rcpp::S4& counts,
                                     const Rcpp::S4& ambient, int first_gene,
                                     const Rcpp::NumericVector& offset,
                                     const Rcpp::NumericVector& rho, bool nb2,
                                     bool zero_collapse, double max_cells, bool fast_density,
                                     int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const std::int64_t n = eta.nrow();
  const std::int64_t n_genes = eta.ncol();
  Rcpp::NumericVector alpha(n_genes);
  std::int64_t n_noninteger = 0;
  const pace::Status status = pace::dispersion_chunk(
      const_span(eta), gene_block(count_holder, first_gene), gene_block(ambient_holder, first_gene),
      double_span(offset), double_span(rho), nb2, zero_collapse, max_cells, r_log_nbinom,
      fast_density, n, n_genes, out_span(alpha), &n_noninteger, n_threads, user_interrupted);
  raise_if_failed(status, "dispersion");
  for (R_xlen_t j = 0; j < alpha.size(); ++j) {
    if (ISNAN(alpha[j])) alpha[j] = NA_REAL;
  }
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("n_noninteger") = static_cast<double>(n_noninteger));
}

// The same MLE for one gene, from its counts and fitted means.
// [[Rcpp::export]]
double pace_dispersion_mle_cpp(const Rcpp::NumericVector& counts, const Rcpp::NumericVector& mu,
                               bool nb2, bool zero_collapse, double max_cells,
                               bool fast_density) {
  const double alpha = pace::dispersion_mle(double_span(counts), double_span(mu), nb2,
                                            zero_collapse, max_cells, r_log_nbinom, fast_density);
  return ISNAN(alpha) ? NA_REAL : alpha;
}

// The prior degrees of freedom of the adaptive tau shrinkage.
// [[Rcpp::export]]
double pace_estimate_d0_cpp(double log_variance, double n_used, double d0_max) {
  return pace::estimate_d0(log_variance, static_cast<std::int64_t>(n_used), d0_max, r_trigamma);
}

// ---------------------------------------------------------------------------
// The kernel preprocessing and the anchor decision (see core/preprocess.hpp).
// Each returns a modified copy so the caller's matrix is left alone.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List pace_drop_sparse_kernel_cpp(const Rcpp::NumericMatrix& kernel,
                                       const Rcpp::IntegerVector& celltype, double min_effective) {
  Rcpp::NumericMatrix out(Rcpp::clone(kernel));
  std::int64_t n_dropped = 0;
  const pace::Status status = pace::drop_sparse_kernel(out_span(out), out.nrow(), out.ncol(),
                                                       int_span(celltype), min_effective,
                                                       &n_dropped);
  raise_if_failed(status, "kernel support");
  return Rcpp::List::create(Rcpp::Named("kernel") = out,
                            Rcpp::Named("n_dropped") = static_cast<double>(n_dropped));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pace_centre_within_groups_cpp(const Rcpp::NumericMatrix& kernel,
                                                  const Rcpp::IntegerVector& image,
                                                  const Rcpp::IntegerVector& celltype,
                                                  int n_images, int n_celltypes) {
  Rcpp::NumericMatrix out(Rcpp::clone(kernel));
  const pace::Status status = pace::centre_within_groups(out_span(out), out.nrow(), out.ncol(),
                                                          int_span(image), int_span(celltype),
                                                          n_images, n_celltypes);
  raise_if_failed(status, "within-group centring");
  return out;
}

// [[Rcpp::export]]
Rcpp::List pace_standardise_cpp(const Rcpp::NumericVector& values) {
  Rcpp::NumericVector out(Rcpp::clone(values));
  double spread = 0;
  const pace::Status status = pace::standardise_column(out_span(out), out.size(), &spread);
  raise_if_failed(status, "standardise");
  return Rcpp::List::create(Rcpp::Named("values") = out, Rcpp::Named("sd") = spread);
}

// [[Rcpp::export]]
Rcpp::List pace_anchor_mask_cpp(const Rcpp::NumericMatrix& core_means, double owner_threshold,
                                double core_threshold) {
  Rcpp::NumericMatrix mask(core_means.nrow(), core_means.ncol());
  Rcpp::IntegerVector n_anchor(core_means.nrow());
  const pace::Status status = pace::anchor_mask(const_span(core_means), core_means.nrow(),
                                                 core_means.ncol(), owner_threshold,
                                                 core_threshold, out_span(mask),
                                                 out_span(n_anchor));
  raise_if_failed(status, "anchor mask");
  return Rcpp::List::create(Rcpp::Named("mask") = mask, Rcpp::Named("n_anchor") = n_anchor);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pace_normalise_rows_to_max_cpp(const Rcpp::NumericMatrix& values) {
  Rcpp::NumericMatrix out(Rcpp::clone(values));
  const pace::Status status = pace::normalise_rows_to_max(out_span(out), out.nrow(), out.ncol());
  raise_if_failed(status, "row normalisation");
  return out;
}

// eta for one gene chunk: X_fixed %*% B[, genes] + Z %*% U[, genes], written
// straight into the returned matrix. Replaces the R closure .eta_block(), which
// allocated three dense n x chunk matrices per call.
//
// `genes` is 1-based, as it comes from R. `x1` is used when p == 1 and may be
// length 0 when it is all ones; `x_fixed` is the dense n x p design otherwise.
// [[Rcpp::export]]
Rcpp::NumericMatrix pace_eta_block_cpp(const Rcpp::NumericVector& x1, bool x1_is_unit,
                                       const Rcpp::NumericMatrix& x_fixed, int p,
                                       const Rcpp::NumericMatrix& b, const Rcpp::S4& z_design,
                                       const Rcpp::NumericMatrix& u,
                                       const Rcpp::IntegerVector& genes, int n_threads) {
  const Rcpp::IntegerVector dims = z_design.slot("Dim");
  const Rcpp::IntegerVector column_pointer = z_design.slot("p");
  const Rcpp::IntegerVector row_index = z_design.slot("i");
  const Rcpp::NumericVector values = z_design.slot("x");
  pace::CscView z;
  z.column_pointer = int_span(column_pointer);
  z.row_index = int_span(row_index);
  z.values = double_span(values);
  z.n_rows = dims[0];
  z.n_cols = dims[1];

  const std::int64_t n = dims[0];
  const std::int64_t n_chunk = genes.size();
  const std::int64_t n_genes_total = b.ncol();

  // R hands us 1-based gene numbers; the core works 0-based.
  std::vector<int> genes0(static_cast<std::size_t>(n_chunk));
  for (std::int64_t j = 0; j < n_chunk; ++j) genes0[static_cast<std::size_t>(j)] = genes[j] - 1;

  Rcpp::NumericMatrix eta(n, n_chunk);
  const pace::Status status = pace::eta_block(
      double_span(x1), x1_is_unit, double_span(x_fixed), p, double_span(b), z, nullptr,
      double_span(u), pace::Span<const int>(genes0.data(), n_chunk), n, n_genes_total,
      pace::Span<double>(eta.begin(), n * n_chunk), 0, n_threads, user_interrupted);
  raise_if_failed(status, "eta block");
  return eta;
}

// The same eta, taken through the block structure instead of through Z's
// compressed columns (see core/linear_predictor.hpp). Bit-identical to
// pace_eta_block_cpp() and much cheaper: the design is read once per gene chunk
// rather than once per gene. Exposed so that the equality can be tested
// directly and the cell panel measured, rather than only through a whole fit.
//
// `blocks` carries col_offset / K_terms / K_groups and `cell_group_list` is
// 1-based, both as build_random_design_multi() produced them. `cell_panel` of 0
// lets the core choose.
// [[Rcpp::export]]
Rcpp::NumericMatrix pace_eta_block_blocked_cpp(const Rcpp::NumericVector& x1, bool x1_is_unit,
                                               const Rcpp::NumericMatrix& x_fixed, int p,
                                               const Rcpp::NumericMatrix& b,
                                               const Rcpp::S4& z_design,
                                               const Rcpp::List& blocks,
                                               const Rcpp::List& terms_list,
                                               const Rcpp::List& cell_group_list,
                                               const Rcpp::NumericMatrix& u,
                                               const Rcpp::IntegerVector& genes, int cell_panel,
                                               int n_threads) {
  const CscHolder z_holder(z_design);
  const std::int64_t n = z_holder.view.n_rows;
  const std::int64_t n_chunk = genes.size();
  const std::int64_t n_genes_total = b.ncol();

  const int n_blocks = blocks.size();
  std::vector<Rcpp::NumericMatrix> terms(static_cast<std::size_t>(n_blocks));
  std::vector<Rcpp::IntegerVector> cell_group(static_cast<std::size_t>(n_blocks));
  std::vector<pace::SolveBlock> core_blocks(static_cast<std::size_t>(n_blocks));
  for (int block_index = 0; block_index < n_blocks; ++block_index) {
    const std::size_t slot = static_cast<std::size_t>(block_index);
    const Rcpp::List block = blocks[block_index];
    core_blocks[slot].col_offset = Rcpp::as<int>(block["col_offset"]);
    core_blocks[slot].n_terms = Rcpp::as<int>(block["K_terms"]);
    core_blocks[slot].n_groups = Rcpp::as<int>(block["K_groups"]);
    terms[slot] = Rcpp::as<Rcpp::NumericMatrix>(terms_list[block_index]);
    cell_group[slot] = Rcpp::as<Rcpp::IntegerVector>(cell_group_list[block_index]);
    core_blocks[slot].terms = const_span(terms[slot]);
    core_blocks[slot].cell_group = int_span(cell_group[slot]);
  }

  // R hands us 1-based gene numbers; the core works 0-based.
  std::vector<int> genes0(static_cast<std::size_t>(n_chunk));
  for (std::int64_t j = 0; j < n_chunk; ++j) genes0[static_cast<std::size_t>(j)] = genes[j] - 1;

  Rcpp::NumericMatrix eta(n, n_chunk);
  const pace::Status status = pace::eta_block(
      double_span(x1), x1_is_unit, double_span(x_fixed), p, double_span(b), z_holder.view,
      &core_blocks, double_span(u), pace::Span<const int>(genes0.data(), n_chunk), n,
      n_genes_total, pace::Span<double>(eta.begin(), n * n_chunk), cell_panel, n_threads,
      user_interrupted);
  raise_if_failed(status, "eta block (blocked)");
  return eta;
}

// The post-solve rho accumulation over every gene chunk, without returning to
// R between them. Per chunk R built TWO n x chunk eta matrices -- the current
// one and the previous iteration's -- handed them in, and took the running
// accumulators back to pass into the next chunk. Here both are buffers reused
// across chunks and the accumulators never leave.
// [[Rcpp::export]]
Rcpp::List pace_rho_pass_cpp(
    const Rcpp::NumericVector& x1, bool x1_is_unit, const Rcpp::NumericMatrix& x_fixed, int p,
    const Rcpp::S4& z_design, const Rcpp::NumericMatrix& b_in, const Rcpp::NumericMatrix& u_in,
    const Rcpp::NumericMatrix& prev_b, const Rcpp::NumericMatrix& prev_u, bool have_previous,
    const Rcpp::S4& counts, const Rcpp::S4& ambient, const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& rho, const Rcpp::NumericVector& alpha,
    const Rcpp::NumericMatrix& mask, const Rcpp::IntegerVector& mask_index, bool nb2,
    bool want_tail_counts, int n_cells, int chunk_size, int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const CscHolder z_holder(z_design);
  const std::int64_t n = n_cells;
  const std::int64_t n_genes_total = u_in.ncol();

  Rcpp::NumericVector num(n);
  Rcpp::NumericVector den(n);
  Rcpp::NumericVector tail_counts(want_tail_counts ? 4 : 0);
  double rel_delta_max = 0;
  double rel_delta_sum = 0;
  std::int64_t n_finite = 0;
  std::int64_t n_nonfinite = 0;
  const pace::Status status = pace::rho_pass(
      double_span(x1), x1_is_unit, const_span(x_fixed), p, z_holder.view, nullptr,
      const_span(b_in), const_span(u_in), const_span(prev_b), const_span(prev_u), have_previous,
      gene_block(count_holder, 1), cached_ambient(ambient_holder, 1), double_span(offset),
      double_span(rho), double_span(alpha), const_span(mask), int_span(mask_index), mask.nrow(),
      nb2, n, n_genes_total, chunk_size, out_span(num), out_span(den), &rel_delta_max,
      &rel_delta_sum, &n_finite, &n_nonfinite, out_span(tail_counts), n_threads,
      user_interrupted);
  raise_if_failed(status, "rho accumulation");
  return Rcpp::List::create(
      Rcpp::Named("num") = num, Rcpp::Named("den") = den,
      Rcpp::Named("rel_delta_max") = rel_delta_max,
      Rcpp::Named("rel_delta_sum") = rel_delta_sum,
      Rcpp::Named("n_finite") = static_cast<double>(n_finite),
      Rcpp::Named("n_nonfinite") = static_cast<double>(n_nonfinite),
      Rcpp::Named("tail_counts") = tail_counts);
}

// The dispersion MLE over every gene chunk, without returning to R between
// them. Per chunk R built eta as an n x chunk matrix, passed it in, took the
// alphas back and dropped the matrix. Here eta is a buffer reused across
// chunks and only the finished alphas cross back.
// The Gaussian path's sigma2 seed: the marginal per-gene variance, before there
// is a fit to take residuals from. Replaces colMeans(Y * Y) - colMeans(Y)^2 in
// R, which materialised a second n x G matrix to do it.
// [[Rcpp::export]]
Rcpp::NumericVector pace_marginal_variance_cpp(const Rcpp::S4& counts, double floor_value,
                                               int n_cells, int n_genes, int n_threads) {
  CscHolder count_holder(counts);
  Rcpp::NumericVector sigma2(n_genes);
  const pace::Status status = pace::marginal_variance(
      gene_block(count_holder, 1), floor_value, n_cells, n_genes, out_span(sigma2), n_threads,
      user_interrupted);
  raise_if_failed(status, "marginal variance");
  return sigma2;
}

// [[Rcpp::export]]
Rcpp::List pace_dispersion_pass_cpp(
    const Rcpp::NumericVector& x1, bool x1_is_unit, const Rcpp::NumericMatrix& x_fixed, int p,
    const Rcpp::S4& z_design, const Rcpp::NumericMatrix& b_in, const Rcpp::NumericMatrix& u_in,
    const Rcpp::S4& counts, const Rcpp::S4& ambient, const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& rho, bool nb2, bool gaussian, bool zero_collapse, double max_cells,
    bool fast_density, int n_cells, int chunk_size, int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const CscHolder z_holder(z_design);
  const std::int64_t n = n_cells;
  const std::int64_t n_genes_total = u_in.ncol();

  Rcpp::NumericVector alpha(n_genes_total);
  std::int64_t n_noninteger = 0;
  const pace::Status status = pace::dispersion_pass(
      double_span(x1), x1_is_unit, const_span(x_fixed), p, z_holder.view, nullptr,
      const_span(b_in), const_span(u_in), gene_block(count_holder, 1),
      cached_ambient(ambient_holder, 1),
      double_span(offset), double_span(rho), nb2, gaussian, zero_collapse, max_cells,
      r_log_nbinom, fast_density, n, n_genes_total, chunk_size, out_span(alpha), &n_noninteger,
      n_threads, user_interrupted);
  raise_if_failed(status, "dispersion");
  for (R_xlen_t j = 0; j < alpha.size(); ++j) {
    if (ISNAN(alpha[j])) alpha[j] = NA_REAL;
  }
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("n_noninteger") = static_cast<double>(n_noninteger));
}

// One logical chunk's working response and weights, built a sub-block at a time
// inside the core.
//
// This replaces an R loop that allocated two dense n x chunk matrices per chunk,
// called pace_working_response_cpp() per sub-block (two more n x sub_genes
// allocations each) and copied every sub-block's result into the chunk matrices.
// At 1.2M cells and a chunk of 64 that was 1.25 GB of R allocation per chunk,
// plus the per-sub-block temporaries and the copies. Here z and w are allocated
// once, as the return value, and each sub-block writes straight into its own
// columns.
//
// eta: when `eta_chunk` has columns it is used as-is -- the fused path already
// built it at chunk width, and a sub-block is then just a pointer into it, with
// no copy. Otherwise, and when not seeding, it is built here per sub-block into
// a scratch buffer, which is what keeps the memory bounded by sub_genes.
//
// `first_gene` and `genes` are 1-based, as they come from R. Sub-blocking does
// not change the result: working_response treats each gene independently.
// One IRLS pass over every gene chunk, without returning to R between them.
//
// This is the chunk loop that used to live in fit_pace_mvpql_streaming(). Per
// chunk R built the working response, handed z and w back, called the solve,
// took three matrices back, and wrote them into B and U. The design was
// reconverted on every one of those solve calls, and z, w and lam were R
// matrices allocated per chunk. Here the design is converted ONCE, z and w are
// plain buffers reused across chunks, and only the finished B, U, re_var and
// standard errors cross back.
//
// The arithmetic is the same calls in the same order, so the fit is unchanged.
// [[Rcpp::export]]
Rcpp::List pace_fit_pass1_cpp(
    const Rcpp::NumericVector& x1, bool x1_is_unit, const Rcpp::NumericMatrix& x_fixed,
    const Rcpp::NumericMatrix& solve_x_fixed, int p,
    const Rcpp::S4& z_design, const Rcpp::List& blocks, const Rcpp::List& terms_list,
    const Rcpp::List& cells_by_group_list, const Rcpp::List& cell_group_list,
    const Rcpp::NumericMatrix& b_in, const Rcpp::NumericMatrix& u_in,
    const Rcpp::NumericMatrix& lam_diag, const Rcpp::S4& counts, const Rcpp::S4& ambient,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& rho,
    const Rcpp::NumericVector& alpha, const Rcpp::NumericVector& sample_weight, bool nb2,
    bool gaussian, bool seed_iteration, int n_cells, int chunk_size, int sub_genes,
    int interior_precision, bool last_iter, int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const CscHolder z_holder(z_design);
  const std::int64_t n = n_cells;
  const std::int64_t q_total = u_in.nrow();
  const std::int64_t n_genes_total = u_in.ncol();

  // Converted once for the whole pass, not once per chunk.
  const SolveDesign design =
      build_solve_design(blocks, terms_list, cells_by_group_list, cell_group_list);

  Rcpp::NumericMatrix beta_out(p, n_genes_total);
  Rcpp::NumericMatrix u_out(q_total, n_genes_total);
  Rcpp::NumericMatrix re_var(q_total, n_genes_total);
  Rcpp::NumericMatrix se_beta(p, n_genes_total);
  Rcpp::NumericMatrix se_u(q_total, n_genes_total);
  std::vector<int> nan_genes;

  const pace::Status status = pace::fit_pass1(
      double_span(x1), x1_is_unit, const_span(x_fixed), const_span(solve_x_fixed), p,
      z_holder.view, design.blocks, const_span(b_in), const_span(u_in), const_span(lam_diag),
      gene_block(count_holder, 1), cached_ambient(ambient_holder, 1), double_span(offset),
      double_span(rho), double_span(alpha), double_span(sample_weight), nb2, gaussian,
      seed_iteration, n, q_total, n_genes_total, chunk_size, sub_genes, interior_precision,
      last_iter, out_span(beta_out), out_span(u_out), out_span(re_var), out_span(se_beta),
      out_span(se_u), &nan_genes, n_threads, user_interrupted);
  raise_if_failed(status, "fit pass 1");

  return Rcpp::List::create(Rcpp::Named("B") = beta_out, Rcpp::Named("U") = u_out,
                            Rcpp::Named("re_var") = re_var, Rcpp::Named("se_B") = se_beta,
                            Rcpp::Named("se_U") = se_u,
                            Rcpp::Named("nan_genes") = Rcpp::wrap(nan_genes));
}

// [[Rcpp::export]]
Rcpp::List pace_working_response_chunk_cpp(
    const Rcpp::NumericMatrix& eta_chunk, const Rcpp::NumericVector& x1, bool x1_is_unit,
    const Rcpp::NumericMatrix& x_fixed, int p, const Rcpp::NumericMatrix& b,
    const Rcpp::S4& z_design, const Rcpp::NumericMatrix& u, const Rcpp::IntegerVector& genes,
    const Rcpp::S4& counts, const Rcpp::S4& ambient, int first_gene,
    const Rcpp::NumericVector& offset, const Rcpp::NumericVector& rho,
    const Rcpp::NumericVector& alpha, const Rcpp::NumericVector& sample_weight, bool nb2,
    bool gaussian, bool seed_iteration, int n_cells, int sub_genes, int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);

  const Rcpp::IntegerVector z_dims = z_design.slot("Dim");
  const Rcpp::IntegerVector z_column_pointer = z_design.slot("p");
  const Rcpp::IntegerVector z_row_index = z_design.slot("i");
  const Rcpp::NumericVector z_values = z_design.slot("x");
  pace::CscView z_view;
  z_view.column_pointer = int_span(z_column_pointer);
  z_view.row_index = int_span(z_row_index);
  z_view.values = double_span(z_values);
  z_view.n_rows = z_dims[0];
  z_view.n_cols = z_dims[1];

  const std::int64_t n = n_cells;
  const std::int64_t m_chunk = genes.size();
  Rcpp::NumericMatrix z(n, m_chunk);
  Rcpp::NumericMatrix w(n, m_chunk);
  Rcpp::NumericVector colsum_w(m_chunk);

  const bool have_eta = eta_chunk.ncol() > 0;
  std::vector<int> genes0(static_cast<std::size_t>(m_chunk));
  for (std::int64_t j = 0; j < m_chunk; ++j) genes0[static_cast<std::size_t>(j)] = genes[j] - 1;

  std::vector<double> eta_scratch;
  const std::int64_t step = sub_genes > 0 ? sub_genes : m_chunk;
  for (std::int64_t start = 0; start < m_chunk; start += step) {
    const std::int64_t len = std::min(step, m_chunk - start);
    pace::Span<const double> eta_span;
    if (!seed_iteration) {
      if (have_eta) {
        eta_span = pace::Span<const double>(eta_chunk.begin() + start * n, n * len);
      } else {
        eta_scratch.resize(static_cast<std::size_t>(n * len));
        const pace::Status eta_status = pace::eta_block(
            double_span(x1), x1_is_unit, double_span(x_fixed), p, double_span(b), z_view,
            nullptr, double_span(u), pace::Span<const int>(genes0.data() + start, len), n,
            b.ncol(), pace::Span<double>(eta_scratch.data(), n * len), 0, n_threads,
            user_interrupted);
        raise_if_failed(eta_status, "eta block");
        eta_span = pace::Span<const double>(eta_scratch.data(), n * len);
      }
    }
    const pace::Status status = pace::working_response(
        eta_span, gene_block(count_holder, first_gene + static_cast<int>(start)),
        gene_block(ambient_holder, first_gene + static_cast<int>(start)), double_span(offset),
        double_span(rho), pace::Span<const double>(alpha.begin() + start, len),
        double_span(sample_weight), nb2, gaussian, seed_iteration, n, len,
        pace::Span<double>(z.begin() + start * n, n * len),
        pace::Span<double>(w.begin() + start * n, n * len),
        pace::Span<double>(colsum_w.begin() + start, len), n_threads, user_interrupted);
    raise_if_failed(status, "working response");
  }
  return Rcpp::List::create(Rcpp::Named("z") = z, Rcpp::Named("w") = w,
                            Rcpp::Named("colsum_w") = colsum_w);
}

// The whole outer PQL iteration, run in the core (see core/irls_driver.hpp).
//
// This replaces the `for (it in seq_len(n_iter))` loop in
// fit_pace_mvpql_streaming(). That loop crossed the boundary once per pass per
// iteration -- about two hundred times on a 32-iteration fit -- and every
// crossing allocated a fresh R copy of each full-panel matrix the pass returned
// before R copied it into the loop's own state. Here the state stays in the
// core and only the finished fit and its history come back.
//
// Everything the loop would print goes out through Rprintf as it happens, so a
// verbose fit still streams its progress. Warnings are COLLECTED and raised
// after the loop returns rather than from inside it: an R warning can longjmp
// (under options(warn = 2) it becomes an error), and a longjmp out of the core
// would skip the destructors of every buffer the loop holds. Their text and
// their order are unchanged.
//
// `tau_shrinkage` is the pace::TauShrinkage enumerator as an integer, in the
// order fit_pace_mvpql_streaming()'s argument lists the modes. A non-finite
// `tau_max` means no cap, as .clamp_tau() read it.
// [[Rcpp::export]]
Rcpp::List pace_irls_driver_cpp(
    const Rcpp::NumericVector& x1, bool x1_is_unit, const Rcpp::NumericMatrix& x_fixed,
    const Rcpp::NumericMatrix& solve_x_fixed, int p, const Rcpp::S4& z_design,
    const Rcpp::List& blocks, const Rcpp::List& terms_list,
    const Rcpp::List& cells_by_group_list, const Rcpp::List& cell_group_list,
    const Rcpp::S4& counts, const Rcpp::S4& ambient, const Rcpp::NumericVector& offset,
    const Rcpp::NumericVector& sample_weight, const Rcpp::NumericMatrix& mask,
    const Rcpp::IntegerVector& mask_index, const Rcpp::NumericMatrix& data_informed_weights,
    const Rcpp::NumericVector& alpha_init, const Rcpp::CharacterVector& gene_names,
    bool ambient_streamed, int n_cells,
    int q_total, int n_genes, int n_iter, int min_iter, double early_stop_tol,
    double alpha_warmup, const std::string& alpha_warmup_label, bool nb2, bool gaussian,
    bool zero_collapse, double alpha_max_cells, bool fast_density, int interior_precision,
    int chunk_size, int sub_genes, double tau_max, int tau_shrinkage, double d0_min,
    bool use_reml, bool rd_diag, bool verbose, int n_threads) {
  CscHolder count_holder(counts);
  CscHolder ambient_holder(ambient);
  const CscHolder z_holder(z_design);
  const std::int64_t n = n_cells;
  const std::int64_t q = q_total;
  const std::int64_t g = n_genes;

  // The random-effect design, converted once for the whole fit rather than once
  // per chunk per iteration.
  const SolveDesign design =
      build_solve_design(blocks, terms_list, cells_by_group_list, cell_group_list);

  // The same blocks as the tau bookkeeping sees them, with the group sizes the
  // hierarchical weights read.
  const int n_blocks = blocks.size();
  std::vector<std::vector<int> > group_sizes(static_cast<std::size_t>(n_blocks));
  std::vector<pace::TauBlock> tau_blocks(static_cast<std::size_t>(n_blocks));
  for (int b = 0; b < n_blocks; ++b) {
    const Rcpp::List block = blocks[b];
    pace::TauBlock& tau_block = tau_blocks[static_cast<std::size_t>(b)];
    tau_block.col_offset = Rcpp::as<int>(block["col_offset"]);
    tau_block.n_terms = Rcpp::as<int>(block["K_terms"]);
    tau_block.n_groups = Rcpp::as<int>(block["K_groups"]);
    tau_block.n_cols = Rcpp::as<int>(block["n_cols"]);
    const Rcpp::List cells = cells_by_group_list[b];
    std::vector<int>& sizes = group_sizes[static_cast<std::size_t>(b)];
    sizes.reserve(static_cast<std::size_t>(cells.size()));
    for (R_xlen_t group = 0; group < cells.size(); ++group) {
      const Rcpp::IntegerVector rows = cells[group];
      sizes.push_back(static_cast<int>(rows.size()));
    }
  }
  // Set after every push_back is done, so no reallocation can leave one dangling.
  for (int b = 0; b < n_blocks; ++b) {
    const std::vector<int>& sizes = group_sizes[static_cast<std::size_t>(b)];
    tau_blocks[static_cast<std::size_t>(b)].group_size =
        pace::Span<const int>(sizes.data(), static_cast<std::int64_t>(sizes.size()));
  }

  Rcpp::NumericMatrix beta_out(p, g);
  Rcpp::NumericMatrix u_out(q, g);
  // The standard errors are filled on the last iteration only, so an unfinished
  // fit reports NA rather than zero, which is what the R loop's
  // matrix(NA_real_, ...) initialisation gave.
  Rcpp::NumericMatrix se_beta(p, g);
  Rcpp::NumericMatrix se_u(q, g);
  std::fill(se_beta.begin(), se_beta.end(), NA_REAL);
  std::fill(se_u.begin(), se_u.end(), NA_REAL);
  Rcpp::NumericVector alpha(g);
  Rcpp::NumericMatrix tau_g_array(q, g);
  Rcpp::NumericVector tau_flat(q);
  Rcpp::NumericVector rho(n);
  Rcpp::NumericMatrix history_tau_flat(q, n_iter);
  Rcpp::NumericMatrix history_alpha(g, n_iter);
  Rcpp::NumericVector history_rel_delta(n_iter);
  Rcpp::NumericVector history_rel_delta_mean(n_iter);
  Rcpp::NumericVector history_n_nan_genes(n_iter);
  Rcpp::NumericVector history_n_nonfinite(n_iter);
  Rcpp::NumericVector history_tau_max_seen(n_iter);
  std::int64_t iterations_run = 0;
  bool converged = false;

  std::vector<std::string> deferred_warnings;
  pace::IrlsLoopReporter reporter;
  reporter.message = [](const std::string& text) { Rprintf("%s", text.c_str()); };
  reporter.warn = [&deferred_warnings](const std::string& text) {
    deferred_warnings.push_back(text);
  };
  // The only message whose text is R's: it names the genes, and the gene names
  // are the column names of the counts.
  reporter.warn_nan_genes = [&deferred_warnings, &gene_names](std::int64_t iteration,
                                                              pace::Span<const int> genes) {
    std::string named;
    const std::int64_t shown = std::min<std::int64_t>(genes.size, 5);
    // k > 0, not !named.empty(): an empty string in colnames(Y) would otherwise
    // swallow the separator before the NEXT name, where R's paste(collapse=", ")
    // keeps it. The index is bounded too -- the old guard established only that
    // gene_names was non-empty, not that it covered every gene.
    for (std::int64_t k = 0; k < shown; ++k) {
      const std::int64_t gene = genes[k] - 1;
      if (gene < 0 || gene >= static_cast<std::int64_t>(gene_names.size())) continue;
      if (k > 0) named += ", ";
      named += Rcpp::as<std::string>(gene_names[gene]);
    }
    char buffer[1024];
    std::snprintf(buffer, sizeof(buffer),
                  "iter %lld: %lld gene(s) carry non-finite coefficients (%s%s); the fit will "
                  "not be reported as converged.",
                  static_cast<long long>(iteration), static_cast<long long>(genes.size),
                  named.c_str(), genes.size > 5 ? ", ..." : "");
    deferred_warnings.push_back(std::string(buffer));
  };

  pace::IrlsLoopInputs inputs;
  inputs.x1 = double_span(x1);
  inputs.x1_is_unit = x1_is_unit;
  inputs.x_fixed = const_span(x_fixed);
  inputs.solve_x_fixed = const_span(solve_x_fixed);
  inputs.p = p;
  inputs.z = z_holder.view;
  inputs.solve_blocks = &design.blocks;
  inputs.tau_blocks = &tau_blocks;
  inputs.counts = gene_block(count_holder, 1);
  // Cached or streamed. In streamed mode `ambient` is W (n x n) rather than the
  // product, and each chunk's columns are computed from W and the counts as the
  // fit walks the panel. The two give identical numbers -- a sparse product is
  // column-independent -- and differ only in holding 5.2 GB or not at the 1.2M
  // by 5,001 scale where this matters.
  std::vector<int> ambient_scratch_p, ambient_scratch_i;
  std::vector<double> ambient_scratch_x;
  pace::CscView ambient_scratch_view;
  if (ambient_streamed) {
    inputs.ambient.streamed = true;
    inputs.ambient.weights = ambient_holder.view;
    inputs.ambient.counts = count_holder.view;
    inputs.ambient.scratch_column_pointer = &ambient_scratch_p;
    inputs.ambient.scratch_row_index = &ambient_scratch_i;
    inputs.ambient.scratch_values = &ambient_scratch_x;
    inputs.ambient.scratch_view = &ambient_scratch_view;
  } else {
    inputs.ambient = cached_ambient(ambient_holder, 1);
  }
  inputs.offset = double_span(offset);
  inputs.sample_weight = double_span(sample_weight);
  inputs.mask = const_span(mask);
  inputs.mask_index = int_span(mask_index);
  inputs.n_mask_rows = mask.nrow();
  inputs.data_informed_weights = const_span(data_informed_weights);
  inputs.alpha_init = double_span(alpha_init);
  inputs.n = n;
  inputs.q = q;
  inputs.n_genes = g;

  pace::IrlsLoopOptions options;
  options.n_iter = n_iter;
  options.min_iter = min_iter;
  options.early_stop_tol = early_stop_tol;
  options.alpha_warmup = alpha_warmup;
  options.alpha_warmup_label = alpha_warmup_label;
  options.nb2 = nb2;
  options.gaussian = gaussian;
  options.zero_collapse = zero_collapse;
  options.alpha_max_cells = alpha_max_cells;
  options.fast_density = fast_density;
  options.interior_precision = interior_precision;
  options.chunk_size = chunk_size;
  options.sub_genes = sub_genes;
  options.tau_max = tau_max;
  options.tau_shrinkage = static_cast<pace::TauShrinkage>(tau_shrinkage);
  options.d0_min = d0_min;
  options.use_reml = use_reml;
  options.rd_diag = rd_diag;
  options.verbose = verbose;
  options.n_threads = n_threads;

  pace::IrlsLoopOutputs outputs;
  outputs.beta = out_span(beta_out);
  outputs.u = out_span(u_out);
  outputs.se_beta = out_span(se_beta);
  outputs.se_u = out_span(se_u);
  outputs.alpha = out_span(alpha);
  outputs.tau_g_array = out_span(tau_g_array);
  outputs.tau_flat = out_span(tau_flat);
  outputs.rho = out_span(rho);
  outputs.history_tau_flat = out_span(history_tau_flat);
  outputs.history_alpha = out_span(history_alpha);
  outputs.history_rel_delta = out_span(history_rel_delta);
  outputs.history_rel_delta_mean = out_span(history_rel_delta_mean);
  outputs.history_n_nan_genes = out_span(history_n_nan_genes);
  outputs.history_n_nonfinite = out_span(history_n_nonfinite);
  outputs.history_tau_max_seen = out_span(history_tau_max_seen);
  outputs.iterations_run = &iterations_run;
  outputs.converged = &converged;

  const pace::Status status = pace::run_irls_loop(inputs, options, r_log_nbinom, r_trigamma,
                                                  reporter, outputs, user_interrupted);
  raise_if_failed(status, "IRLS loop");

  // The warnings are RETURNED, not raised here. Rf_warningcall() under
  // options(warn = 2) -- or any calling handler promoting a warning to an
  // error -- longjmps out of this function, and a longjmp is not a C++ unwind,
  // so the destructors of count_holder, ambient_holder, z_holder and the design
  // buffers never run. Those holders release their SEXPs from R's precious list
  // in their destructors, so the counts, the ambient field and Z would stay
  // protected for the life of the session, along with the design vectors (some
  // 10 MB a block at 1.2M cells). Buffering in the core moved that longjmp one
  // frame out; raising from R removes it.
  return Rcpp::List::create(
      Rcpp::Named("warnings") = Rcpp::wrap(deferred_warnings),
      Rcpp::Named("B") = beta_out, Rcpp::Named("U") = u_out, Rcpp::Named("se_B") = se_beta,
      Rcpp::Named("se_U") = se_u, Rcpp::Named("alpha") = alpha,
      Rcpp::Named("tau_g_array") = tau_g_array, Rcpp::Named("tau_flat") = tau_flat,
      Rcpp::Named("rho") = rho, Rcpp::Named("history_tau_flat") = history_tau_flat,
      Rcpp::Named("history_alpha") = history_alpha,
      Rcpp::Named("history_rel_delta") = history_rel_delta,
      Rcpp::Named("history_rel_delta_mean") = history_rel_delta_mean,
      Rcpp::Named("history_n_nan_genes") = history_n_nan_genes,
      Rcpp::Named("history_n_nonfinite") = history_n_nonfinite,
      Rcpp::Named("history_tau_max_seen") = history_tau_max_seen,
      Rcpp::Named("n_iter") = static_cast<double>(iterations_run),
      Rcpp::Named("converged") = converged);
}
