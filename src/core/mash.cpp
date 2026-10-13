#include "mash.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pace {
namespace {

// mashr's normalize_cov: divide by the largest diagonal entry, unless every
// diagonal entry is zero (the null component, which must stay all zeros).
void normalise_in_place(std::vector<double>* values, std::int64_t n) {
  double largest_diagonal = 0.0;
  for (std::int64_t i = 0; i < n; ++i) largest_diagonal = std::max(largest_diagonal, (*values)[i * n + i]);
  if (largest_diagonal == 0.0) return;
  for (double& value : *values) value /= largest_diagonal;
}

CovarianceComponent make_component(std::string name, std::int64_t n) {
  CovarianceComponent component;
  component.name = std::move(name);
  component.values.assign(static_cast<std::size_t>(n * n), 0.0);
  return component;
}

}  // namespace

std::vector<CovarianceComponent> canonical_covariances(
    std::int64_t n_conditions, const std::vector<std::string>& condition_names) {
  const std::int64_t n = n_conditions;
  std::vector<CovarianceComponent> out;
  out.reserve(static_cast<std::size_t>(n + 5));

  // identity
  {
    CovarianceComponent identity = make_component("identity", n);
    for (std::int64_t i = 0; i < n; ++i) identity.values[i * n + i] = 1.0;
    out.push_back(std::move(identity));
  }

  // one singleton per condition: a single 1 on the diagonal
  const bool have_names = static_cast<std::int64_t>(condition_names.size()) == n;
  for (std::int64_t r = 0; r < n; ++r) {
    CovarianceComponent singleton =
        make_component(have_names ? condition_names[static_cast<std::size_t>(r)]
                                  : "singleton_" + std::to_string(r + 1), n);
    singleton.values[r * n + r] = 1.0;
    out.push_back(std::move(singleton));
  }

  // equal effects: every entry 1
  {
    CovarianceComponent equal = make_component("equal_effects", n);
    std::fill(equal.values.begin(), equal.values.end(), 1.0);
    out.push_back(std::move(equal));
  }

  // simple heterogeneous: off-diagonal rho, diagonal 1
  const double correlations[3] = {0.25, 0.5, 0.75};
  for (int k = 0; k < 3; ++k) {
    CovarianceComponent het = make_component("simple_het_" + std::to_string(k + 1), n);
    std::fill(het.values.begin(), het.values.end(), correlations[k]);
    for (std::int64_t i = 0; i < n; ++i) het.values[i * n + i] = 1.0;
    out.push_back(std::move(het));
  }

  for (CovarianceComponent& component : out) normalise_in_place(&component.values, n);
  return out;
}

Status autoselect_grid(Span<const double> bhat, Span<const double> shat,
                       std::int64_t n_effects, std::int64_t n_conditions, double mult,
                       std::vector<double>* grid) {
  if (grid == nullptr) return Status::failure(StatusCode::invalid_argument, "grid must not be null");
  const std::int64_t n_values = n_effects * n_conditions;
  if (bhat.size != n_values || shat.size != n_values) {
    return Status::failure(StatusCode::invalid_argument,
                           "bhat and shat must each hold n_effects * n_conditions values");
  }
  if (!(mult > 1.0)) {
    return Status::failure(StatusCode::invalid_argument, "gridmult must be greater than 1");
  }

  // mashr excludes a value when its standard error is zero or non-finite, or
  // its effect is NA, and takes both ends over what remains.
  double smallest_se = std::numeric_limits<double>::infinity();
  double largest_excess = -std::numeric_limits<double>::infinity();
  std::int64_t n_used = 0;
  for (std::int64_t i = 0; i < n_values; ++i) {
    const double se = shat[i];
    const double effect = bhat[i];
    if (se == 0.0 || !std::isfinite(se) || std::isnan(effect)) continue;
    ++n_used;
    smallest_se = std::min(smallest_se, se);
    largest_excess = std::max(largest_excess, effect * effect - se * se);
  }
  if (n_used == 0) {
    return Status::failure(StatusCode::invalid_argument, "every effect was excluded from the grid");
  }

  const double grid_min = smallest_se / 10.0;
  const double grid_max = largest_excess <= 0.0 ? 8.0 * grid_min : 2.0 * std::sqrt(largest_excess);

  const std::int64_t n_point =
      static_cast<std::int64_t>(std::ceil(std::log2(grid_max / grid_min) / std::log2(mult)));
  grid->clear();
  grid->reserve(static_cast<std::size_t>(n_point + 1));
  for (std::int64_t k = -n_point; k <= 0; ++k) {
    grid->push_back(std::pow(mult, static_cast<double>(k)) * grid_max);
  }
  return Status::success();
}


std::int64_t n_mixture_components(std::int64_t n_covariances, std::int64_t n_grid,
                                  bool use_point_mass) {
  return (use_point_mass ? 1 : 0) + n_covariances * n_grid;
}

namespace {

constexpr double kLog2Pi = 1.8378770664093454835606594728112;

// In-place lower Cholesky of a symmetric n x n matrix held column-major.
// Returns false when the matrix is not positive definite, leaving `a` partly
// overwritten; the caller discards it in that case.
bool cholesky_lower(double* a, std::int64_t n) {
  for (std::int64_t j = 0; j < n; ++j) {
    double diagonal = a[j * n + j];
    for (std::int64_t k = 0; k < j; ++k) {
      const double value = a[k * n + j];
      diagonal -= value * value;
    }
    if (!(diagonal > 0.0)) return false;
    const double root = std::sqrt(diagonal);
    a[j * n + j] = root;
    for (std::int64_t i = j + 1; i < n; ++i) {
      double sum = a[j * n + i];
      for (std::int64_t k = 0; k < j; ++k) sum -= a[k * n + i] * a[k * n + j];
      a[j * n + i] = sum / root;
    }
  }
  return true;
}

// Log N(b ; 0, S) with S already decomposed into its lower Cholesky factor.
// scratch must hold n doubles.
double log_density_from_chol(const double* chol, const double* b, std::int64_t n, double* scratch) {
  double log_determinant = 0.0;
  for (std::int64_t i = 0; i < n; ++i) log_determinant += std::log(chol[i * n + i]);
  // forward-solve L z = b, then the quadratic form is z'z
  for (std::int64_t i = 0; i < n; ++i) {
    double sum = b[i];
    for (std::int64_t k = 0; k < i; ++k) sum -= chol[k * n + i] * scratch[k];
    scratch[i] = sum / chol[i * n + i];
  }
  double quadratic = 0.0;
  for (std::int64_t i = 0; i < n; ++i) quadratic += scratch[i] * scratch[i];
  return -0.5 * (static_cast<double>(n) * kLog2Pi + quadratic) - log_determinant;
}

}  // namespace

Status mixture_loglik(Span<const double> bhat, Span<const double> shat,
                      Span<const double> v, std::int64_t n_effects, std::int64_t n_conditions,
                      const std::vector<CovarianceComponent>& covariances,
                      const std::vector<double>& grid, bool use_point_mass, int n_threads,
                      Span<double> loglik) {
  (void)n_threads;  // the per-effect loop is trivially parallel; serial for now
  const std::int64_t n = n_conditions;
  const std::int64_t n_values = n_effects * n;
  if (bhat.size != n_values || shat.size != n_values) {
    return Status::failure(StatusCode::invalid_argument,
                           "bhat and shat must each hold n_effects * n_conditions values");
  }
  if (v.size != 0 && v.size != n * n) {
    return Status::failure(StatusCode::invalid_argument, "v must be empty or n_conditions^2");
  }
  const std::int64_t n_cov = static_cast<std::int64_t>(covariances.size());
  const std::int64_t n_grid = static_cast<std::int64_t>(grid.size());
  const std::int64_t n_comp = n_mixture_components(n_cov, n_grid, use_point_mass);
  if (loglik.size != n_effects * n_comp) {
    return Status::failure(StatusCode::invalid_argument,
                           "loglik must hold n_effects * n_components values");
  }
  for (const CovarianceComponent& component : covariances) {
    if (static_cast<std::int64_t>(component.values.size()) != n * n) {
      return Status::failure(StatusCode::invalid_argument,
                             "every covariance must be n_conditions x n_conditions");
    }
  }

  std::vector<double> error_cov(static_cast<std::size_t>(n * n));
  std::vector<double> sigma(static_cast<std::size_t>(n * n));
  std::vector<double> effect(static_cast<std::size_t>(n));
  std::vector<double> scratch(static_cast<std::size_t>(n));

  for (std::int64_t j = 0; j < n_effects; ++j) {
    for (std::int64_t i = 0; i < n; ++i) effect[i] = bhat[i * n_effects + j];

    // V_j = diag(shat_j) * V * diag(shat_j); with V empty this is diag(shat_j^2)
    for (std::int64_t c = 0; c < n; ++c) {
      const double sc = shat[c * n_effects + j];
      for (std::int64_t r = 0; r < n; ++r) {
        const double sr = shat[r * n_effects + j];
        const double v_rc = v.size == 0 ? (r == c ? 1.0 : 0.0) : v[c * n + r];
        error_cov[c * n + r] = sr * v_rc * sc;
      }
    }

    if (use_point_mass) {
      sigma = error_cov;
      loglik[0 * n_effects + j] = cholesky_lower(sigma.data(), n)
                                      ? log_density_from_chol(sigma.data(), effect.data(), n,
                                                              scratch.data())
                                      : -std::numeric_limits<double>::infinity();
    }
    const std::int64_t offset = use_point_mass ? 1 : 0;
    for (std::int64_t l = 0; l < n_grid; ++l) {
      const double scale = grid[static_cast<std::size_t>(l)] * grid[static_cast<std::size_t>(l)];
      for (std::int64_t k = 0; k < n_cov; ++k) {
        const std::vector<double>& u = covariances[static_cast<std::size_t>(k)].values;
        for (std::int64_t i = 0; i < n * n; ++i) sigma[static_cast<std::size_t>(i)] =
            scale * u[static_cast<std::size_t>(i)] + error_cov[static_cast<std::size_t>(i)];
        const std::int64_t component = offset + l * n_cov + k;
        loglik[component * n_effects + j] =
            cholesky_lower(sigma.data(), n)
                ? log_density_from_chol(sigma.data(), effect.data(), n, scratch.data())
                : -std::numeric_limits<double>::infinity();
      }
    }
  }
  return Status::success();
}


Status optimize_mixture_weights(Span<const double> loglik, std::int64_t n_effects,
                                std::int64_t n_components, double null_weight,
                                std::int64_t max_iterations, double tolerance,
                                std::vector<double>* weights, std::int64_t* iterations_used) {
  if (weights == nullptr) return Status::failure(StatusCode::invalid_argument, "weights must not be null");
  if (loglik.size != n_effects * n_components) {
    return Status::failure(StatusCode::invalid_argument,
                           "loglik must hold n_effects * n_components values");
  }
  const std::size_t n_cells = static_cast<std::size_t>(n_effects * n_components);

  // Relative likelihood: subtract each effect's maximum before exponentiating.
  // classprob is a ratio within an effect, so the subtraction cancels exactly.
  std::vector<double> lik(n_cells);
  for (std::int64_t j = 0; j < n_effects; ++j) {
    double row_max = -std::numeric_limits<double>::infinity();
    for (std::int64_t c = 0; c < n_components; ++c)
      row_max = std::max(row_max, loglik[c * n_effects + j]);
    if (!std::isfinite(row_max)) {
      return Status::failure(StatusCode::invalid_argument,
                             "an effect has no finite likelihood under any component");
    }
    for (std::int64_t c = 0; c < n_components; ++c)
      lik[static_cast<std::size_t>(c * n_effects + j)] =
          std::exp(loglik[c * n_effects + j] - row_max);
  }

  std::vector<double> prior(static_cast<std::size_t>(n_components), 1.0);
  if (n_components > 0) prior[0] = null_weight;

  weights->assign(static_cast<std::size_t>(n_components), 1.0 / static_cast<double>(n_components));
  std::vector<double> next(static_cast<std::size_t>(n_components));

  std::int64_t iteration = 0;
  for (; iteration < max_iterations; ++iteration) {
    std::fill(next.begin(), next.end(), 0.0);
    for (std::int64_t j = 0; j < n_effects; ++j) {
      double row_sum = 0.0;
      for (std::int64_t c = 0; c < n_components; ++c)
        row_sum += (*weights)[static_cast<std::size_t>(c)] * lik[static_cast<std::size_t>(c * n_effects + j)];
      if (!(row_sum > 0.0)) continue;  // an effect that no component explains adds nothing
      for (std::int64_t c = 0; c < n_components; ++c)
        next[static_cast<std::size_t>(c)] +=
            (*weights)[static_cast<std::size_t>(c)] *
            lik[static_cast<std::size_t>(c * n_effects + j)] / row_sum;
    }
    double total = 0.0;
    for (std::int64_t c = 0; c < n_components; ++c) {
      next[static_cast<std::size_t>(c)] =
          std::max(0.0, next[static_cast<std::size_t>(c)] + prior[static_cast<std::size_t>(c)] - 1.0);
      total += next[static_cast<std::size_t>(c)];
    }
    if (!(total > 0.0)) return Status::failure(StatusCode::internal_error,
                                               "mixture weights collapsed to zero");
    double shift = 0.0;
    for (std::int64_t c = 0; c < n_components; ++c) {
      next[static_cast<std::size_t>(c)] /= total;
      shift = std::max(shift, std::fabs(next[static_cast<std::size_t>(c)] -
                                        (*weights)[static_cast<std::size_t>(c)]));
    }
    weights->swap(next);
    if (shift < tolerance) { ++iteration; break; }
  }
  if (iterations_used != nullptr) *iterations_used = iteration;
  return Status::success();
}


namespace {

// Solve A X = B in place for a general n x n A, via Gaussian elimination with
// partial pivoting. Returns false if A is singular. `a` and `b` are both
// column-major, b holding n right-hand sides.
bool solve_general(double* a, double* b, std::int64_t n) {
  for (std::int64_t col = 0; col < n; ++col) {
    std::int64_t pivot = col;
    double best = std::fabs(a[col * n + col]);
    for (std::int64_t r = col + 1; r < n; ++r) {
      const double candidate = std::fabs(a[col * n + r]);
      if (candidate > best) { best = candidate; pivot = r; }
    }
    if (!(best > 0.0)) return false;
    if (pivot != col) {
      for (std::int64_t c = 0; c < n; ++c) std::swap(a[c * n + col], a[c * n + pivot]);
      for (std::int64_t c = 0; c < n; ++c) std::swap(b[c * n + col], b[c * n + pivot]);
    }
    const double diagonal = a[col * n + col];
    for (std::int64_t r = col + 1; r < n; ++r) {
      const double factor = a[col * n + r] / diagonal;
      if (factor == 0.0) continue;
      for (std::int64_t c = col; c < n; ++c) a[c * n + r] -= factor * a[c * n + col];
      for (std::int64_t c = 0; c < n; ++c) b[c * n + r] -= factor * b[c * n + col];
    }
  }
  for (std::int64_t col = n - 1; col >= 0; --col) {
    const double diagonal = a[col * n + col];
    for (std::int64_t c = 0; c < n; ++c) {
      double sum = b[c * n + col];
      for (std::int64_t k = col + 1; k < n; ++k) sum -= a[k * n + col] * b[c * n + k];
      b[c * n + col] = sum / diagonal;
    }
  }
  return true;
}

// Standard normal CDF, as R's pnorm computes it for finite arguments.
// 1/sqrt(2) is spelled out rather than taken from M_SQRT1_2, which is POSIX
// rather than standard C++ and is absent on MinGW without _USE_MATH_DEFINES.
constexpr double kInvSqrt2 = 0.70710678118654752440084436210485;
double standard_normal_cdf(double x) { return 0.5 * std::erfc(-x * kInvSqrt2); }

}  // namespace

Status posterior_summaries(Span<const double> bhat, Span<const double> shat,
                           Span<const double> v, std::int64_t n_effects,
                           std::int64_t n_conditions,
                           const std::vector<CovarianceComponent>& covariances,
                           const std::vector<double>& grid, bool use_point_mass,
                           const std::vector<double>& weights, Span<const double> loglik,
                           Span<double> mean, Span<double> sd, Span<double> lfsr) {
  const std::int64_t n = n_conditions;
  const std::int64_t n_cov = static_cast<std::int64_t>(covariances.size());
  const std::int64_t n_grid = static_cast<std::int64_t>(grid.size());
  const std::int64_t n_comp = n_mixture_components(n_cov, n_grid, use_point_mass);
  if (static_cast<std::int64_t>(weights.size()) != n_comp) {
    return Status::failure(StatusCode::invalid_argument, "weights must hold one value per component");
  }
  if (loglik.size != n_effects * n_comp) {
    return Status::failure(StatusCode::invalid_argument, "loglik must be n_effects x n_components");
  }
  const std::int64_t n_out = n_effects * n;
  if (mean.size != n_out || sd.size != n_out || lfsr.size != n_out) {
    return Status::failure(StatusCode::invalid_argument,
                           "mean, sd and lfsr must each be n_effects x n_conditions");
  }

  const std::size_t nn = static_cast<std::size_t>(n * n);
  std::vector<double> error_cov(nn), v_inverse(nn), lhs(nn), u1(nn), scaled_u(nn);
  std::vector<double> effect(static_cast<std::size_t>(n)), v_inv_b(static_cast<std::size_t>(n));
  std::vector<double> post_weight(static_cast<std::size_t>(n_comp));
  std::vector<double> acc_mean(static_cast<std::size_t>(n)), acc_mean2(static_cast<std::size_t>(n));
  std::vector<double> acc_neg(static_cast<std::size_t>(n)), acc_zero(static_cast<std::size_t>(n));

  for (std::int64_t j = 0; j < n_effects; ++j) {
    for (std::int64_t i = 0; i < n; ++i) effect[static_cast<std::size_t>(i)] = bhat[i * n_effects + j];
    for (std::int64_t c = 0; c < n; ++c) {
      const double sc = shat[c * n_effects + j];
      for (std::int64_t r = 0; r < n; ++r) {
        const double sr = shat[r * n_effects + j];
        const double v_rc = v.size == 0 ? (r == c ? 1.0 : 0.0) : v[c * n + r];
        error_cov[static_cast<std::size_t>(c * n + r)] = sr * v_rc * sc;
      }
    }
    // V^-1 by solving V X = I
    lhs = error_cov;
    std::fill(v_inverse.begin(), v_inverse.end(), 0.0);
    for (std::int64_t i = 0; i < n; ++i) v_inverse[static_cast<std::size_t>(i * n + i)] = 1.0;
    if (!solve_general(lhs.data(), v_inverse.data(), n)) {
      return Status::failure(StatusCode::internal_error,
                             "measurement covariance is singular for an effect");
    }
    for (std::int64_t i = 0; i < n; ++i) {
      double sum = 0.0;
      for (std::int64_t k = 0; k < n; ++k)
        sum += v_inverse[static_cast<std::size_t>(k * n + i)] * effect[static_cast<std::size_t>(k)];
      v_inv_b[static_cast<std::size_t>(i)] = sum;
    }

    // posterior weights: pi_c L_jc normalised over c, using the relative
    // likelihood so the exponential stays in range
    double row_max = -std::numeric_limits<double>::infinity();
    for (std::int64_t c = 0; c < n_comp; ++c) row_max = std::max(row_max, loglik[c * n_effects + j]);
    double weight_total = 0.0;
    for (std::int64_t c = 0; c < n_comp; ++c) {
      const double w = weights[static_cast<std::size_t>(c)] *
                       std::exp(loglik[c * n_effects + j] - row_max);
      post_weight[static_cast<std::size_t>(c)] = w;
      weight_total += w;
    }
    if (!(weight_total > 0.0)) {
      return Status::failure(StatusCode::internal_error, "posterior weights vanished for an effect");
    }
    for (double& w : post_weight) w /= weight_total;

    std::fill(acc_mean.begin(), acc_mean.end(), 0.0);
    std::fill(acc_mean2.begin(), acc_mean2.end(), 0.0);
    std::fill(acc_neg.begin(), acc_neg.end(), 0.0);
    std::fill(acc_zero.begin(), acc_zero.end(), 0.0);

    for (std::int64_t c = 0; c < n_comp; ++c) {
      const double w = post_weight[static_cast<std::size_t>(c)];
      if (w == 0.0) continue;
      // the component's covariance: the null is all zeros, the rest grid^2 * U
      if (use_point_mass && c == 0) {
        std::fill(scaled_u.begin(), scaled_u.end(), 0.0);
      } else {
        const std::int64_t index = use_point_mass ? c - 1 : c;
        const std::int64_t l = index / n_cov, k = index % n_cov;
        const double scale = grid[static_cast<std::size_t>(l)] * grid[static_cast<std::size_t>(l)];
        const std::vector<double>& u = covariances[static_cast<std::size_t>(k)].values;
        for (std::size_t i = 0; i < nn; ++i) scaled_u[i] = scale * u[i];
      }
      // U1 = U (V^-1 U + I)^-1, obtained by solving (V^-1 U + I)' X' = U'
      // Here we form M = V^-1 U + I and solve M^T Y = U^T, so U1 = Y^T.
      std::fill(lhs.begin(), lhs.end(), 0.0);
      for (std::int64_t col = 0; col < n; ++col)
        for (std::int64_t row = 0; row < n; ++row) {
          double sum = 0.0;
          for (std::int64_t k2 = 0; k2 < n; ++k2)
            sum += v_inverse[static_cast<std::size_t>(k2 * n + row)] *
                   scaled_u[static_cast<std::size_t>(col * n + k2)];
          lhs[static_cast<std::size_t>(row * n + col)] = sum + (row == col ? 1.0 : 0.0);
        }
      for (std::int64_t col = 0; col < n; ++col)
        for (std::int64_t row = 0; row < n; ++row)
          u1[static_cast<std::size_t>(row * n + col)] = scaled_u[static_cast<std::size_t>(col * n + row)];
      if (!solve_general(lhs.data(), u1.data(), n)) continue;  // degenerate component contributes nothing

      for (std::int64_t i = 0; i < n; ++i) {
        double mu = 0.0;
        for (std::int64_t k2 = 0; k2 < n; ++k2)
          mu += u1[static_cast<std::size_t>(i * n + k2)] * v_inv_b[static_cast<std::size_t>(k2)];
        const double variance = std::max(0.0, u1[static_cast<std::size_t>(i * n + i)]);
        acc_mean[static_cast<std::size_t>(i)] += w * mu;
        acc_mean2[static_cast<std::size_t>(i)] += w * (mu * mu + variance);
        if (variance == 0.0) {
          acc_zero[static_cast<std::size_t>(i)] += w;
        } else {
          acc_neg[static_cast<std::size_t>(i)] += w * standard_normal_cdf(-mu / std::sqrt(variance));
        }
      }
    }

    for (std::int64_t i = 0; i < n; ++i) {
      const double m = acc_mean[static_cast<std::size_t>(i)];
      const double m2 = acc_mean2[static_cast<std::size_t>(i)];
      const double neg = acc_neg[static_cast<std::size_t>(i)];
      const double zero = acc_zero[static_cast<std::size_t>(i)];
      mean[i * n_effects + j] = m;
      sd[i * n_effects + j] = std::sqrt(std::max(0.0, m2 - m * m));
      const double value = neg > 0.5 * (1.0 - zero) ? 1.0 - neg : neg + zero;
      lfsr[i * n_effects + j] = value < 0.0 ? 0.0 : value;
    }
  }
  return Status::success();
}

}  // namespace pace
