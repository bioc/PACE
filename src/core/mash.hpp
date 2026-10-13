// mash.hpp -- canonical multivariate adaptive shrinkage, R-free.
//
// Reproduces the subset of mashr that PACE uses: canonical covariances, the
// autoselected grid, the mixture likelihood, EM mixture weights, and posterior
// means, standard deviations and lfsr. Data-driven covariances (cov_pca,
// cov_ed) are deliberately absent -- dropping them costs 11 of 1637 breast
// cancer calls and 0 of 46 melanoma calls, Jaccard 0.993 and 1.000.
//
// The mixture-weight optimiser is EM rather than mashr's default mixSQP. The
// two give byte-identical call sets on both canonical cohorts, Jaccard 1.0000,
// which is what makes this port tractable: mixSQP is a sequential quadratic
// program, EM is a fixed-point iteration.
#ifndef PACE_MASH_HPP
#define PACE_MASH_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "core_types.hpp"

namespace pace {

// One R x R covariance, stored column-major like every other matrix in the core.
struct CovarianceComponent {
  std::string name;
  std::vector<double> values;  // R * R, column-major
};

// The canonical set, in mashr's order: identity, then one singleton per
// condition, then equal effects, then simple heterogeneous at rho = 0.25, 0.5
// and 0.75. That is R + 5 components for R conditions.
//
// Each is normalised the way mashr normalises: divided by the largest diagonal
// entry, unless every diagonal entry is zero. For this set the largest diagonal
// is already 1, so normalisation is an identity, but it is applied so that a
// caller-supplied component would behave the same way.
//
// `condition_names` may be empty; when it has R entries the singletons take
// those names, as mashr names them after the columns of Bhat.
std::vector<CovarianceComponent> canonical_covariances(
    std::int64_t n_conditions, const std::vector<std::string>& condition_names);

// The grid of scaling factors, on the STANDARD DEVIATION scale: a component's
// covariance is grid[l]^2 * U_k.
//
//   gmin   = min(Shat) / 10
//   gmax   = 2 * sqrt(max(Bhat^2 - Shat^2)), or 8 * gmin when no effect exceeds
//            its own standard error
//   npoint = ceil(log2(gmax / gmin) / log2(mult))
//   grid   = mult^(-npoint .. 0) * gmax
//
// `mult` is mashr's gridmult and defaults to sqrt(2). Entries of Shat that are
// zero or non-finite, and of Bhat that are NaN, are excluded from both ends,
// exactly as mashr excludes them.
Status autoselect_grid(Span<const double> bhat, Span<const double> shat,
                       std::int64_t n_effects, std::int64_t n_conditions, double mult,
                       std::vector<double>* grid);

// Expanded mixture components, in mashr's order: the null (all zeros) first
// when `usepointmass` is true, then grid-major -- for each grid value in turn,
// every canonical covariance. So component (l, k) sits at 1 + l * K + k.
//
// The covariance of component (l, k) is grid[l]^2 * U_k. `grid` is on the
// standard deviation scale, which is why it enters squared.
std::int64_t n_mixture_components(std::int64_t n_covariances, std::int64_t n_grid,
                                  bool use_point_mass);

// Log density of every effect under every mixture component.
//
//   loglik[j, c] = log N(bhat_j ; 0, Sigma_c + V_j)
//   V_j          = diag(shat_j) * V * diag(shat_j)
//
// with Sigma_c the expanded component above and V the shared correlation of the
// measurement error. Pass an empty `v` for V = I, which is what PACE uses
// (`null_correlation = FALSE`, as published).
//
// `loglik` is filled column-major, n_effects x n_mixture_components. An effect
// whose covariance is not positive definite under some component yields
// -infinity for that component rather than failing the whole call, matching
// mashr's treatment of a degenerate slice.
Status mixture_loglik(Span<const double> bhat, Span<const double> shat,
                      Span<const double> v, std::int64_t n_effects, std::int64_t n_conditions,
                      const std::vector<CovarianceComponent>& covariances,
                      const std::vector<double>& grid, bool use_point_mass, int n_threads,
                      Span<double> loglik);

// EM for the mixture weights, on the penalised likelihood ashr maximises:
//
//   pi <- normalise( colSums(classprob) + prior - 1 )
//   classprob[j, c] = pi_c L_jc / sum_c' pi_c' L_jc'
//
// `L` is formed from `loglik` by subtracting each effect's maximum and
// exponentiating, which leaves classprob unchanged and keeps the exponential in
// range. The prior is mashr's "nullbiased": `null_weight` on the first
// component, 1 on the rest, so the default of 10 pulls weight toward the null.
//
// ashr accelerates this fixpoint with SQUAREM; plain iteration reaches the same
// fixed point, which is what the weights are defined by.
Status optimize_mixture_weights(Span<const double> loglik, std::int64_t n_effects,
                                std::int64_t n_components, double null_weight,
                                std::int64_t max_iterations, double tolerance,
                                std::vector<double>* weights, std::int64_t* iterations_used);

// Posterior summaries, mixed over the fitted components.
//
// Per effect j and component c, with V_j the measurement covariance:
//
//   U1   = U_c (V_j^-1 U_c + I)^-1          (posterior covariance)
//   mu1  = U1 V_j^-1 bhat_j                 (posterior mean)
//   var  = max(0, diag(U1))
//   neg  = P(beta < 0) = Phi(-mu1 / sqrt(var)), and 0 when var is 0
//   zero = 1 when var is 0, else 0
//
// mixed with the posterior weights w_jc proportional to pi_c L_jc, then
//
//   sd   = sqrt( E[beta^2] - E[beta]^2 )
//   lfsr = neg > 0.5 (1 - zero) ? 1 - neg : neg + zero,  floored at 0
//
// `mean`, `sd` and `lfsr` are each filled column-major, n_effects x
// n_conditions. This assumes mashr's default alpha = 0, which is what PACE
// uses, so the Shat_alpha rescaling is the identity.
Status posterior_summaries(Span<const double> bhat, Span<const double> shat,
                           Span<const double> v, std::int64_t n_effects,
                           std::int64_t n_conditions,
                           const std::vector<CovarianceComponent>& covariances,
                           const std::vector<double>& grid, bool use_point_mass,
                           const std::vector<double>& weights, Span<const double> loglik,
                           Span<double> mean, Span<double> sd, Span<double> lfsr);

}  // namespace pace

#endif  // PACE_MASH_HPP
