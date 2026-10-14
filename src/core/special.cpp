#include "special.hpp"

#include <cmath>
#include <limits>

namespace pace {
namespace {

// psi'(x) ~ 1/x + 1/(2x^2) + sum_k B_2k / x^(2k+1), accurate once x is large
// enough that the series is well past its smallest term.
constexpr double kAsymptoticThreshold = 20.0;

double trigamma_asymptotic(double x) {
  const double inv = 1.0 / x;
  const double inv2 = inv * inv;
  // Bernoulli coefficients B_2k / (2k)! * (2k)! = B_2k, to x^-15
  double series = 1.0;                       // leading 1/x term, factored below
  series += 0.5 * inv;                       // + 1/(2x^2)
  double power = inv2;                       // x^-2
  series += (1.0 / 6.0) * power;             // + 1/(6x^3)
  power *= inv2;
  series += (-1.0 / 30.0) * power;           // - 1/(30x^5)
  power *= inv2;
  series += (1.0 / 42.0) * power;            // + 1/(42x^7)
  power *= inv2;
  series += (-1.0 / 30.0) * power;           // - 1/(30x^9)
  power *= inv2;
  series += (5.0 / 66.0) * power;            // + 5/(66x^11)
  power *= inv2;
  series += (-691.0 / 2730.0) * power;       // - 691/(2730x^13)
  power *= inv2;
  series += (7.0 / 6.0) * power;             // + 7/(6x^15)
  return inv * series;
}

constexpr double kPi = 3.14159265358979323846264338327950288;

}  // namespace

double trigamma(double x) {
  if (std::isnan(x)) return x;

  if (x <= 0.0) {
    // A non-positive integer is a pole.
    if (x == std::floor(x)) return std::numeric_limits<double>::infinity();
    // Reflection: psi'(x) + psi'(1 - x) = pi^2 / sin^2(pi x)
    const double s = std::sin(kPi * x);
    return (kPi * kPi) / (s * s) - trigamma(1.0 - x);
  }

  // Recurrence psi'(x) = psi'(x + 1) + 1/x^2, applied until the asymptotic
  // series is accurate, then unwound.
  double result = 0.0;
  double value = x;
  while (value < kAsymptoticThreshold) {
    result += 1.0 / (value * value);
    value += 1.0;
  }
  return result + trigamma_asymptotic(value);
}

}  // namespace pace
