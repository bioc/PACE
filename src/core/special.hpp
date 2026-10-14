// special.hpp -- the special functions the core needs without R.
//
// The R bindings hand the solver R's own Rf_trigamma. A Python binding has no
// such thing, so the core carries its own, which must agree with R's closely
// enough that a fit does not depend on which language called it.
#ifndef PACE_SPECIAL_HPP
#define PACE_SPECIAL_HPP

namespace pace {

// The trigamma function, psi'(x), for finite x.
//
//   x > 0  : recurrence up to the asymptotic regime, then the Bernoulli series
//   x <= 0 : the reflection formula, except at the poles (non-positive
//            integers) where the value is infinite
//
// Agrees with R's Rf_trigamma to better than 1e-15 relative over the range
// estimate_d0() searches; see tests/testthat/test-special.R.
double trigamma(double x);

}  // namespace pace

#endif  // PACE_SPECIAL_HPP
