
#ifndef KB_SHARED_H
#define KB_SHARED_H

#include <RcppArmadillo.h>

// sample from standard uniform distribution

inline double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}

// sample integers from discrete distribution
// cumul_weights is a vector of cumulative weights

inline unsigned int sample_int(const arma::vec& cumul_weights) {
  double u = rng_unif();
  unsigned int j;
  for (j = 0; j < cumul_weights.n_elem; j++) {
    if (cumul_weights[j] >= u)
      break;
  }
  return j;
}

// inline bool is_diag(const arma::mat& x, double tol = 1e-12) {
//   return (x.n_cols == x.n_rows) && arma::approx_equal(arma::diagmat(x), x, "absdiff", tol);
// }

#endif
