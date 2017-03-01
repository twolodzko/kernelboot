
#ifndef KB_SHARED_H
#define KB_SHARED_H

#include <RcppArmadillo.h>

double rng_unif();

unsigned int sample_int(const arma::vec& cumul_weights);

// inline bool is_diag(const arma::mat& x, double tol = 1e-12) {
//   return (x.n_cols == x.n_rows) && arma::approx_equal(arma::diagmat(x), x, "absdiff", tol);
// }

#endif
