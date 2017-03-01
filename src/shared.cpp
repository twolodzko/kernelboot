

#include <RcppArmadillo.h>
#include "shared.h"


double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}

unsigned int sample_int(arma::vec cumul_weights) {
  unsigned int j;
  double u = rng_unif();
  for (j = 0; j < cumul_weights.n_elem; j++) {
    if (cumul_weights[j] >= u)
      break;
  }
  return j;
}
