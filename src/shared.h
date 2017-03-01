
#ifndef KB_SHARED_H
#define KB_SHARED_H

#include <RcppArmadillo.h>

double rng_unif();
unsigned int sample_int(const arma::vec& cumul_weights);

#endif
