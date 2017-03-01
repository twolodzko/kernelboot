
#ifndef SHARED_H
#define SHARED_H

#include <RcppArmadillo.h>

double rng_unif();

unsigned int sample_int(arma::vec cumul_weights);

#endif
