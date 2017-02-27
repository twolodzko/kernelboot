
#ifndef MVN_H
#define MVN_H

#include <RcppArmadillo.h>

arma::vec cpp_dmvn(arma::mat x, arma::rowvec mu, arma::mat sigma,
                   bool log_prob = false);

arma::mat cpp_rmvn(int n, arma::vec mu, arma::mat sigma);

#endif
