
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "kernels.h"
#include "shared.h"


// [[Rcpp::export]]
Rcpp::List cpp_ruvk(
    const unsigned int& n,
    const arma::vec& y,
    const double& bandwidth,
    const arma::vec& weights,
    const std::string& kernel = "gaussian",
    const bool& shrinked = false
  ) {

  double (*rng_kern)();

  if (kernel == "rectangular") {
    rng_kern = rng_rect;
  } else if (kernel == "triangular") {
    rng_kern = rng_triang;
  } else if (kernel == "biweight") {
    rng_kern = rng_biweight;
  } else if (kernel == "cosine") {
    rng_kern = rng_cosine;
  } else if (kernel == "optcosine") {
    rng_kern = rng_optcos;
  } else if (kernel == "epanechnikov") {
    rng_kern = rng_epan;
  } else {
    rng_kern = R::norm_rand;
  }

  const unsigned int k = y.n_elem;
  arma::vec samp(n);
  arma::vec c_weights(k);
  std::vector<unsigned int> idx(n);

  if (bandwidth <= 0.0)
    Rcpp::stop("bandwidth needs to be positive");

  if (!R_FINITE(bandwidth))
    Rcpp::stop("inappropriate value of bandwidth");

  if (arma::any(weights < 0.0))
    Rcpp::stop("weights need to be non-negative");

  if (!arma::is_finite(weights))
    Rcpp::stop("inappropriate values of weights");

  if (weights.n_elem == 1) {
    c_weights.fill( 1.0/static_cast<double>(k) );
  } else {
    if (weights.n_elem != k)
      Rcpp::stop("dimmensions of weights and y do not match");
    c_weights = weights;
  }

  for (unsigned int i = 1; i < k; i++)
    c_weights[i] += c_weights[i-1];
  c_weights /= c_weights[k-1];

  if (k == 1) {

    for (unsigned int i = 0; i < n; i++) {
      samp[i] = y[0] + rng_kern() * bandwidth;
    }

  } else if (!shrinked) {

    unsigned int j;
    for (unsigned int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      samp[i] = y[j] + rng_kern() * bandwidth;
    }

  } else {

    const double my = arma::mean(y);
    const double sy = arma::var(y);
    const double c = std::sqrt(1.0 + (bandwidth*bandwidth)/sy);

    unsigned int j;
    for (unsigned int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      samp[i] = my + (y[j] - my + rng_kern() * bandwidth) / c;
    }

  }

  for (unsigned int i = (k-1); i > 0; i--)
    c_weights[i] -= c_weights[i-1];

  return Rcpp::List::create(
    Rcpp::Named("samples") = samp,
    Rcpp::Named("boot_index") = idx,
    Rcpp::Named("data") = y,
    Rcpp::Named("bandwidth") = bandwidth,
    Rcpp::Named("weights") = c_weights,
    Rcpp::Named("kernel") = kernel,
    Rcpp::Named("shrinked") = shrinked
  );

}

