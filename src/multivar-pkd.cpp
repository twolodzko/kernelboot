
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "kernels.h"
#include "shared.h"


// [[Rcpp::export]]
Rcpp::List cpp_dmvpk(
    const arma::mat& x,
    const arma::mat& y,
    const arma::vec& bandwidth,
    const arma::vec& weights,
    const std::string& kernel = "gaussian",
    const bool& log_prob = false
  ) {

  double (*dens_kern)(double, double);

  if (kernel == "rectangular") {
    dens_kern = dens_rect;
  } else if (kernel == "triangular") {
    dens_kern = dens_triang;
  } else if (kernel == "biweight") {
    dens_kern = dens_biweight;
  } else if (kernel == "triweight") {
    dens_kern = dens_triweight;
  } else if (kernel == "cosine") {
    dens_kern = dens_cosine;
  } else if (kernel == "optcosine") {
    dens_kern = dens_optcos;
  } else if (kernel == "epanechnikov") {
    dens_kern = dens_epan;
  } else {
    dens_kern = dens_gauss;
  }

  const unsigned int n = x.n_rows;
  const unsigned int m = y.n_cols;
  const unsigned int k = y.n_rows;
  arma::vec p(n), c_weights(k);

  if (x.n_cols != m || bandwidth.n_elem != m)
    Rcpp::stop("dimmensions of x, y and bandwidth do not match");

  if (arma::any(bandwidth <= 0.0))
    Rcpp::stop("bandwidth needs to be positive");

  if (arma::any(weights < 0.0))
    Rcpp::stop("weights need to be non-negative");

  try {

    if (weights.n_elem == 1) {
      c_weights.fill( 1.0/static_cast<double>(k) );
    } else {
      if (weights.n_elem != k)
        Rcpp::stop("dimmensions of weights and y do not match");
      c_weights = weights;
    }

    c_weights /= arma::sum(c_weights);

    double tmp;
    for (unsigned int i = 0; i < n; i++) {
      p[i] = 0.0;
      for (unsigned int j = 0; j < k; j++) {
        tmp = 0.0;
        for (unsigned int l = 0; l < m; l++)
          tmp += std::log(dens_kern(x(i, l) - y(j, l), bandwidth[l]));
        p[i] += std::exp(tmp) * c_weights[j];
      }
    }

    if (log_prob)
      p = arma::log(p);

  } catch ( std::exception& __ex__ ) {
    forward_exception_to_r(__ex__);
  } catch (...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return Rcpp::List::create(
    Rcpp::Named("density") = p,
    Rcpp::Named("data") = y,
    Rcpp::Named("bandwidth") = bandwidth,
    Rcpp::Named("weights") = c_weights,
    Rcpp::Named("kernel") = kernel,
    Rcpp::Named("log_prob") = log_prob
  );

}


// [[Rcpp::export]]
Rcpp::List cpp_rmvpk(
    const unsigned int& n,
    const arma::mat& y,
    const arma::vec& bandwidth,
    const arma::vec& weights,
    const std::string& kernel = "gaussian"
  ) {

  double (*rng_kern)();

  if (kernel == "rectangular") {
    rng_kern = rng_rect;
  } else if (kernel == "triangular") {
    rng_kern = rng_triang;
  } else if (kernel == "biweight") {
    rng_kern = rng_biweight;
  } else if (kernel == "triweight") {
    rng_kern = rng_triweight;
  } else if (kernel == "cosine") {
    rng_kern = rng_cosine;
  } else if (kernel == "optcosine") {
    rng_kern = rng_optcos;
  } else if (kernel == "epanechnikov") {
    rng_kern = rng_epan;
  } else {
    rng_kern = R::norm_rand;
  }

  const unsigned int m = y.n_cols;
  const unsigned int k = y.n_rows;
  arma::mat samp(n, m);
  arma::vec c_weights(k);
  std::vector<unsigned int> idx(n);

  if (bandwidth.n_elem != m)
    Rcpp::stop("dimmensions of y and bandwidth do not match");

  if (arma::any(bandwidth <= 0.0))
    Rcpp::stop("bandwidth needs to be positive");

  if (arma::any(weights < 0.0))
    Rcpp::stop("weights need to be non-negative");

  try {

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

    unsigned int j;
    for (unsigned int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      for (unsigned int l = 0; l < m; l++)
        samp(i, l) = y(j, l) + rng_kern() * bandwidth[l];
    }

    for (unsigned int i = (k-1); i > 0; i--)
      c_weights[i] -= c_weights[i-1];

  } catch ( std::exception& __ex__ ) {
    forward_exception_to_r(__ex__);
  } catch (...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return Rcpp::List::create(
    Rcpp::Named("samples") = samp,
    Rcpp::Named("boot_index") = idx,
    Rcpp::Named("data") = y,
    Rcpp::Named("bandwidth") = bandwidth,
    Rcpp::Named("weights") = c_weights,
    Rcpp::Named("kernel") = kernel
  );

}

