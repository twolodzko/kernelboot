
#include <RcppArmadillo.h>
#include "kernels.h"
#include "shared.h"


// [[Rcpp::export]]
Rcpp::List cpp_duvkd(
    const arma::vec& x,
    const arma::vec& y,
    const double& bandwidth,
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

  const unsigned int k = y.n_elem;
  const unsigned int n = x.n_elem;
  arma::vec p(n), c_weights(k);

  if (bandwidth <= 0.0)
    Rcpp::stop("bandwidth needs to be greater than zero");

  if (any(weights < 0.0))
    Rcpp::stop("weights need to be non-negative");

  if (weights.n_elem == 1) {
    c_weights.fill( 1.0/static_cast<double>(k) );
  } else {
    if (weights.n_elem != k)
      Rcpp::stop("length(weights) != length(y)");
    c_weights = weights;
  }

  c_weights /= sum(c_weights);

  for (unsigned int i = 0; i < n; i++) {
    if (ISNAN(y[i])) {
      p[i] = y[i];
      continue;
    }
    p[i] = 0.0;
    for (unsigned int j = 0; j < k; j++)
      p[i] += dens_kern(x[i] - y[j], bandwidth) * c_weights[j];
  }

  if (log_prob)
    p = log(p);

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
Rcpp::List cpp_ruvkd(
    const unsigned int& n,
    const arma::vec& y,
    const double& bandwidth,
    const arma::vec& weights,
    const std::string& kernel = "gaussian",
    const bool& preserve_var = false
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

  const unsigned int k = y.n_elem;
  arma::vec samp(n);
  arma::vec c_weights(k);
  std::vector<unsigned int> idx(n);

  if (bandwidth < 0.0)
    Rcpp::stop("bandwidth needs to be non-negative");

  if (any(weights < 0.0))
    Rcpp::stop("weights need to be non-negative");

  if (weights.n_elem == 1) {
    c_weights.fill( 1.0/static_cast<double>(k) );
  } else {
    if (weights.n_elem != k)
      Rcpp::stop("length(weights) != length(y)");
    c_weights = weights;
  }

  for (unsigned int i = 1; i < k; i++)
    c_weights[i] += c_weights[i-1];
  c_weights /= c_weights[k-1];

  if (!preserve_var) {

    unsigned int j;
    for (unsigned int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      samp[i] = y[j] + rng_kern() * bandwidth;
    }

  } else {

    double my, sy, c;
    my = mean(y);
    sy = var(y);
    c = sqrt(1.0 + pow(bandwidth, 2.0)/sy);

    unsigned int j;
    for (unsigned int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      samp[i] = my + (y[j] - my + rng_kern() * bandwidth) / c;
    }

  }

  for (unsigned int i = k; i > 0; i--)
    c_weights[i] -= c_weights[i-1];

  return Rcpp::List::create(
    Rcpp::Named("samples") = samp,
    Rcpp::Named("boot_index") = idx,
    Rcpp::Named("data") = y,
    Rcpp::Named("bandwidth") = bandwidth,
    Rcpp::Named("weights") = c_weights,
    Rcpp::Named("kernel") = kernel,
    Rcpp::Named("preserve_var") = preserve_var
  );

}

