#include <Rcpp.h>
#include "kernels.h"

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


// [[Rcpp::export]]
NumericVector cpp_rkernel(
    int n,
    const NumericVector &data,
    const std::string &kernel = "gaussian",
    const NumericVector &bw = 1.0,
    const NumericVector &adjust = 1.0
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

  NumericVector x(n);
  int k = data.length();
  int nbwh = bw.length();
  int nadj = adjust.length();

  for (int i = 0; i < n; i++) {
    x[i] = rng_kern() * bw[i % nbwh] * adjust[i % nadj];
    x[i] += data[ static_cast<int>(rng_unif() * static_cast<double>(k)) ];
  }

  return x;

}



// [[Rcpp::export]]
NumericVector cpp_dkernel(
    const NumericVector &x,
    const NumericVector &data,
    const std::string &kernel = "gaussian",
    const NumericVector &bw = 1.0,
    const NumericVector &adjust = 1.0,
    bool log_prob = false
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

  int n = x.length();
  int k = data.length();
  int nbwh = bw.length();
  int nadj = adjust.length();
  NumericVector p(n, 0.0);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++)
      p[i] += dens_kern(x[i] - data[j], bw[i % nbwh] * adjust[i % nadj]);
    p[i] /= static_cast<double>(k);
  }

  if (log_prob)
    p = Rcpp::log(p);

  return p;

}
