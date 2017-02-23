#include <Rcpp.h>
#include "kernels.h"

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


// [[Rcpp::export]]
NumericMatrix add_noise(
    const NumericMatrix &x,
    const std::string &kernel = "gaussian",
    const NumericVector &bandwidth = 1.0,
    const NumericVector &mean = 1.0,
    const NumericVector &var = 0.0,
    const bool &preserve_var = false
  ) {

  double (*rng_fun)();

  if (kernel == "rectangular") {
    rng_fun = rng_rect;
  } else if (kernel == "triangular") {
    rng_fun = rng_triang;
  } else if (kernel == "biweight") {
    rng_fun = rng_biweight;
  } else if (kernel == "triweight") {
    rng_fun = rng_triweight;
  } else if (kernel == "cosine") {
    rng_fun = rng_cosine;
  } else if (kernel == "optcosine") {
    rng_fun = rng_optcos;
  } else if (kernel == "epanechnikov") {
    rng_fun = rng_epan;
  } else {
    rng_fun = R::norm_rand;
  }

  int n = x.nrow();
  int k = x.ncol();
  NumericMatrix out(n, k);
  NumericVector bw(k);
  int bn = bandwidth.length();

  if (!preserve_var) {

    for (int j = 0; j < k; j++) {
      bw[j] = bandwidth[j % bn];
      if (bw[j] < 0.0)
        Rcpp::stop("bandwidth < 0");
    }

    for (int i = 0; i < n; i++)
      for (int j = 0; j < k; j++)
        out(i, j) = x(i, j) + rng_fun() * bw[j];

  } else {

    NumericVector bw(k);
    NumericVector m(k);
    NumericVector s(k);
    NumericVector norm_const(k);

    int mn = mean.length();
    int sn = var.length();

    for (int j = 0; j < k; j++) {
      bw[j] = bandwidth[j % bn];
      if (bw[j] < 0.0)
        Rcpp::stop("bandwidth < 0");
      m[j] = mean[j % mn];
      s[j] = var[j % sn];
      if (s[j] < 0.0)
        Rcpp::stop("variance < 0");
      norm_const[j] = sqrt(1.0 + pow(bw[j], 2.0)/s[j]);
    }

    for (int i = 0; i < n; i++)
      for (int j = 0; j < k; j++)
        out(i, j) = m[j] + (x(i, j) - m[j] + rng_fun() * bw[j]) / norm_const[j];

  }

  return out;

}


// [[Rcpp::export]]
NumericVector rng_kern(
    int n,
    const std::string &kernel = "gaussian"
  ) {

  double (*rng_fun)();

  if (kernel == "rectangular") {
    rng_fun = rng_rect;
  } else if (kernel == "triangular") {
    rng_fun = rng_triang;
  } else if (kernel == "biweight") {
    rng_fun = rng_biweight;
  } else if (kernel == "triweight") {
    rng_fun = rng_triweight;
  } else if (kernel == "cosine") {
    rng_fun = rng_cosine;
  } else if (kernel == "optcosine") {
    rng_fun = rng_optcos;
  } else if (kernel == "epanechnikov") {
    rng_fun = rng_epan;
  } else {
    rng_fun = R::norm_rand;
  }

  NumericVector x(n);

  for (int i = 0; i < n; i++)
    x[i] = rng_fun();

  return x;


}
