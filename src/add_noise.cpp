
#include "kernels.h"
#include "kernel-class.h"
#include <Rcpp.h>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::_;


// [[Rcpp::export]]
NumericMatrix add_noise(
    const NumericMatrix& x,
    const NumericVector& weights,
    const std::string& kernel = "gaussian",
    const NumericVector& bandwidth = 1.0,
    const bool& preserve_var = false
  ) {

  if (weights.length() != x.nrow())
    Rcpp::stop("");

  if (bandwidth.length() != x.ncol())
    Rcpp::stop("");

  MultivarKernel kern(x, weights);
  kern.set_rng_fun(get_rng_fun(kernel));
  kern.set_bandwidth(bandwidth);

  return kern.rng(x.nrow(), preserve_var);

  /*
  std::vector<UnivarKernel> kern(k);
  for (int i = 0; i < k; i++) {
    kern[i] = UnivarKernel(x(_, i), weights); // !!!!
    kern[i].set_rng_fun(rng_kern);
    kern[i].set_bandwidth(bandwidth[i]);
    xx(_, i) = kern[i].rng(n, preserve_var);
  }
   */

}

