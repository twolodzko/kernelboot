#include <Rcpp.h>
#include "kernels.h"
#include "kernel-class.h"

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


// [[Rcpp::export]]
NumericVector cpp_ruvkern(
    int n,
    const NumericVector& data,
    const NumericVector& weights,
    const std::string& kernel = "gaussian",
    double bw = 1.0,
    bool preserve_var = false
  ) {

  UnivarKernel kern(data, weights);
  kern.set_bandwidth(bw);
  kern.set_rng_fun(get_rng_fun(kernel));
  return kern.rng(n, preserve_var);

}


// [[Rcpp::export]]
NumericVector cpp_duvkern(
    const NumericVector& x,
    const NumericVector& data,
    const NumericVector& weights,
    const std::string& kernel = "gaussian",
    double bw = 1.0,
    bool log_prob = false
  ) {

  UnivarKernel kern(data, weights);
  kern.set_bandwidth(bw);
  kern.set_pdf_fun(get_pdf_fun(kernel));

  if (log_prob)
    return Rcpp::log(kern.pdf(x));

  return kern.pdf(x);

}
