
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "kernels.h"
#include "shared.h"


// [[Rcpp::export]]
Rcpp::List cpp_rmvk(
    const unsigned int& n,
    const arma::mat& y,
    const arma::vec& bandwidth,
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

  const unsigned int m = y.n_cols;
  const unsigned int k = y.n_rows;
  arma::mat samp(n, m);
  arma::vec c_weights(k);
  std::vector<unsigned int> idx(n);

  if (bandwidth.n_elem != m)
    Rcpp::stop("dimmensions of y and bandwidth do not match");

  if (arma::any(bandwidth <= 0.0))
    Rcpp::stop("bandwidth needs to be positive");

  if (!arma::is_finite(bandwidth))
    Rcpp::stop("inappropriate values of bandwidth");

  if (arma::any(weights < 0.0))
    Rcpp::stop("weights need to be non-negative");

  if (!arma::is_finite(weights))
    Rcpp::stop("inappropriate values of weights");

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

    if (k == 1) {

      for (unsigned int i = 0; i < n; i++) {
        for (unsigned int l = 0; l < m; l++)
          samp(i, l) = y(0, l) + rng_kern() * bandwidth[l];
      }

    } else if (!shrinked) {

      unsigned int j;
      for (unsigned int i = 0; i < n; i++) {
        j = sample_int(c_weights);
        idx[i] = j + 1;
        for (unsigned int l = 0; l < m; l++)
          samp(i, l) = y(j, l) + rng_kern() * bandwidth[l];
      }

    } else {

      const arma::rowvec my = arma::mean(y);
      const arma::rowvec sy = arma::var(y);
      arma::rowvec bw_sq(m);
      for (unsigned int i = 0; i < m; i++)
        bw_sq[i] = bandwidth[i] * bandwidth[i];
      const arma::rowvec c = arma::sqrt(1.0 + bw_sq/sy);

      unsigned int j;
      for (unsigned int i = 0; i < n; i++) {
        j = sample_int(c_weights);
        idx[i] = j + 1;
        for (unsigned int l = 0; l < m; l++)
          samp(i, l) = my[l] + (y(j, l) - my[l] + rng_kern() * bandwidth[l]) / c[l];

      }

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
    Rcpp::Named("kernel") = kernel,
    Rcpp::Named("shrinked") = shrinked
  );

}

