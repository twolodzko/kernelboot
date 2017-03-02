
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "kernels.h"
#include "shared.h"


// [[Rcpp::export]]
Rcpp::List cpp_dmvk(
    const arma::mat& x,
    const arma::mat& y,
    const arma::mat& bandwidth,
    const arma::vec& weights,
    const bool& log_prob = false,
    const bool& is_chol = false    // bw is passed as Cholesky decomposition
  ) {

  const unsigned int n = x.n_rows;
  const unsigned int m = y.n_cols;
  const unsigned int k = y.n_rows;
  arma::vec p(n), c_weights(k);

  if (x.n_cols != m || bandwidth.n_cols != m)
    Rcpp::stop("dimmensions of x, y and bandwidth do not match");

  if (bandwidth.n_cols != bandwidth.n_rows)
    Rcpp::stop("bandwidth is not a square matrix");

  if (any(weights < 0.0))
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

    arma::mat bw_chol;
    if (is_chol) {
      bw_chol = bandwidth;
    } else {
      bw_chol = arma::chol(bandwidth);
    }

    const arma::mat rooti = arma::trans(arma::inv(arma::trimatu(bw_chol)));
    const double rootisum = arma::sum(arma::log(rooti.diag()));
    const double c = -(static_cast<double>(m) / 2.0) * M_LN_2PI;

    arma::vec z;
    double tmp;
    for (unsigned int i = 0; i < n; i++) {
      p[i] = 0.0;
      for (unsigned int j = 0; j < k; j++) {
        z = rooti * arma::trans( x.row(i) - y.row(j) ) ;
        tmp = c - 0.5 * arma::sum(z % z) + rootisum;
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
    Rcpp::Named("log_prob") = log_prob
  );

}


// [[Rcpp::export]]
Rcpp::List cpp_rmvk(
    const unsigned int& n,
    const arma::mat& y,
    const arma::mat& bandwidth,
    const arma::vec& weights,
    const bool& is_chol = false    // bw is passed as Cholesky decomposition
  ) {

  const unsigned int m = y.n_cols;
  const unsigned int k = y.n_rows;
  arma::mat samp(n, m);
  arma::vec c_weights(k);
  std::vector<unsigned int> idx(n);

  if (bandwidth.n_cols != m)
    Rcpp::stop("dimmensions of y and bandwidth do not match");

  if (bandwidth.n_cols != bandwidth.n_rows)
    Rcpp::stop("bandwidth is not a square matrix");

  if (any(weights < 0.0))
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

    arma::mat bw_chol;
    if (is_chol) {
      bw_chol = bandwidth;
    } else {
      bw_chol = arma::chol(bandwidth);
    }

    samp = arma::randn(n, m) * bw_chol;
    arma::mat means(n, m);

    unsigned int j;
    for (unsigned int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      means.row(i) = y.row(j);
    }

    samp += means;

    for (unsigned int i = k; i > 0; i--)
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
    Rcpp::Named("weights") = c_weights
  );

}

