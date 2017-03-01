
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include "kernels.h"
#include "shared.h"


// [[Rcpp::export]]
Rcpp::List cpp_dmvkde(
    const arma::mat& y,
    const arma::mat& x,
    const arma::mat& bandwidth,
    const arma::vec& weights,
    const bool& log_prob = false,
    const bool& is_chol = false
  ) {

  const unsigned int n = y.n_rows;
  const unsigned int m = x.n_cols;
  const unsigned int k = x.n_rows;
  arma::vec p(n), c_weights(k);

  if (y.n_cols != m)
    Rcpp::stop("wrong dimmensions of y");

  if (bandwidth.n_cols != bandwidth.n_rows || bandwidth.n_cols != m)
    Rcpp::stop("wrong dimmensions of bandwidth");

  if (any(weights < 0.0))
    Rcpp::stop("weights need to be non-negative");

  try {

    if (weights.n_elem == 1) {
      c_weights.fill( 1.0/static_cast<double>(k) );
    } else {
      if (weights.n_elem != k)
        Rcpp::stop("length(weights) != nrow(x)");
      c_weights = weights;
    }

    c_weights /= sum(c_weights);

    arma::mat bw_chol;
    if (is_chol) {
      bw_chol = bandwidth;
    } else {
      bw_chol = arma::chol(bandwidth);
    }

    arma::mat rooti = arma::trans(arma::inv(arma::trimatu(bw_chol)));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(k) / 2.0) * M_LN_2PI;

    arma::vec z;
    double tmp;
    for (unsigned int i = 0; i < n; i++) {
      p[i] = 0.0;
      for (unsigned int j = 0; j < k; j++) {
        z = rooti * arma::trans( y.row(i) - x.row(j) ) ;
        tmp = constants - 0.5 * arma::sum(z % z) + rootisum;
        tmp += log(c_weights[j]);
        p[i] += exp(tmp);
      }
    }

    if (log_prob)
      p = log(p);

  } catch ( std::exception& __ex__ ) {
    forward_exception_to_r(__ex__);
  } catch (...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return Rcpp::List::create(
    Rcpp::Named("density") = p,
    Rcpp::Named("x") = x,
    Rcpp::Named("bandwidth") = bandwidth,
    Rcpp::Named("weights") = c_weights,
    Rcpp::Named("log_prob") = log_prob
  );

}


// [[Rcpp::export]]
Rcpp::List cpp_rmvkde(
    const unsigned int& n,
    const arma::mat& x,
    const arma::mat& bandwidth,
    const arma::vec& weights,
    const bool& is_chol = false
  ) {

  const unsigned int m = x.n_cols;
  const unsigned int k = x.n_rows;
  arma::mat samp(n, m);
  arma::vec c_weights(k);
  std::vector<unsigned int> idx(n);

  if (bandwidth.n_cols != bandwidth.n_rows || bandwidth.n_cols != m)
    Rcpp::stop("wrong dimmensions of bandwidth");

  if (any(weights < 0.0))
    Rcpp::stop("weights need to be non-negative");

  try {

    if (weights.n_elem == 1) {
      c_weights.fill( 1.0/static_cast<double>(k) );
    } else {
      if (weights.n_elem != k)
        Rcpp::stop("length(weights) != nrow(x)");
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
      means.row(i) = x.row(j);
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
    Rcpp::Named("sample") = samp,
    Rcpp::Named("boot_index") = idx,
    Rcpp::Named("x") = x,
    Rcpp::Named("bandwidth") = bandwidth,
    Rcpp::Named("weights") = c_weights
  );

}

