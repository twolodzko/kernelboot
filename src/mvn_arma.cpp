
#include <RcppArmadillo.h>

/*
 * The following code comes from Rcpp Gallery articles:
 *
 * Nino Hardt, Dicko Ahmadou (Jul 13, 2013).
 * Faster Multivariate Normal densities with RcppArmadillo and OpenMP.
 * http://gallery.rcpp.org/articles/dmvnorm_arma/
 *
 * and
 *
 * Ahmadou Dicko (Mar 12, 2013).
 * Generating a multivariate gaussian distribution using RcppArmadillo.
 * http://gallery.rcpp.org/articles/simulate-multivariate-normal/
 *
 * The code comes under GPLv2 License
 * (https://www.gnu.org/licenses/gpl-2.0.html)
 *
 */


// [[Rcpp::export]]
arma::vec cpp_dmvn(
    const arma::mat& x,
    const arma::rowvec& mu,
    const arma::mat& sigma,
    const bool& log_prob = false
  ) {

  const unsigned int n = x.n_rows;
  const unsigned int m = x.n_cols;
  arma::vec p(n);

  try {

    const arma::mat chol_sigma = arma::chol(sigma);
    const arma::mat rooti = arma::trans(arma::inv(arma::trimatu(chol_sigma)));
    const double rootisum = arma::sum(arma::log(rooti.diag()));
    const double c = -(static_cast<double>(m) / 2.0) * M_LN_2PI;

    arma::vec z;
    for (unsigned int i = 0; i < n; i++) {
      z = rooti * arma::trans( x.row(i) - mu ) ;
      p[i] = c - 0.5 * arma::sum(z % z) + rootisum;
    }

    if (!log_prob)
      p = arma::exp(p);

    return p;

  } catch ( std::exception& __ex__ ) {
    forward_exception_to_r(__ex__);
  } catch (...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  p.fill(NAN);
  return p;

}


// [[Rcpp::export]]
arma::mat cpp_rmvn(
    const unsigned int& n,
    const arma::vec& mu,
    const arma::mat& sigma
  ) {

  const unsigned int m = sigma.n_cols;
  arma::mat res(n, m);

  try {

    const arma::mat Y = arma::randn(n, m);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);

  } catch ( std::exception& __ex__ ) {
    forward_exception_to_r(__ex__);
  } catch (...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  res.fill(NA_REAL);
  return res;

}

