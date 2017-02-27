
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]


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
    arma::mat x,  
    arma::rowvec mu,  
    arma::mat sigma, 
    bool log_prob = false
  ) { 
  
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(k) / 2.0) * M_LN_2PI;
  
  for (int i = 0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mu) ;    
    out(i) = constants - 0.5 * arma::sum(z % z) + rootisum;     
  }  
  
  if (!log_prob)
    out = exp(out);

  return out;
}


// [[Rcpp::export]]
arma::mat cpp_rmvn(
    int n,
    arma::vec mu,
    arma::mat sigma
  ) {
  
  int k = sigma.n_cols;
  arma::mat Y = arma::randn(n, k);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

