#include <Rcpp.h>

#define SQRT_3 1.732050807568877193177
#define SQRT_5 2.236067977499789805051
#define SQRT_6 2.449489742783177881336
#define SQRT_7 2.645751311064590716171

using std::pow;
using std::sqrt;
using std::abs;
using std::exp;
using std::log;
using std::floor;
using std::ceil;
using std::sin;
using std::cos;
using std::tan;
using std::atan;
using std::acos;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}


//' Random generation from kernels
//'
//' @param n number of observations.
//'
//' @name KernelRNG
//' @aliases KernelRNG
//' @aliases rempan
//' @aliases rtriang
//' @aliases rrect
//' @aliases roptcos
//' @aliases rbiweight
//' @aliases rtriweight

// [[Rcpp::export]]
NumericVector rempan(int n) {

  NumericVector x(n);
  double u1, u2, u3, au1, au2, au3;

  for (int i = 0; i < n; i++) {

    u1 = R::runif(-1.0, 1.0);
    u2 = R::runif(-1.0, 1.0);
    u3 = R::runif(-1.0, 1.0);
    au1 = abs(u1);
    au2 = abs(u2);
    au3 = abs(u3);

    if (au3 >= au2 && au3 >= au1) {
      x[i] = u2;
    } else {
      x[i] = u3;
    }

    x[i] *= SQRT_5;
  }

  return x;
}

//' @rdname KernelRNG
// [[Rcpp::export]]
NumericVector roptcos(int n) {
  NumericVector x(n);
  double u;
  for (int i = 0; i < n; i++) {
    u = R::runif(-1.0, 1.0);
    x[i] = (2.0 * acos(u)/M_PI - 1.0) * 2.297603117487196922042;
  }
  return x;
}


//' @rdname KernelRNG
// [[Rcpp::export]]
NumericVector rtriang(int n) {
  NumericVector x(n);
  double u, v;
  for (int i = 0; i < n; i++) {
    u = rng_unif();
    v = rng_unif();
    x[i] = (u-v) * SQRT_6;
  }
  return x;
}


//' @rdname KernelRNG
// [[Rcpp::export]]
NumericVector rrect(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = R::runif(-1.0, 1.0) * SQRT_3;
  return x;
}


//' @rdname KernelRNG
// [[Rcpp::export]]
NumericVector rbiweight(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = (R::rbeta(3.0, 3.0) * 2.0 - 1.0) * SQRT_7;
  return x;
}


//' @rdname KernelRNG
// [[Rcpp::export]]
NumericVector rtriweight(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = (R::rbeta(4.0, 4.0) * 2.0 - 1.0) * 3.0;
  return x;
}


// [[Rcpp::export]]
NumericMatrix rsmvnorm(int n, NumericVector sigma) {

  int k = sigma.length();
  NumericMatrix x(n, k);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < k; j++)
      x(i, j) = R::rnorm(0.0, sigma[j]);

  return x;
}

