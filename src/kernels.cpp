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

// kernel RNGs

double rng_epan() {
  double u1, u2, u3, au1, au2, au3;
  u1 = (rng_unif() * 2.0 - 1.0);
  u2 = (rng_unif() * 2.0 - 1.0);
  u3 = (rng_unif() * 2.0 - 1.0);
  au1 = abs(u1);
  au2 = abs(u2);
  au3 = abs(u3);

  if (au3 >= au2 && au3 >= au1) {
    return u2 * SQRT_5;
  } else {
    return u3 * SQRT_5;
  }
}

double rng_cosine() {
  return (R::rbeta(3.3575, 3.3575) * 2.0 - 1.0) * 2.766159483867713042571;
}

double rng_optcos() {
  double u;
  u = (rng_unif() * 2.0 - 1.0);
  return (2.0 * acos(u)/M_PI - 1.0) * 2.297603117487196922042;
}

double rng_triang() {
  double u, v;
  u = rng_unif();
  v = rng_unif();
  return (u-v) * SQRT_6;
}

double rng_rect() {
  return (rng_unif() * 2.0 - 1.0) * SQRT_3;
}

double rng_biweight() {
  return (R::rbeta(3.0, 3.0) * 2.0 - 1.0) * SQRT_7;
}

double rng_triweight() {
  return (R::rbeta(4.0, 4.0) * 2.0 - 1.0) * 3.0;
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
//' @aliases rcosine
//' @aliases rbiweight
//' @aliases rtriweight
//' @export
// [[Rcpp::export]]
NumericVector rempan(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = rng_epan();
  return x;
}

//' @rdname KernelRNG
//' @export
// [[Rcpp::export]]
NumericVector rcosine(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = rng_cosine();
  return x;
}

//' @rdname KernelRNG
//' @export
// [[Rcpp::export]]
NumericVector roptcos(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = rng_optcos();
  return x;
}

//' @rdname KernelRNG
//' @export
// [[Rcpp::export]]
NumericVector rtriang(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = rng_triang();
  return x;
}

//' @rdname KernelRNG
//' @export
// [[Rcpp::export]]
NumericVector rrect(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = rng_rect();
  return x;
}

//' @rdname KernelRNG
//' @export
// [[Rcpp::export]]
NumericVector rbiweight(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = rng_biweight();
  return x;
}

//' @rdname KernelRNG
//' @export
// [[Rcpp::export]]
NumericVector rtriweight(int n) {
  NumericVector x(n);
  for (int i = 0; i < n; i++)
    x[i] = rng_triweight();
  return x;
}


// [[Rcpp::export]]
NumericMatrix add_noise(
    const NumericMatrix &x,
    const std::string &kernel = "gaussian",
    const NumericVector &bandwidth = 1.0,
    const NumericVector &mean = 1.0,
    const NumericVector &var = 0.0,
    const bool &preserve_var = false
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

  int n = x.nrow();
  int k = x.ncol();
  NumericMatrix out(n, k);
  NumericVector bw(k);
  int bn = bandwidth.length();

  if (!preserve_var) {

    for (int j = 0; j < k; j++) {
      bw[j] = bandwidth[j % bn];
      if (bw[j] < 0)
        Rcpp::stop("bandwidth < 0");
    }

    for (int i = 0; i < n; i++)
      for (int j = 0; j < k; j++)
        out(i, j) = x(i, j) * rng_kern() * bw[j];

  } else {

    NumericVector bw(k);
    NumericVector m(k);
    NumericVector s(k);
    NumericVector norm_const(k);

    int mn = mean.length();
    int sn = var.length();

    for (int j = 0; j < k; j++) {
      bw[j] = bandwidth[j % bn];
      if (bw[j] < 0)
        Rcpp::stop("bandwidth < 0");
      m[j] = mean[j % mn];
      s[j] = var[j % sn];
      if (s[j] < 0)
        Rcpp::stop("variance < 0");
      norm_const[j] = sqrt(1.0 + pow(bw[j], 2.0)/s[j]);
    }

    for (int i = 0; i < n; i++)
      for (int j = 0; j < k; j++)
        out(i, j) = m[j] + (x(i, j) - m[j] + rng_kern() * bw[j]) / norm_const[j];

  }

  return out;

}

