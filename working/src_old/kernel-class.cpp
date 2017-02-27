

#include "kernels.h"
#include "kernel-class.h"
#include <Rcpp.h>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::_;


UnivarKernel::UnivarKernel() {
  data = 0.0;
  n = 1;
  mean = 0.0;
  var = 1.0;
  bandwidth = 1.0;
  weights = 1.0;
  cumul_weights = 1.0;
  pdf_fun = dens_gauss;
  rng_fun = R::norm_rand;
}


UnivarKernel::UnivarKernel(const NumericVector& x) {
  data = x;
  n = x.length();
  mean = Rcpp::mean(x);
  var = Rcpp::var(x);
  bandwidth = 1.06 * std::sqrt(var) * std::pow(n, -1.0/5.0); // rule of thumb
  double p = 1.0/static_cast<double>(n);
  weights = NumericVector(n, p);
  cumul_weights[0] = p;
  for (int i = 1; i < n; i++)
    cumul_weights[i] = cumul_weights[i-1] + p;
  pdf_fun = dens_gauss;
  rng_fun = R::norm_rand;
}

UnivarKernel::UnivarKernel(const NumericVector& x, const NumericVector& w) {
  data = x;
  n = x.length();
  mean = Rcpp::mean(x);
  var = Rcpp::var(x);
  bandwidth = 1.06 * std::sqrt(var) * std::pow(n, -1.0/5.0); // rule of thumb
  set_weights(w);
  pdf_fun = dens_gauss;
  rng_fun = R::norm_rand;
}

void UnivarKernel::set_pdf_fun(double (*fun)(double x, double bw)) {
  pdf_fun = fun;
}

void UnivarKernel::set_rng_fun(double (*fun)()) {
  rng_fun = fun;
}

void UnivarKernel::set_bandwidth(double bw) {
  if (bw <= 0.0)
    Rcpp::stop("bandwidth needs to be greater than zero");
  bandwidth = bw;
}

void UnivarKernel::set_weights(const NumericVector& w) {

  if (w.length() != n)
    Rcpp::stop("weights should be of the same size as data");

  weights = Rcpp::clone(w);
  cumul_weights = Rcpp::clone(w);

  for (int i = 1; i < n; i++)
    cumul_weights[i] += cumul_weights[i-1];

  double total_w = cumul_weights[n-1];
  if (std::abs(total_w - 1.0) > 1e-8) {
    for (int i = 0; i < n; i++) {
      weights[i] /= total_w;
      cumul_weights[i] /= total_w;
    }
  }

}

NumericVector UnivarKernel::pdf(const NumericVector& x) {
  int m = x.length();
  NumericVector p(m);
  for (int i = 0; i < m; i++) {
    if (ISNAN(x[i])) {
      p[i] = x[i];
      continue;
    }
    p[i] = 0.0;
    for (int j = 0; j < n; j++)
      p[i] += pdf_fun(x[i] - data[j], bandwidth) * weights[j];
  }
  return p;
}

NumericVector UnivarKernel::rng(int m, bool preserve_var) {
  int j;
  double u;
  NumericVector x(m);

  if (preserve_var) {

    double c = sqrt(1.0 + pow(bandwidth, 2.0)/var);

    for (int i = 0; i < m; i++) {
      u = rng_unif();
      for (j = 0; j < n; j++) {
        if (cumul_weights[j] >= u)
          break;
      }
      x[i] = mean + (data[j] - mean + rng_fun() * bandwidth) / c;
    }

  } else {

    for (int i = 0; i < m; i++) {
      u = rng_unif();
      for (j = 0; j < n; j++) {
        if (cumul_weights[j] >= u)
          break;
      }
      x[i] = data[j] + rng_fun() * bandwidth;
    }

  }

  return x;
}


// Multivariate Kenrnel



MultivarKernel::MultivarKernel(const NumericMatrix& x) {
  data = x;
  n = x.nrow();
  k = x.ncol();
  NumericVector mean(k), var(k), bandwidth(k),
                weights(n), cumul_weights(n);

  for (int i = 0; i < k; i++) {
    mean[i] = Rcpp::mean(x(_, i));
    var[i] = Rcpp::var(x(_, i));
    bandwidth[i] = 1.06 * std::sqrt(var[i]) * std::pow(n, -1.0/5.0); // rule of thumb
  }

  double p = 1.0/static_cast<double>(n);
  weights = NumericVector(n, p);
  cumul_weights[0] = p;
  for (int i = 1; i < n; i++)
    cumul_weights[i] = cumul_weights[i-1] + p;

  pdf_fun = dens_gauss;
  rng_fun = R::norm_rand;
}

MultivarKernel::MultivarKernel(const NumericMatrix& x, const NumericVector& w) {
  data = x;
  n = x.nrow();
  k = x.ncol();
  NumericVector mean(k), var(k), bandwidth(k),
                weights(n), cumul_weights(n);

  for (int i = 0; i < k; i++) {
    mean[i] = Rcpp::mean(x(_, i));
    var[i] = Rcpp::var(x(_, i));
    bandwidth[i] = 1.06 * std::sqrt(var[i]) * std::pow(n, -1.0/5.0); // rule of thumb
  }

  set_weights(w);

  pdf_fun = dens_gauss;
  rng_fun = R::norm_rand;
}

void MultivarKernel::set_pdf_fun(double (*fun)(double x, double bw)) {
  pdf_fun = fun;
}

void MultivarKernel::set_rng_fun(double (*fun)()) {
  rng_fun = fun;
}

void MultivarKernel::set_bandwidth(const NumericVector& bw) {
  if (bw.length() != k)
    Rcpp::stop("");
  // if (Rcpp::any(bw <= 0.0))
  //   Rcpp::stop("bandwidth needs to be greater than zero");
  bandwidth = Rcpp::clone(bw);
}

void MultivarKernel::set_weights(const NumericVector& w) {

  if (w.length() != n)
    Rcpp::stop("weights should be of the same size as data");

  weights = Rcpp::clone(w);
  cumul_weights = Rcpp::clone(w);

  for (int i = 1; i < n; i++)
    cumul_weights[i] += cumul_weights[i-1];

  double total_w = cumul_weights[n-1];
  if (std::abs(total_w - 1.0) > 1e-8) {
    for (int i = 0; i < n; i++) {
      weights[i] /= total_w;
      cumul_weights[i] /= total_w;
    }
  }

}

NumericMatrix MultivarKernel::rng(int m, bool preserve_var) {
  int j;
  double u;
  NumericMatrix x(n, k);

  if (preserve_var) {

    NumericVector c(k);
    for (int i = 0; i < k; i++)
      c[i] = sqrt(1.0 + pow(bandwidth[i], 2.0)/var[i]);

    for (int i = 0; i < m; i++) {
      u = rng_unif();
      for (j = 0; j < n; j++) {
        if (cumul_weights[j] >= u)
          break;
      }
      for (int h = 0; h < k; h++)
        x(i, h) = mean[h] + (data(j, h) - mean[h] + rng_fun() * bandwidth[h]) / c[h];
    }

  } else {

    for (int i = 0; i < m; i++) {
      u = rng_unif();
      for (j = 0; j < n; j++) {
        if (cumul_weights[j] >= u)
          break;
      }
      for (int h = 0; h < k; h++)
        x(i, h) = data(j, h) + rng_fun() * bandwidth[h];
    }

  }

  return x;
}
