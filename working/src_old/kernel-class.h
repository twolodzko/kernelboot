
#ifndef KERNEL_CLASS_H
#define KERNEL_CLASS_H

#include "kernels.h"
#include <Rcpp.h>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::_;


class UnivarKernel
{

private:
  NumericVector data;
  int n;
  double mean, var, bandwidth;
  NumericVector weights;
  NumericVector cumul_weights;
  double (*pdf_fun)(double x, double bw);
  double (*rng_fun)(void);

public:
  UnivarKernel();
  UnivarKernel(const NumericVector& x);
  UnivarKernel(const NumericVector& x, const NumericVector& w);
  void set_pdf_fun(double (*fun)(double x, double bw));
  void set_rng_fun(double (*fun)());
  void set_bandwidth(double bw);
  void set_weights(const NumericVector& w);
  NumericVector pdf(const NumericVector& x);
  NumericVector rng(int m, bool preserve_var = false);

};


class MultivarKernel
{

private:
  NumericMatrix data;
  int n;
  int k;
  NumericVector mean, var, bandwidth, weights, cumul_weights;
  double (*pdf_fun)(double x, double bw);
  double (*rng_fun)(void);

public:
  MultivarKernel(const NumericMatrix& x);
  MultivarKernel(const NumericMatrix& x, const NumericVector& w);
  void set_pdf_fun(double (*fun)(double x, double bw));
  void set_rng_fun(double (*fun)());
  void set_bandwidth(const NumericVector& bw);
  void set_weights(const NumericVector& w);
  //NumericMatrix pdf(const NumericMatrix& x);
  NumericMatrix rng(int m, bool preserve_var = false);

};


#endif


