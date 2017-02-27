
#ifndef KDE_CLASS_H
#define KDE_CLASS_H

#include <RcppArmadillo.h>
#include "kernels.h"


template <class T>
class kde
{

private:
  T data;
  int n, k;
  arma::vec weights;
  arma::vec bandwidth;
  double (*pdf_fun)(double x, double bw);
  double (*rng_fun)(void);

public:
  kde();
  kde(const arma::vec& x);
  kde(const arma::mat& x);
  void change_pdf_fun(double (*fun)(double x, double bw));
  void change_rng_fun(double (*fun)());
  void change_bandwidth(const arma::vec& bw);
  void change_weights(const arma::vec& w);
  arma::vec pdf(const arma::vec& x) const;
  arma::vec pdf(const arma::mat& x) const;
  T rng(int m, bool preserve_var = false) const;

};

# endif

