#include <Rcpp.h>

using Rcpp::NumericVector;

typedef void (*fun_ptr)(void);


double dens_gauss(double x, double bw) {
  return R::dnorm(x, 0.0, bw, false);
}



class UvKernel
{
  
private:
  NumericVector data;
  int n;
  double bandwidth;
  NumericVector weights;
  double (*pdf_fun)(double x, double bw);
  double (*rng_fun)(void);
  
public:
  UvKernel();
  UvKernel(const NumericVector& x);
  void set_pdf_fun(double (*fun)(double x, double bw));
  void set_rng_fun(double (*fun)());
  void set_bandwidth(double bw);
  NumericVector pdf(const NumericVector& x);
  NumericVector rng(int m);
  
};

UvKernel::UvKernel() {
  data = 0.0;
  n = 1;
  bandwidth = 1.0;
  weights = 1;
  pdf_fun = dens_gauss;
  rng_fun = R::norm_rand;
}


UvKernel::UvKernel(const NumericVector& x) {
  data = x;
  n = x.length();
  bandwidth = 1.06 * Rcpp::sd(x) * std::pow(n, -1.0/5.0); // rule of thumb
  weights = NumericVector(n, 1.0/n);
  pdf_fun = dens_gauss;
  rng_fun = R::norm_rand;
}

void UvKernel::set_pdf_fun(double (*fun)(double x, double bw)) {
  pdf_fun = fun;
}

void UvKernel::set_rng_fun(double (*fun)()) {
  rng_fun = fun;
}

void UvKernel::set_bandwidth(double bw) {
  if (bw <= 0.0)
    Rcpp::stop("bandwidth need to be >0");
  bandwidth = bw;
}

NumericVector UvKernel::pdf(const NumericVector& x) {
  NumericVector p(n);
  for (int i = 0; i < n; i++) {
    if (ISNAN(x[1]))
      p[i] = x[i];
    else
      p[i] = pdf_fun(x[i], bandwidth);
  }
  return p;
}

NumericVector UvKernel::rng(int m) {
  NumericVector x(n);
  for (int i = 0; i < m; i++) {
    
  }
  return x;
}
