
#include "kernels.h"
#include "shared.h"

// Constants

#define SQRT_3 1.732050807568877193177
#define SQRT_5 2.236067977499789805051
#define SQRT_6 2.449489742783177881336
#define SQRT_7 2.645751311064590716171

// kernel RNGs

double rng_epan() {
  double u1, u2, u3, au1, au2, au3;
  u1 = (rng_unif() * 2.0 - 1.0);
  u2 = (rng_unif() * 2.0 - 1.0);
  u3 = (rng_unif() * 2.0 - 1.0);
  au1 = std::abs(u1);
  au2 = std::abs(u2);
  au3 = std::abs(u3);

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
  return (2.0 * std::acos(u)/M_PI - 1.0) * 2.297603117487196922042;
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


// kernel density functions

double dens_epan(double x, double bw) {
  double a = bw * SQRT_5;
  double ax = std::abs(x);
  if (ax > a)
    return 0.0;
  return 0.75 * (1.0 - std::pow(ax/a, 2)) / a;
}

double dens_cosine(double x, double bw) {
  double a = bw * 0.36151205519132795; // sqrt(1/3 - 2/pi^2)
  if (x < -a || x > a)
    return 0.0;
  return (1.0 + std::cos(M_PI * x/a)) / (2.0 * a);
}

double dens_optcos(double x, double bw) {
  double a = bw * 0.43523617825417249; // sqrt(1 - 8/pi^2)
  if (x < -a || x > a)
    return 0.0;
  return M_PI/4.0 * std::cos(M_PI * x/(2.0 * a)) / a;
}

double dens_triang(double x, double bw) {
  double a = bw * SQRT_6;
  double ax = std::abs(x);
  if (ax > a)
    return 0.0;
  return (1.0 - ax/a) / a;
}

double dens_rect(double x, double bw) {
  double a = bw * SQRT_3;
  if (x < -a || x > a)
    return 0.0;
  return 0.5 / a;
}

double dens_biweight(double x, double bw) {
  double a = bw * SQRT_7;
  double ax = std::abs(x);
  if (ax > a)
    return 0.0;
  return 0.9375 * std::pow(1.0 - std::pow(ax/a, 2.0), 2.0) / a; // 15/16
}

double dens_triweight(double x, double bw) {
  double a = bw * 3.0;
  if (x < -a || x > a)
    return 0.0;
  return 1.09375 * std::pow(1 - std::pow(x, 2.0), 3.0) / a; // 35/32
}

double dens_gauss(double x, double bw) {
  return R::dnorm(x, 0.0, bw, false);
}

