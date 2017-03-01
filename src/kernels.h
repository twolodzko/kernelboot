
#ifndef KERNELS_H
#define KERNELS_H

#include "shared.h"

// kernel RNGs

double rng_epan();
double rng_cosine();
double rng_optcos();
double rng_triang();
double rng_rect();
double rng_biweight();
double rng_triweight();

// kernel density functions

double dens_epan(double x, double bw);
double dens_cosine(double x, double bw);
double dens_optcos(double x, double bw);
double dens_triang(double x, double bw);
double dens_rect(double x, double bw);
double dens_biweight(double x, double bw);
double dens_triweight(double x, double bw);
double dens_gauss(double x, double bw);


#endif


