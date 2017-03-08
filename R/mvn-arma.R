
#' Multivariate normal distribution
#'
#' @param x         \eqn{n \times m}{n*m} numeric matrix.
#' @param n         number of observations. If length(n) > 1,
#'                  the length is taken to be the number required.
#' @param mu        vector of means of length \eqn{m}.
#' @param sigma     \eqn{m \times m}{m*m} covariance matrix.
#' @param log.prob  if \code{TRUE}, probabilities p are given as log(p).
#'
#'
#' @details
#'
#' Multivariate normal probability density function is
#'
#' \deqn{
#' f(x_1,\dots,x_n) = \frac{1}{ (2\pi)^{m/2} \sqrt{\mathrm{det}(\Sigma)} } \,
#' \exp\left\{ -\frac{1}{2}(\mathbf{x}-\boldsymbol{\mu})' \Sigma^{-1} (\mathbf{x}-\boldsymbol{\mu}) \right\}
#' }{
#' f(x) = 1/[(2\pi)^(m/2) * sqrt(det(\Sigma))] * exp( -1/2 (x-\mu)' inv(\Sigma) (x-\mu) )
#' }
#'
#' Random generation from this distribution is possible by taking
#'
#' \deqn{
#' x = A' z + \mu
#' }{
#' x = A' z + \mu
#' }
#'
#' where \eqn{z} is a vector of \eqn{m} i.i.d. standard normal deviates,
#' \eqn{\mu} is a vector of means and \eqn{A} is a \eqn{m \times m}{m*m}
#' matrix such that \eqn{A'A=\Sigma}{A'A=\Sigma} (\eqn{A} is a Cholesky
#' factor of \eqn{\Sigma}).
#'
#' RcppArmadillo implementation of probability density function and random generation
#' is based on examples from the papers by Nino Hardt and Dicko Ahmadou.
#'
#'
#' @references
#' Gentle, J.E. (2006). Random number generation and Monte Carlo methods. Springer.
#'
#' @references
#' Nino Hardt, Dicko Ahmadou (Jul 13, 2013).
#' Faster Multivariate Normal densities with RcppArmadillo and OpenMP.
#' \url{http://gallery.rcpp.org/articles/dmvnorm_arma/}
#'
#' @references
#' Ahmadou Dicko (Mar 12, 2013).
#' Generating a multivariate gaussian distribution using RcppArmadillo.
#' \url{http://gallery.rcpp.org/articles/simulate-multivariate-normal/}
#'
#'
#' @examples
#'
#' partmp <- par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
#'
#' for (r in c(0, -0.25, 0.5, 0.8)) {
#'
#'   mu <- c(0, 0)
#'   Sigma <- matrix(c(1, r, r, 1), ncol = 2)
#'
#'   pal <- colorRampPalette(c("chartreuse4", "yellow", "orange", "brown"))
#'   x <- rmvn(5000, mu, Sigma)
#'   col <- pal(10)[cut(dmvn(x, mu, Sigma), breaks = 10)]
#'
#'   plot(x, col = col, pch = 16, axes = FALSE, main = paste0("correlation = ", r))
#'   axis(1); axis(2)
#'
#'   grid <- seq(-4, 4, by = 0.1)
#'   contour(grid, grid, z = outer(grid, grid, function(x, y) {
#'                              dmvn(cbind(x, y), mu, Sigma)) }, add = TRUE)
#'
#' }
#'
#' par(partmp)
#'
#'
#' @importFrom stats dnorm rnorm
#' @export

dmvn <- function(x, mu, sigma, log.prob = FALSE) {
  drop(cpp_dmvn(as.matrix(x), mu, as.matrix(sigma), log.prob))
}

#' @rdname dmvn
#' @export

rmvn <- function(n, mu, sigma) {
  if (length(n) > 1L) n <- length(n)
  cpp_rmvn(n, mu, as.matrix(sigma))
}

