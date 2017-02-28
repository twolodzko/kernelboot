
#' Multivariate normal distribution
#'
#' @param x        numeric matrix.
#' @param n        number of observations. If length(n) > 1,
#'                 the length is taken to be the number required.
#' @param mu       numeric vector.
#' @param sigma    numeric matrix.
#' @param log.prob logical; if TRUE, probabilities p are given as log(p).
#'
#' @author
#' Nino Hardt, Dicko Ahmadou, Tymoteusz Wolodzko
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
#' @export

dmvn <- function(x, mu, sigma, log.prob = FALSE) {
  cpp_dmvn(x, mu, sigma, log.prob)
}

#' @rdname dmvn
#' @export

rmvn <- function(n, mu, sigma) {
  if (length(n) > 1) n <- length(n)
  cpp_rmvn(n, mu, sigma)
}

