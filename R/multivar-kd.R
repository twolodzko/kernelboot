
#' Multivariate kernel density with Gaussian kernel
#'
#' @param x            \eqn{k \times m}{k*m} numeric matrix; kernel density
#'                     is evaluated on those values.
#' @param y            \eqn{n \times m}{n*m} numeric matrix; kernel density
#'                     is estimated using those values.
#' @param n            number of observations. If length(n) > 1,
#'                     the length is taken to be the number required.
#' @param bw           \eqn{m \times m}{m*m} numeric matrix.
#'                     \emph{Notice:} this is a \emph{covariance} matrix of
#'                     multivariate normal distribution (see \code{\link{dmvn}}),
#'                     while the \code{\link{duvkd}} and \code{\link{dmvpkd}}
#'                     functions are parametrized with standard deviations
#'                     as in \code{\link{density}}.
#' @param weights      numeric vector of length \eqn{n}; must be non-negative.
#' @param adjust       scalar; the bandwidth used is actually \code{adjust*bw}.
#'                     This makes it easy to specify values like 'half the default'
#'                     bandwidth.
#' @param log.prob     logical; if \code{TRUE}, probabilities p are given as log(p).
#'
#'
#' @details
#'
#' Multivariate kernel density estimator is defined as
#'
#' \deqn{
#' \hat{f_H}(x_1,\dots,x_n) = \sum_{i=1}^n w_i \, K_H \left( \mathbf{x}-\boldsymbol{y}_i \right)
#' }{
#' f(x) = sum[i](w[i] * KH(x-y[i]))
#' }
#'
#' where \eqn{w} is a vector of weights such that \eqn{\sum_i w_i = 1}{sum(w) = 1}
#' (by default \eqn{w_i=1/n}{w[i]=1/n} for all \eqn{i}), \eqn{K_H}{KH} is
#' kernel \eqn{K} parametrized by bandwidth matrix \eqn{H} and \eqn{\boldsymbol{y}}{y}
#' is a matrix of data points used for estimating the kernel density.
#'
#' This function uses Gaussian (multiariate normal) kernel (see \code{\link{dmvn}}).
#'
#'
#' @references
#' Silverman, B. W. (1986). Density estimation for statistics and data analysis.
#' Chapman and Hall/CRC.
#'
#' @references
#' Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman and Hall/CRC.
#'
#' @references
#' Scott, D. W. (1992). Multivariate density estimation: theory, practice,
#' and visualization. John Wiley & Sons.
#'
#'
#' @seealso \code{\link{kernelboot}}, \code{\link{dmvn}}
#'
#'
#' @export

dmvkd <- function(x, y, bw = bw.silv(y), weights = NULL,
                  adjust = 1, log.prob = FALSE) {
  if (is.null(weights)) weights <- 1
  bw <- bw * adjust[1L]
  drop(cpp_dmvkd(x, y, bw, weights, log.prob, FALSE)$density)
}

#' @rdname dmvkd
#' @export

rmvkd <- function(n, y, bw = bw.silv(y), weights = NULL,
                  adjust = 1) {
  if (length(n) > 1L) n <- length(n)
  if (is.null(weights)) weights <- 1
  bw <- bw * adjust[1L]
  cpp_rmvkd(n, y, bw, weights, FALSE)$samples
}

