
#' Multivariate product kernel density
#'
#' @param x            numeric vector; kernel density is evaluated on those values.
#' @param y            numeric vector; kernel density is estimated using those values.
#' @param n            number of observations. If length(n) > 1,
#'                     the length is taken to be the number required.
#' @param bw           numeric matrix
#' @param weights      numeric vector
#' @param adjust       scalar; the bandwidth used is actually \code{adjust*bw}.
#'                     This makes it easy to specify values like 'half the default'
#'                     bandwidth.
#' @param kernel       a character string giving the smoothing kernel to be used.
#'                     This must partially match one of "gaussian", "rectangular",
#'                     "triangular", "epanechnikov", "biweight", "cosine" or
#'                     "optcosine", with default "gaussian", and may be abbreviated
#'                     to a unique prefix (single letter).
#' @param log.prob     logical; if \code{TRUE}, probabilities p are given as log(p).
#'
#'
#' @details
#'
#' Multivariate product kernel density estimator is defined as a product of univariate kernels
#'
#' \deqn{
#' \hat{f_H}(x_1,\dots,x_n) = \sum_{i=1}^n w_i \prod_{j=1}^m
#' K_{h_j} \left( \frac{x_i - y_{ij}}{h_j} \right)
#' }{
#' f(x) = sum[i](w[i] * prod[j]( Khj((x[i]-y[i,j])/h[j]) ))
#' }
#'
#' where \eqn{w} is a vector of weights such that \eqn{\sum_i w_i = 1}{sum(w) = 1}
#' (by default \eqn{w_i=1/n}{w[i]=1/n} for all \eqn{i}), \eqn{K_{h_j}}{Khj} is
#' kernel \eqn{K} parametrized by bandwidth \eqn{h_j} and \eqn{\boldsymbol{y}}{y}
#' is a matrix of data points used for estimating the kernel density.
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
#' @seealso \code{\link{kernelboot}}
#'
#'
#' @export

dmvpkd <- function(x, y, bw = sqrt(diag(bw.silv(y))), weights = NULL,
                   kernel = c("gaussian", "epanechnikov", "rectangular",
                              "triangular", "biweight", "triweight",
                              "cosine", "optcosine"),
                   adjust = 1, log.prob = FALSE) {
  kernel <- match.arg(kernel)
  if (is.null(weights)) weights <- 1
  bw <- bw * adjust[1L]
  drop(cpp_dmvpkd(x, y, bw, weights, kernel, log.prob)$density)
}

#' @rdname dmvpkd
#' @export

rmvpkd <- function(n, y, bw = sqrt(diag(bw.silv(y))), weights = NULL,
                   kernel = c("gaussian", "epanechnikov", "rectangular",
                              "triangular", "biweight", "triweight",
                              "cosine", "optcosine"), adjust = 1) {
  kernel <- match.arg(kernel)
  if (length(n) > 1) n <- length(n)
  if (is.null(weights)) weights <- 1
  bw <- bw * adjust[1L]
  cpp_rmvpkd(n, y, bw, weights, kernel)$samples
}

