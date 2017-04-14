
#' Random generation from multivariate Gaussian kernel density
#'
#' @param y         numeric matrix.
#' @param n         number of observations. If \code{length(n) > 1},
#'                  the length is taken to be the number required.
#' @param bw        numeric matrix with number of rows and columns equal to
#'                  \code{ncol(y)}; the smoothing bandwidth to be used. This is the
#'                  \emph{covariance matrix} of the smoothing kernel. If provided as
#'                  a single value, the same bandwidth is used for each variable.
#'                  If provided as a single value, or as a vector, variables are
#'                  considered as uncorrelated.
#' @param weights   numeric vector of length equal to \code{nrow(y)}; must be non-negative.
#' @param adjust    scalar; the bandwidth used is actually \code{adjust*bw}.
#'                  This makes it easy to specify values like 'half the default'
#'                  bandwidth.
#'
#'
#' @details
#'
#' Multivariate kernel density estimator with multivariate Gaussian (normal) kernels
#' \eqn{K_H}{KH} is defined as
#'
#' \deqn{
#' \hat{f_H}(\mathbf{x}) = \sum_{i=1}^n w_i \, K_H \left( \mathbf{x}-\boldsymbol{y}_i \right)
#' }{
#' f(x) = sum[i](w[i] * KH(x-y[i]))
#' }
#'
#' where \eqn{w} is a vector of weights such that \eqn{\sum_i w_i = 1}{sum(w) = 1}
#' (by default uniform weights are used), \eqn{K_H}{KH} is kernel \eqn{K} parametrized by
#' bandwidth matrix \eqn{H} and \eqn{\boldsymbol{y}}{y} is a matrix of data points used for
#' estimating the kernel density.
#'
#' Random generation from multivariate normal distribution is possible by taking
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
#' factor of \eqn{\Sigma}). In the case of multivariate Gaussian kernel
#' density, \eqn{\mu}, is the \eqn{i}-th row of \eqn{y}, where \eqn{i}
#' is drawn randomly with replacement with probability proportional to
#' \eqn{w_i}{w[i]}, and \eqn{\Sigma} is the bandwidth matrix \eqn{H}.
#'
#' For functions estimating kernel densities please check \pkg{KernSmooth},
#' \pkg{ks}, or other packages reviewed by Deng and Wickham (2011).
#'
#'
#' @references
#' Deng, H. and Wickham, H. (2011). Density estimation in R.
#' \url{http://vita.had.co.nz/papers/density-estimation.pdf}
#'
#' @examples
#'
#' dat <- mtcars[, c("mpg", "disp")]
#'
#' plot(rmvg(5000, dat), col = "#ADD8E640", pch = 16,
#'      xlim = c(0, 45), ylim = c(-200, 800),
#'      main = "Multivariate Gaussian kernel")
#' points(dat, pch = 2, lwd = 2, col = "red")
#'
#' @seealso \code{\link{kernelboot}}
#'
#' @export

rmvg <- function(n, y, bw = bw.silv(y), weights = NULL, adjust = 1) {

  if (length(n) > 1L) n <- length(n)

  if (is.simple.vector(y)) {
    y <- matrix(y, nrow = 1L)
  } else {
    y <- as.matrix(y)
  }

  if (is.matrix(bw) || is.data.frame(bw)) {
    if (!is.square(bw))
      stop("bw is not a square matrix")
    bw <- as.matrix(bw)
  } else if (is.simple.vector(bw)) {
    if (length(bw) == 1L)
      bw <- diag(bw, ncol(y))
    else
      bw <- diag(bw)
  }
  if (ncol(bw) != ncol(y))
    stop("bw has incorrect dimensions")

  if (!is.null(weights) && length(weights) == 1L)
    weights <- NULL

  bw <- bw * adjust[1L]

  if (!all(is.finite(bw)))
    stop("inappropriate values of bw")

  idx <- sample.int(nrow(y), n, replace = TRUE, prob = weights)
  mu <- y[idx, , drop = FALSE]
  Az <- matrix(rnorm(n*ncol(y)), n, ncol(y)) %*% chol(bw)
  return(Az + mu)

}

