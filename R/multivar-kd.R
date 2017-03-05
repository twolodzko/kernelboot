
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
#'                     while the univariate (\code{\link{duvk}}) and product
#'                     kernels (\code{\link{dmvpk}}) are parametrized in terms
#'                     of standard deviations as in \code{\link{density}} function.
#'                     If provided as a single value, the same bandwidth is used
#'                     for each variable.
#' @param weights      numeric vector of length \eqn{n}; must be non-negative.
#' @param adjust       scalar; the bandwidth used is actually \code{adjust*bw}.
#'                     This makes it easy to specify values like 'half the default'
#'                     bandwidth.
#' @param log.prob     logical; if \code{TRUE}, probabilities p are given as log(p).
#'
#'
#' @details
#'
#' Multivariate Gaussian kernel density estimator is defined as
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
#' Silverman, B.W. (1986). Density estimation for statistics and data analysis. Chapman and Hall/CRC.
#'
#' @references
#' Wand, M.P. and Jones, M.C. (1995). Kernel Smoothing. Chapman and Hall/CRC.
#'
#' @references
#' Scott, D.W. (1992). Multivariate density estimation: theory, practice,
#' and visualization. John Wiley & Sons.
#'
#'
#' @examples
#'
#' # Comparison of multivariate and product kernels
#'
#' dat <- as.matrix(mtcars[, c("mpg", "disp")])
#' pal <- colorRampPalette(c("chartreuse4", "yellow", "orange", "brown"))
#'
#' partmp <- par(mfrow = c(1, 2), mar = c(3,3,3,3))
#'
#' samp1 <- rmvpk(5000, dat)
#' col1 <- pal(10)[cut(dmvpk(samp1, dat), breaks = 10)]
#'
#' plot(samp1, col = col1, pch = 16, axes = FALSE)
#' points(dat, pch = 2, lwd = 2)
#' axis(1); axis(2)
#' title("Product kernel", cex.sub = 0.5)
#' legend("topright", pch = c(2, 16), col = c("black", "chartreuse4"),
#'        legend = c("actual data", "bootstrap samples"), bty = "n", cex = 0.8 )
#'
#'
#' samp2 <- rmvk(5000, dat)
#' col2 <- pal(10)[cut(dmvk(samp2, dat), breaks = 10)]
#'
#' plot(samp2, col = col2, pch = 16, axes = FALSE)
#' points(dat, pch = 2, lwd = 2)
#' axis(1); axis(2)
#' title("Multivariate Gaussian kernel", cex.sub = 0.5)
#' legend("topright", pch = c(2, 16), col = c("black", "chartreuse4"),
#'        legend = c("actual data", "bootstrap samples"), bty = "n", cex = 0.8 )
#'
#' par(partmp)
#'
#'
#' @seealso \code{\link{kernelboot}}, \code{\link{dmvn}}
#'
#' @importFrom grDevices colorRampPalette
#' @export

dmvk <- function(x, y, bw = bw.silv(y), weights = NULL,
                  adjust = 1, log.prob = FALSE) {

  if (is.null(weights)) weights <- 1
  bw <- bw * adjust[1L]
  if (is.vector(bw)) {
    if (length(bw) == 1L)
      bw <- diag(ncol(y)) * bw
    else
      bw <- diag(bw)
  }

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  } else {
    x <- as.matrix(x)
  }

  if (is.vector(y)) {
    y <- matrix(y, nrow = 1)
  } else {
    y <- as.matrix(y)
  }

  drop(cpp_dmvk(x, y, bw, weights, log.prob, FALSE)$density)
}


#' @rdname dmvk
#' @export

rmvk <- function(n, y, bw = bw.silv(y), weights = NULL,
                  adjust = 1) {

  if (length(n) > 1L) n <- length(n)
  if (is.null(weights)) weights <- 1
  bw <- bw * adjust[1L]
  if (is.vector(bw)) {
    if (length(bw) == 1L)
      bw <- diag(ncol(y)) * bw
    else
      bw <- diag(bw)
  }

  if (is.vector(y)) {
    y <- matrix(y, nrow = 1)
  } else {
    y <- as.matrix(y)
  }

  out <- cpp_rmvk(n, y, bw, weights, FALSE)$samples
  colnames(out) <- colnames(y)
  out
}

