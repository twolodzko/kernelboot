
#' Random generation from univariate kernel density
#'
#' @param y         numeric vector of length \eqn{n}.
#' @param n         number of observations. If length(n) > 1,
#'                  the length is taken to be the number required.
#' @param bw        the smoothing bandwidth to be used. The kernels are scaled
#'                  such that this is the standard deviation of the smoothing
#'                  kernel (see \code{\link[stats]{density}} for details).
#' @param weights   numeric vector of length \eqn{n}; must be non-negative.
#' @param kernel    a character string giving the smoothing kernel to be used.
#'                  This must partially match one of "gaussian", "rectangular",
#'                  "triangular", "epanechnikov", "biweight", "cosine"
#'                  or "optcosine", with default "gaussian", and may be abbreviated.
#' @param adjust    scalar; the bandwidth used is actually \code{adjust*bw}.
#'                  This makes it easy to specify values like 'half the default'
#'                  bandwidth.
#' @param shrinked  if \code{TRUE} random generation algorithm preserves
#'                  mean and variance of the original sample.
#'
#' @details
#'
#' Univariate kernel density estimator is defined as
#'
#' \deqn{
#' \hat{f_h}(x) = \sum_{i=1}^n w_i \, K_h\left(\frac{x-y_i}{h}\right)
#' }{
#' f(x) = sum[i](w[i] * Kh((x-y[i])/h))
#' }
#'
#' where \eqn{w} is a vector of weights such that \eqn{\sum_i w_i = 1}{sum(w) = 1}
#' (by default \eqn{w_i=1/n}{w[i]=1/n} for all \eqn{i}), \eqn{K_h = K(x/h)/h}{Kh = K(x/h)/h} is
#' kernel \eqn{K} parametrized by bandwidth \eqn{h} and \eqn{y} is a vector of
#' data points used for estimating the kernel density.
#'
#' For functions estimating kernel densities please check \pkg{KernSmooth},
#' \pkg{ks}, or other packages reviewed by Deng and Wickham (2011).
#'
#' For random generation the algorithm described in \code{\link{kernelboot}} is used.
#'
#'
#' @references
#' Deng, H. and Wickham, H. (2011). Density estimation in R.
#' \url{http://vita.had.co.nz/papers/density-estimation.pdf}
#'
#'
#' @examples
#'
#' hist(ruvk(1e5, mtcars$mpg), 100, freq = FALSE, xlim = c(5, 40))
#' lines(density(mtcars$mpg, bw = bw.nrd0(mtcars$mpg)), col = "red")
#'
#' hist(ruvk(1e5, mtcars$mpg, shrinked = TRUE), 100, freq = FALSE, xlim = c(5, 40))
#' lines(density(mtcars$mpg, bw = bw.nrd0(mtcars$mpg)), col = "red")
#'
#' # Comparison of different univariate kernels under standard parametrization
#'
#' kernels <- c("gaussian", "epanechnikov", "rectangular", "triangular",
#'              "biweight", "cosine", "optcosine")
#'
#' partmp <- par(mfrow = c(2, 4), mar = c(3, 3, 3, 3))
#' for (k in kernels) {
#'   hist(ruvk(1e5, 0, 1, kernel = k), 25, freq = FALSE, main = k)
#'   lines(density(0, 1, kernel = k), col = "red")
#' }
#' par(partmp)
#'
#' @seealso \code{\link{kernelboot}}
#'
#' @importFrom stats density bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv
#' @export

ruvk <- function(n, y, bw = bw.nrd0(y),
                 kernel = c("gaussian", "epanechnikov", "rectangular",
                            "triangular", "biweight", "cosine", "optcosine"),
                 weights = NULL, adjust = 1, shrinked = FALSE) {

  kernel <- match.arg(kernel)
  if (length(n) > 1L) n <- length(n)
  if (is.null(weights)) weights <- 1
  if (length(bw) > 1L) {
    bw <- bw[1L]
    message("bw has length > 1 and only the first element will be used")
  }
  bw <- bw * adjust[1L]

  drop(cpp_ruvk(n, y, bw, weights, kernel, shrinked)$samples)
}

