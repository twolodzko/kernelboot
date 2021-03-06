
#' Random generation from univariate kernel density
#'
#' @param y         numeric vector.
#' @param n         number of observations. If \code{length(n) > 1},
#'                  the length is taken to be the number required.
#' @param bw        the smoothing bandwidth to be used. The kernels are scaled
#'                  such that this is the standard deviation of the smoothing
#'                  kernel (see \code{\link[stats]{density}} for details).
#' @param weights   numeric vector of length equal to \code{length(y)}; must be
#'                  non-negative.
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
#' \hat{f_h}(x) = \sum_{i=1}^n w_i \, K_h(x-y_i)
#' }{
#' f(x) = sum[i](w[i] * Kh(x-y[i]))
#' }
#'
#' where \eqn{w} is a vector of weights such that all \eqn{w_i \ge 0}{w[i] \ge 0}
#' and \eqn{\sum_i w_i = 1}{sum(w) = 1} (by default uniform \eqn{1/n} weights are used),
#' \eqn{K_h = K(x/h)/h}{Kh = K(x/h)/h} is kernel \eqn{K} parametrized by bandwidth
#' \eqn{h} and \eqn{y} is a vector of data points used for estimating the kernel density.
#'
#' For estimating kernel densities use the \code{\link[stats]{density}} function.
#'
#' The random generation algorithm is described in the documentation of
#' \code{\link{kernelboot}} function.
#'
#'
#' @references
#' Deng, H. and Wickham, H. (2011). Density estimation in R.
#' \url{http://vita.had.co.nz/papers/density-estimation.pdf}
#'
#'
#' @examples
#'
#' # ruvk() produces samples from kernel densities as estimated using
#' # density() function from base R
#'
#' hist(ruvk(1e5, mtcars$mpg), 100, freq = FALSE, xlim = c(5, 40))
#' lines(density(mtcars$mpg, bw = bw.nrd0(mtcars$mpg)), col = "red")
#'
#' # when using 'shrinked = TRUE', the samples differ from density() estimates
#' # since they are shrinked to have the same variance as the underlying data
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
#' @seealso \code{\link{kernelboot}}, \code{\link[stats]{density}}
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
  as.vector(cpp_ruvk(n, y, bw, weights, kernel, shrinked))

}
