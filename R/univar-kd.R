
#' Univariate kernel density
#'
#' @param x            numeric vector of length \eqn{k}; kernel density is
#'                     evaluated on those values.
#' @param y            numeric vector of length \eqn{n}; kernel density is
#'                     evaluated on those values.
#' @param n            number of observations. If length(n) > 1,
#'                     the length is taken to be the number required.
#' @param bw           scalar; must be greater than zero.
#' @param weights      numeric vector of length \eqn{n}; must be non-negative.
#' @param kernel       a character string giving the smoothing kernel to be used.
#'                     This must partially match one of "gaussian", "rectangular",
#'                     "triangular", "epanechnikov", "biweight", "cosine" or
#'                     "optcosine", with default "gaussian", and may be abbreviated
#'                     to a unique prefix (single letter).
#' @param adjust       scalar; the bandwidth used is actually \code{adjust*bw}.
#'                     This makes it easy to specify values like 'half the default'
#'                     bandwidth.
#' @param preserve.var logical; if \code{TRUE} random generation algorithm preserves
#'                     mean and variance of the original sample (see
#'                     \code{\link{kernelboot}} for details).
#' @param log.prob     logical; if \code{TRUE}, probabilities p are given as log(p).
#'
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
#' @seealso \code{\link[stats]{density}}, \code{\link{kernelboot}}
#'
#'
#' @examples
#'
#' hist(ruvkd(1e5, mtcars$mpg), 100, freq = FALSE)
#' curve(duvkd(x, mtcars$mpg), from = 0, to = 50, col = "red", add = TRUE)
#'
#' hist(ruvkd(1e5, mtcars$mpg, preserve.var = TRUE), 100, freq = FALSE)
#' curve(duvkd(x, mtcars$mpg), from = 0, to = 50, col = "red", add = TRUE)
#'
#'
#' @importFrom stats bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv
#' @export

duvkd <- function(x, y, bw = bw.nrd0(y), weights = NULL,
                  kernel = c("gaussian", "epanechnikov", "rectangular",
                              "triangular", "biweight", "triweight",
                              "cosine", "optcosine"),
                  adjust = 1, log.prob = FALSE) {
  kernel <- match.arg(kernel)
  if (is.null(weights)) weights <- 1
  bw <- bw * adjust[1L]
  drop(cpp_duvkd(x, y, bw, weights, kernel, log.prob)$density)
}

#' @rdname duvkd
#' @export

ruvkd <- function(n, y, bw = bw.nrd0(y), weights = NULL,
                  kernel = c("gaussian", "epanechnikov", "rectangular",
                             "triangular", "biweight", "triweight",
                             "cosine", "optcosine"),
                  adjust = 1, preserve.var = FALSE) {
  kernel <- match.arg(kernel)
  if (length(n) > 1L) n <- length(n)
  if (is.null(weights)) weights <- 1
  bw <- bw * adjust[1L]
  drop(cpp_ruvkd(n, y, bw, weights, kernel, preserve.var)$samples)
}

