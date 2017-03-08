
#' Multivariate product kernel
#'
#' @param x         \eqn{k \times m}{k*m} numeric matrix; kernel density
#'                  is evaluated on those values.
#' @param y         \eqn{n \times m}{n*m} numeric matrix; kernel density
#'                  is estimated using those values.
#' @param n         number of observations. If length(n) > 1,
#'                  the length is taken to be the number required.
#' @param bw        numeric vector of length \eqn{m}; the smoothing bandwidth
#'                  to be used. The kernels are scaled such that this is the
#'                  standard deviation of the smoothing kernel (see
#'                  \code{\link[stats]{density}} for details). If provided as
#'                  a single value, the same bandwidth is used for each variable.
#' @param weights   numeric vector of length \eqn{n}; must be non-negative.
#' @param adjust    scalar; the bandwidth used is actually \code{adjust*bw}.
#'                  This makes it easy to specify values like 'half the default'
#'                  bandwidth.
#' @param kernel    a character string giving the smoothing kernel to be used.
#'                  This must partially match one of "gaussian", "rectangular",
#'                  "triangular", "epanechnikov", "biweight", "triweight", "cosine"
#'                  or "optcosine", with default "gaussian", and may be abbreviated.
#' @param shrinked  if \code{TRUE} random generation algorithm preserves mean and
#'                  variances of the individual variables (see \code{\link{ruvk}}).
#'                  Shrinking is applied to each of the variables individually.
#' @param log.prob  if \code{TRUE}, probabilities p are given as log(p).
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
#' Random generation from product kernel is done by drawing with replacement
#' rows of \code{y}, and then adding random noise from univariate kernel \eqn{K}
#' (see \code{\link{duvk}}), parametrized by corresponding bandwidth parameter
#' \eqn{h}, to the sampled values.
#'
#' When \code{shrinked = TRUE}, product kernel density is a product of shrinked univariate
#' kernels (see \code{\link{duvk}}).
#'
#'
#' @references
#' Silverman, B.W. (1986). Density estimation for statistics and data analysis. Chapman and Hall/CRC.
#'
#' @references
#' Wand, M.P. and Jones, M.C. (1995). Kernel smoothing. Chapman and Hall/CRC.
#'
#' @references
#' Scott, D.W. (1992). Multivariate density estimation: theory, practice,
#' and visualization. John Wiley & Sons.
#'
#' @references
#' Epanechnikov, V.A. (1969). Non-parametric estimation of a multivariate orobability density.
#' Theory of Probability & Its Applications, 14(1): 153-158.
#'
#'
#' @examples
#'
#' dat <- as.matrix(mtcars[, c("mpg", "disp")])
#' pal <- colorRampPalette(c("chartreuse4", "yellow", "orange", "brown"))
#'
#' gridx <- seq(0, 45, length.out = 200)
#' gridy <- seq(-200, 800, length.out = 200)
#'
#' partmp <- par(mfrow = c(1, 2), mar = c(3,3,3,3))
#'
#' samp1 <- rmvpk(5000, dat, shrink = FALSE)
#' col1 <- pal(10)[cut(dmvpk(samp1, dat, shrink = FALSE), breaks = 10)]
#'
#' plot(samp1, col = col1, pch = 16, axes = FALSE, xlim = c(0, 45), ylim = c(-200, 800))
#' points(dat, pch = 2, lwd = 2)
#' fx1 <- outer(gridx, gridy, function(x, y) dmvpk(cbind(x, y), dat, shrink = FALSE))
#' contour(gridx, gridy, z = fx1, add = TRUE)
#' axis(1)
#' axis(2)
#' title("Product kernel", cex.sub = 0.5)
#' legend("topright", pch = c(2, 16), col = c("black", "chartreuse4"),
#'        legend = c("actual data", "bootstrap samples"), bty = "n", cex = 0.8 )
#'
#'
#' samp2 <- rmvpk(5000, dat, shrink = TRUE)
#' col2 <- pal(10)[cut(dmvpk(samp2, dat, shrink = TRUE), breaks = 10)]
#'
#' plot(samp2, col = col2, pch = 16, axes = FALSE, xlim = c(0, 45), ylim = c(-200, 800))
#' points(dat, pch = 2, lwd = 2)
#' fx2 <- outer(gridx, gridy, function(x, y) dmvpk(cbind(x, y), dat, shrink = TRUE))
#' contour(gridx, gridy, z = fx2, add = TRUE)
#' axis(1)
#' axis(2)
#' title("Product kernel (shrinked)", cex.sub = 0.5)
#' legend("topright", pch = c(2, 16), col = c("black", "chartreuse4"),
#'        legend = c("actual data", "bootstrap samples"), bty = "n", cex = 0.8 )
#'
#' par(partmp)
#'
#' cov(dat)
#' cov(samp1)
#' cov(samp2)
#'
#'
#' @seealso \code{\link{kernelboot}}, \code{\link{duvk}}, \code{\link{dmvk}}
#'
#'
#' @export

dmvpk <- function(x, y, bw = sqrt(diag(bw.silv(y))),
                  kernel = c("gaussian", "epanechnikov", "rectangular",
                             "triangular", "biweight", "triweight",
                             "cosine", "optcosine"),
                  weights = NULL, adjust = 1, shrinked = FALSE,
                  log.prob = FALSE) {

  kernel <- match.arg(kernel)

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1L)
  } else {
    x <- as.matrix(x)
  }

  if (is.vector(y)) {
    y <- matrix(y, nrow = 1L)
  } else {
    y <- as.matrix(y)
  }

  if (is.null(weights)) weights <- 1

  if (is.matrix(bw) || is.data.frame(bw)) {
    if (!is.square(bw))
      stop("bw is not a square matrix")
    bw <- diag(bw)
  }
  bw <- bw * adjust[1L]
  if (length(bw) == 1L)
    bw <- rep(bw, ncol(y))

  drop(cpp_dmvpk(x, y, bw, weights, kernel, shrinked, log.prob)$density)
}


#' @rdname dmvpk
#' @export

rmvpk <- function(n, y, bw = sqrt(diag(bw.silv(y))),
                  kernel = c("gaussian", "epanechnikov", "rectangular",
                             "triangular", "biweight", "triweight",
                             "cosine", "optcosine"),
                  weights = NULL, adjust = 1, shrinked = FALSE) {

  kernel <- match.arg(kernel)
  if (length(n) > 1L) n <- length(n)

  if (is.vector(y)) {
    y <- matrix(y, nrow = 1L)
  } else {
    y <- as.matrix(y)
  }

  if (is.null(weights)) weights <- 1

  if (is.matrix(bw) || is.data.frame(bw)) {
    if (!is.square(bw))
      stop("bw is not a square matrix")
    bw <- diag(bw)
  }
  bw <- bw * adjust[1L]
  if (length(bw) == 1L)
    bw <- rep(bw, ncol(y))

  out <- cpp_rmvpk(n, y, bw, weights, kernel, shrinked)$samples
  colnames(out) <- colnames(y)
  out
}

