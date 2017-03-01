

#' Bandwidth Selector for Multivariate Kernel Density Estimation
#'
#' Rule of thumb bandwidth selectors for Gaussian kernels as described by
#' Scott (1992) and Silverman (1986).
#'
#' @param x numeric matrix or data.frame.
#'
#' @details
#'
#' Scott's (1992) rule is defined as
#'
#' \deqn{
#' H = \left(n^{-1/(k+4)}\right)^2 \mathrm{diag}(S)
#' }{
#' H = [n^(-1/(k+4))]^2 * diag(S)
#' }
#'
#' Silverman's (1986) rule is defined as
#'
#' \deqn{
#' H = \left(\left(\frac{4}{(k+2)n}\right)^{1/(k+4)}\right)^2 \mathrm{diag}(S)
#' }{
#' H = [(4/((k+2)*n))^(1/(k+4))]^2 * diag(S)
#' }
#'
#' where \eqn{k} is number of variables, \eqn{n} is sampel size and \eqn{S}
#' is the empirical covariance matrix.
#'
#' @references
#' Silverman, B. W. (1986). Density estimation for statistics and data analysis. Chapman and Hall/CRC.
#'
#' @references
#' Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman and Hall/CRC.
#'
#' @references
#' Scott, D. W. (1992). Multivariate density estimation: theory, practice,
#' and visualization. John Wiley & Sons.
#'
#' @seealso \code{\link[stats]{bandwidth}}
#'
#' @importFrom stats cov
#' @export

bw.scott <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  d <- ncol(x)
  n <- nrow(x)
  S <- cov(x) * diag(d)
  (n^(-1/(d+4)))^2 * S
}


#' @rdname bw.scott
#' @export

bw.silv <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  d <- ncol(x)
  n <- nrow(x)
  S <- cov(x) * diag(d)
  (4/(d+2))^(1/(d+4)) * n^(-1/(d+4))^2 * S
}
