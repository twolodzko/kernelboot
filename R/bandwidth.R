

#' Bandwidth Selectors for Multivariate Kernel Density Estimation
#'
#' Rule of thumb bandwidth selectors for Gaussian kernels as described by
#' Silverman (1986) and Scott (1992).
#'
#' @param x numeric vector.
#'
#' @details
#'
#' Scott's rule is defined as
#'
#' \deqn{
#' h_i = n^{-1/(d+4)} \sigma_i
#' }{
#' h[i] = n^(-1/(d+4)) * \sigma[i]
#' }
#'
#' Silverman's rule is defined as
#'
#' \deqn{
#' h_i = \left(\frac{4}{d+2}\right)^{1/(d+4)} n^{-1/(d+4)} \sigma_i
#' }{
#' h[i] = (4/(d+2))^(1/(d+4)) * n^(-1/(d+4)) * \sigma[i]
#' }
#'
#' where \eqn{d} is number of variables and \eqn{n} is sampel size.
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
#' @seealso \code{\link[stats]{bw.nrd}}
#'
#' @name bandwidth
#' @aliases bw.silv86
#' @aliases bw.scott
#' @export

bw.silv86 <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  S <- diag(cov(x))
  d <- ncol(x)
  n <- nrow(x)
  (4/(d+2))^(1/(d+4)) * n^(-1/(d+4)) * S
}


#' @rdname bandwidth
#' @export

bw.scott <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  S <- diag(cov(x))
  d <- ncol(x)
  n <- nrow(x)
  n^(-1/(d+4)) * S
}
