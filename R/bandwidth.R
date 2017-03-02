

#' Bandwidth selector for multivariate kernel density estimation
#'
#' Rule of thumb bandwidth selectors for Gaussian kernels as described by
#' Scott (1992) and Silverman (1986).
#'
#' @param x      numeric matrix or data.frame.
#' @param diag   if \code{TRUE} returns the diagonal of bandwidth matrix.
#'
#'
#' @details
#'
#' Scott's (1992) rule is defined as
#'
#' \deqn{
#' H = n^{-2/(m+4)} \hat\Sigma
#' }{
#' H = n^(-2/(m+4)) * S
#' }
#'
#' Silverman's (1986; see Chacon, Duong and Wand, 2011) rule is defined as
#'
#' \deqn{
#' H = \left(\frac{4}{n(m+2)}\right)^{2/(m+4)} \hat\Sigma
#' }{
#' H = (4/(n*(m+2)))^(2/(m+4)) * S
#' }
#'
#' where \eqn{m} is number of variables, \eqn{n} is sample size, \eqn{\hat\Sigma}{S}
#' is the empirical covariance matrix.
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
#' @references
#' Chacon J.E., Duong, T. and Wand, M.P. (2011). Asymptotics for general multivariate kernel density
#' derivative estimators. Statistica Sinica, 21, 807-840.
#'
#'
#' @seealso \code{\link[stats]{bandwidth}}
#'
#' @importFrom stats var
#'
#' @name bandwidth
#' @aliases bw.silv
#' @export

bw.silv <- function(x, diag = FALSE) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  m <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  if (diag) S <- S * diag(m)
  (4/(n*(m + 2)))^(2/(m + 4)) * S
}


#' @rdname bandwidth
#' @export

bw.scott <- function(x, diag = FALSE) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  m <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  if (diag) S <- S * diag(m)
  n^(-2/(m + 4)) * S
}
