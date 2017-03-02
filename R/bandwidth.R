

#' Bandwidth selector for multivariate kernel density estimation
#'
#' Rule of thumb bandwidth selectors for Gaussian kernels as described by
#' Scott (1992), Silverman (1986) and Chacon, Duong and Wand (2011).
#'
#' @param x      numeric matrix or data.frame.
#' @param diag   if \code{TRUE} returns the diagonal of bandwidth matrix.
#' @param r      derivative order in \code{bw.ns} (zero by default).
#'
#' @details
#'
#' Scott's (1992) rule is defined as
#'
#' \deqn{
#' H = \left(n^{-1/(m+4)}\right)^2 \hat\Sigma
#' }{
#' H = [n^(-1/(m+4))]^2 * S
#' }
#'
#' Silverman's (1986) rule is defined as
#'
#' \deqn{
#' H = \left(\left(\frac{4}{(m+2)n}\right)^{1/(m+4)}\right)^2 \hat\Sigma
#' }{
#' H = [(4/((m+2)*n))^(1/(m+4))]^2 * S
#' }
#'
#' Normal scale bandwidth selector (Chacon, Duong and Wand, 2011) is defined as
#'
#' \deqn{
#' H = \left(\frac{4}{n (m+2r+2)}\right)^{2/(m+2r+4)} \hat\Sigma
#' }{
#' H = (4/(n*(m+2*r+2)))^(2/(m+2*r+4)) * S
#' }
#'
#' where \eqn{m} is number of variables, \eqn{n} is sampel size, \eqn{\hat\Sigma}{S}
#' is the empirical covariance matrix and \eqn{r} is derivative order.
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
#'
#' @seealso \code{\link[stats]{bandwidth}}
#'
#' @importFrom stats var
#'
#' @name bandwidth
#' @aliases bw.ns
#' @export

bw.ns <- function(x, diag = FALSE, r = 0) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  if (is.vector(x)) {
    n <- 1L
    m <- length(x)
  } else {
    n <- nrow(x)
    m <- ncol(x)
  }
  S <- var(x)
  if (diag) S <- S * diag(m)
  (4/(n * (m + 2*r + 2)))^(2/(m + 2*r + 4)) * S
}


#' @rdname bandwidth
#' @export

bw.silv <- function(x, diag = FALSE) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  m <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  if (diag) S <- S * diag(m)
  ((4/(m + 2))^(1/(m + 4)) * n^(-1/(m + 4)))^2 * S
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
  (n^(-1/(m + 4)))^2 * S
}
