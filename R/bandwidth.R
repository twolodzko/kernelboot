

#' Bandwidth selector for multivariate kernel density estimation
#'
#' Rule of thumb bandwidth selectors for Gaussian kernels as described by
#' Scott (1992) and Silverman (1986).
#'
#' @param x      numeric matrix or data.frame.
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
#' is the empirical covariance matrix. The bandwidth is returned as a \emph{covariance matrix},
#' so to use it for a product kernel, take square root of it's diagonal: \code{sqrt(diag(H))}.
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
#' Chacon J.E., Duong, T. and Wand, M.P. (2011). Asymptotics for general multivariate kernel density
#' derivative estimators. Statistica Sinica, 21, 807-840.
#'
#' @references
#' Epanechnikov, V.A. (1969). Non-parametric estimation of a multivariate orobability density.
#' Theory of Probability & Its Applications, 14(1): 153-158.
#'
#'
#' @seealso \code{\link[stats]{bandwidth}}
#'
#' @importFrom stats var
#'
#' @export

bw.silv <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  if (!all(is_numeric(x)))
    stop("all columns need to be numeric-valued")
  m <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  (4/(n*(m + 2)))^(2/(m + 4)) * S
}


#' @rdname bw.silv
#' @export

bw.scott <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  if (!all(is_numeric(x)))
    stop("all columns need to be numeric-valued")
  m <- ncol(x)
  n <- nrow(x)
  S <- var(x)
  n^(-2/(m + 4)) * S
}
