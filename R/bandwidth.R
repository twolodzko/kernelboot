

#' Bandwidth Selectors for Multivariate Kernel Density Estimation
#'
#' Bandwidth selectors for Gaussian kernels
#'
#' @param x numeric vector.
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
  sigma <- apply(x, 2, sd)
  d <- ncol(x)
  n <- nrow(x)
  (4/(d+2))^(1/(d+4)) * n^(-1/(d+4)) * sigma
}


#' @rdname bandwidth
#' @export

bw.scott <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x)))
    stop("this method works only for matrix, or data.frame objects")
  sigma <- apply(x, 2, sd)
  d <- ncol(x)
  n <- nrow(x)
  n^(-1/(d+4)) * sigma
}
