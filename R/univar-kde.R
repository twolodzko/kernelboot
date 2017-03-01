
#' Univariate kernel density
#'
#' @param y,x          numeric vector.
#' @param n            number of observations. If length(n) > 1,
#'                     the length is taken to be the number required.
#' @param bw           numeric value.
#' @param weights      numeric vector.
#' @param kernel       a character string giving the smoothing kernel to be used.
#'                     This must partially match one of "gaussian", "rectangular",
#'                     "triangular", "epanechnikov", "biweight", "cosine" or
#'                     "optcosine", with default "gaussian", and may be abbreviated
#'                     to a unique prefix (single letter).
#' @param preserve.var logical;
#' @param log.prob     logical; if TRUE, probabilities p are given as log(p).
#'
#' @export

duvkd <- function(y, x, bw = bw.nrd0(x), weights = rep(1, length(x)),
                   kernel = c("gaussian", "epanechnikov", "rectangular",
                              "triangular", "biweight", "triweight",
                              "cosine", "optcosine"), log.prob = FALSE) {
  kernel <- match.arg(kernel)
  drop(cpp_duvkde(y, x, bw, weights, kernel, log.prob, FALSE)$density)
}

#' @rdname duvkd
#' @export

ruvkd <- function(n, x, bw = bw.nrd0(x), weights = rep(1, length(x)),
                   kernel = c("gaussian", "epanechnikov", "rectangular",
                              "triangular", "biweight", "triweight",
                              "cosine", "optcosine"), preserve.var = FALSE) {
  kernel <- match.arg(kernel)
  if (length(n) > 1) n <- length(n)
  drop(cpp_ruvkde(n, x, bw, weights, kernel, preserve.var)$sample)
}

