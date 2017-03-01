
#' Univariate kernel density
#'
#' @param y,x          numeric vector.
#' @param n            number of observations. If length(n) > 1,
#'                     the length is taken to be the number required.
#' @param bw           numeric value.
#' @param weights      numeric vector.
#' @param preserve.var logical;
#' @param log.prob     logical; if TRUE, probabilities p are given as log(p).
#'
#' @export

duvkde <- function(y, x, bw = bw.nrd0(x), weights = rep(1, length(x)),
                   kernel = "gaussian", log.prob = FALSE) {
  drop(cpp_duvkde(y, x, bw, weights, kernel, log.prob, FALSE)$density)
}

#' @rdname duvkde
#' @export

ruvkde <- function(n, x, bw = bw.nrd0(x), weights = rep(1, length(x)),
                   kernel = "gaussian", preserve.var = FALSE) {
  if (length(n) > 1) n <- length(n)
  drop(cpp_ruvkde(n, x, bw, weights, kernel, preserve.var)$sample)
}

