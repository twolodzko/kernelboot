
#' Multivariate kernel density
#'
#' @param y,x          numeric matrix.
#' @param n            number of observations. If length(n) > 1,
#'                     the length is taken to be the number required.
#' @param bw           numeric matrix
#' @param weights      numeric vector
#' @param preserve.var logical;
#' @param log.prob     logical; if TRUE, probabilities p are given as log(p).
#'
#' @export

dmvkd <- function(y, x, bw = bw.scott(x), weights = rep(1, nrow(x)), log.prob = FALSE) {
  drop(cpp_dmvkde(y, x, bw, weights, log.prob, FALSE)$density)
}

#' @rdname dmvkd
#' @export

rmvkd <- function(n, x, bw = bw.scott(x), weights = rep(1, nrow(x)), preserve.var = FALSE) {
  if (length(n) > 1) n <- length(n)
  cpp_rmvkde(n, x, bw, weights, preserve.var)$sample
}

