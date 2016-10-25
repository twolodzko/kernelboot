

is.square <- function(x) {
  NCOL(x) == NROW(x)
}


is.diag <- function(x) {
  diag(x) <- 0
  isTRUE(all.equal(x, matrix(0, nrow(x), ncol(x)),
                   check.attributes = FALSE,
                   check.names = FALSE))
}
