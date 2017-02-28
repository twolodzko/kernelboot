
# check if different objects are numeric
# for data.frames and matrix objects check the individual columns

is_numeric <- function(x) UseMethod("is_numeric")

is_numeric.default <- function(x) {
  is.numeric(x)
}

is_numeric.matrix <- function(x) {
  structure(rep(is.numeric(x), ncol(x)), names = colnames(x))
}

is_numeric.data.frame <- function(x) {
  unlist(lapply(x, is.numeric))
}

# matrix rank

mtxrank <- function(x) {
  sum(svd(x)$d > 1e-12)
}

# check for square matrix

is.square <- function(x) {
  NCOL(x) == NROW(x)
}

# check for diagonal matrix

is.diag <- function(x) {
  diag(x) <- 0
  isTRUE(all.equal(x, matrix(0, nrow(x), ncol(x)),
                   check.attributes = FALSE,
                   check.names = FALSE))
}
