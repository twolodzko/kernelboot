
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

# check for square matrix

is.square <- function(x) {
  NCOL(x) == NROW(x)
}

# check for diagonal matrix

is.diag <- function(x, tol = 1e-12) {
  dx <- diag(diag(x))
  is.square(x) && all(abs(dx - x) < tol)
}
