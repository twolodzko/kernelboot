

silv86 <- function(x) {
  d <- ncol(x)
  H <- (4/(d+2))^(1/(d+4))  * n^(-1/(d+4))
  H^2
}

scott <- function(x) {
  d <- ncol(x)
  H <- n^(-1/(d+4))
  H^2
}
