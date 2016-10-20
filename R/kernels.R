

rempan <- function(n) {
  replicate(n, {
    u <- runif(3) * 2 - 1
    au <- abs(u)
    if (au[3] >= au[2] && au[3] >= au[1]) u[2] else u[3]
  })
}

rtriang <- function(n) {
  runif(n) - runif(n)
}

rrect <- function(n) {
  runif(n) * 2 - 1
}

rbiweight <- function(n) {
  rbeta(n, 3, 3)
}

rtriweight <- function(n) {
  rbeta(n, 4, 4)
}

# pi/4 * cos(pi/2 * x)
rcosine <- function(n) {
  u <- runif(n)
  2*acos(2*u-1)/pi
}


rsmvnorm <- function(n, k = length(sd), sd = rep(1, k)) {
  k <- length(sd)
  matrix(rnorm(n*k, sd = sd), n, k, byrow = TRUE)
}
