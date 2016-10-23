


if (kernel == "mvn") {
  if (is.vector(bw) && is.numeric(bw)) {
    if (length(bw) != 1)
      stop("'bw' is not a scalar")
  } else if (is.matrix(bw)) {
    if (ncol(bw) != nrow(bw) || ncol(bw) != d)
      stop("'bw' is not a square matrix")
  } else {
    stop("'bw' should be scalar or matrix")
  }
} else {



} else {
  if (kernel == "mvn") {
    bw <- (4/(d+2))^(1/(d+4)) * n^(-1/(d+4)) * cov(data)
  } else {
    bw <- bw.silv86(data)
  }
}



if (kernel == "mvn" && preserve.var) {

  mx <- apply(data, 2, mean)
  sx <- cov(data)
  mx <- matrix(rep(mx, n), n, d, byrow = TRUE)

  res <- repeatFun(1:R, function(i) {

    idx <- sample.int(n, n, replace = TRUE, prob = weights)
    boot.data <- data[idx, ]
    boot.data <- mx + (boot.data - mx + rmvnorm(n, mx[1,], bw))/sqrt(1 + bw^2)

    statistic(boot.data, ...)

  }, mc.cores = mc.cores)

} else if(kernel == "mvn") {

  res <- repeatFun(1:R, function(i) {

    idx <- sample.int(n, n, replace = TRUE, prob = weights)
    boot.data <- data[idx, ]
    boot.data <- boot.data + rmvnorm(n, mx[1,], bw)

    statistic(boot.data, ...)

  }, mc.cores = mc.cores)

} else
