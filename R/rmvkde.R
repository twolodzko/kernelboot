
library(ks)
library(mvtnorm)

expand.vector <- function(x, n = 1) {
  matrix(x, ncol = length(x), nrow = n, byrow = TRUE)
}

n <- 10
x <- mtcars
nx <- nrow(x)
kx <- ncol(x)
mx <- colMeans(x)
covariance <- cov(x)
bw <- 0.5
H <- sqrt(bw*diag(covariance))

idx <- sample.int(nx, n, replace = TRUE)
x[idx,] + rsmvnorm(n, kx, sd = H)


rmvkde <- function(n, x, bw = "silverman", covariance = cov(x),
                   kernel = c("gaussian", "mvn"),
                   method = c("default", "shrinked"),
                   weights = NULL) {

  method <- match.arg(method)
  kernel <- match.arg(kernel)

  if (anyNA(x))
    stop("'x' contains missing values")

  nx <- nrow(x)
  kx <- ncol(x)

  if (is.character(bw)) {
    bw <- switch (tolower(bw),
                  silverman = silv86(x),
                  scott     = scott(x),
                  hns       = ks::Hns(x),
                  stop("unknown bandwidth rule"))
  } else if (is.numeric(bw)) {
    if (!is.matrix(bw) && length(bw) > ncol(x))
      stop("incorrect 'bw' format")
    if (!all(is.finite(bw)))
      stop("non-finite 'bw'")
  }
  H <- bw * covariance * diag(kx)
  idx <- sample.int(nx, n, replace = TRUE)

  eps <- switch (tolower(kernel),
                 "gaussian" = rsmvnorm(n, kx, sd = sqrt(diag(H))),
                 mvtnorm::rmvnorm(n, sigma = H))

  if (method == "shrinked") {
    mx <- expand.vector(colMeans(x), n)
    sx <- sqrt(expand.vector(diag(covariance), n))
    h <- if (is.numeric(h) & length(h == 1)) h else expand.vector(diag(H), n)
    (x[idx,] - mx + eps)/sqrt(1+h) + mx
  } else {
    x[idx,] + eps
  }

}










library(mvtnorm)
library(ks)

n <- 100
x <- mtcars
h <- 1

# http://vita.had.co.nz/papers/density-estimation.pdf
# http://www.math.uah.edu/stat/special/MultiNormal.html
# https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation#cite_note-WJ1995-5

# + dodatkowy parametr h

# ?KernSmooth::bkde
# kernel
# character string which determines the smoothing kernel. kernel can be:
#   "normal" - the Gaussian density function (the default).
# "box" - a rectangular box.
# "epanech" - the centred beta(2,2) density.
# "biweight" - the centred beta(3,3) density.
# "triweight" - the centred beta(4,4) density. This can be abbreviated to any unique abbreviation.

# Kernels RNG

#' Wand, M.P and Jones, M.C. (1995). Kernel Smoothing. London: Chapman & Hall/CRC.
#'
#' Silverman, B.W. (1986). Density Estimation for Statistics and Data Analysis. Chapman & Hall/CRC.
#'
#' Devroye, L. and Gyorfi, L. (1985) Nonparametric Density Estimation: the L1 view. Wiley.

# selecting bw






x <- as.matrix(mtcars)
adjust <- 1
bw <- "hns"
method <- "smoothed"

n <- nrow(x)
k <- ncol(x)
idx <- sample.int(n, n, replace = TRUE)
mx <- colMeans(x)
S <- cov(x)
sx <- sqrt(diag(S))
mxMtx <- matrix(mx, n, k, byrow = TRUE)
#sxMtx <- matrix(s, n, K, byrow = TRUE)

H <- switch(tolower(bw),
            silverman = silv86(x) * sx,
            scott     = scott(x)  * sx,
            hns       = ks::Hns(x),
            stop("unknown bandwidth rule"))

H <- adjust * H

if (NCOL(H) == 1) {
  hMtx <- matrix(H, n, k, byrow = TRUE)
  heps <- rsmvnorm(n, sd = H)
} else {
  hMtx <- matrix(diag(H), n, k, byrow = TRUE)
  heps <- mvtnorm::rmvnorm(n, sigma = H)
}

y <- switch (method,
             "smoothed" = mxMtx + (x[idx,] - mxMtx + heps)/sqrt(1 + hMtx^2),
             x[idx,] + heps
)



#
# kernelboot <- function(data, statistic, R = 500,
#                        bw, adjust = 1, method = c("default", "smoothed"), ...) {
#
#   call <- match.call()
#   method <- match.arg(method)
#
#   tryCatch(stat.result <- statistic(data, ...),
#            error = function(e) {
#              message("Applying the statistic on the original data resulted in an error")
#              stop(e)
#            }
#   )
#
#   if (is.matrix(data) || is.data.frame(data)) {
#
#     if (missing(bw))
#       bw <- "silverman"
#
#     if (is.character(bw)) {
#       h <- switch(tolower(bw),
#                    silverman = silv86(data),
#                    scott     = scott(data),
#                    hns       = ks::Hns(data),
#                    stop("unknown bandwidth rule"))
#     } else {
#       h <- bw
#     }
#
#     h <- adjust * h
#
#     n <- nrow(data)
#     k <- ncol(data)
#     mx <- colMeans(data)
#     mxMtx <- matrix(mx, n, k, byrow = TRUE)
#
#     if (bw == "hns" || is.numeric(bw)) {
#
#       hMtx <- h
#
#       res <- replicate(R, {
#         idx <- sample.int(n, n, replace = TRUE)
#         heps <- mvtnorm::rmvnorm(n, sigma = hMtx)
#         new.data <- switch (method,
#                             "smoothed" = mxMtx + (data[idx, , drop = FALSE] - mxMtx + heps)/sqrt(1 + hMtx^2),
#                             data[idx, , drop = FALSE] + heps
#         )
#         statistic(new.data, ...)
#       })
#
#     } else {
#
#       hMtx <- matrix(diag(h), n, k, byrow = TRUE)
#
#       res <- replicate(R, {
#         idx <- sample.int(n, n, replace = TRUE)
#         heps <- rsmvnorm(n, sd = diag(h))
#         new.data <- switch (method,
#                             "smoothed" = mxMtx + (data[idx, , drop = FALSE] - mxMtx + heps)/sqrt(1 + hMtx^2),
#                             data[idx, , drop = FALSE] + heps
#         )
#         statistic(new.data, ...)
#       })
#
#     }
#
#
#
#   } else {
#
#     if (missing(bw))
#       bw <- "nrd0"
#
#     if (is.character(bw)) {
#       bw <- switch(tolower(bw),
#                    nrd0 = bw.nrd0(data),
#                    nrd = bw.nrd(data),
#                    ucv = bw.ucv(data),
#                    bcv = bw.bcv(data), sj = ,
#                    `sj-ste` = bw.SJ(data, method = "ste"),
#                    `sj-dpi` = bw.SJ(data, method = "dpi"),
#                    stop("unknown bandwidth rule"))
#     } else {
#       h <- bw
#     }
#
#     h <- adjust * h
#
#     n <- length(data)
#     mx <- mean(data)
#     sx <- var(data)
#
#     res <- replicate(R, {
#       idx <- sample.int(n, n, replace = TRUE)
#       heps <- rnorm(n, sd = h)
#       new.data <- switch (method,
#               "smoothed" = mx + (data[idx] - mx + heps)/sqrt(1 + h^2/sx),
#               data[idx] + heps
#       )
#       statistic(new.data, ...)
#     })
#
#   }
#
#   structure(
#     res,
#     class = "kernelboot",
#     R = R,
#     bw = bw,
#     adjust = adjust,
#     h = h,
#     method = method,
#     call = call
#   )
#
# }
#
