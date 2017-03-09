
set.seed(123)
rm(list = ls())

dat <- mtcars
n <- nrow(dat)
k <- ncol(dat)
w <- rep(1/n, n)
H <- cov(dat)*0.5

expect_silent(rmvk(10, dat))

expect_silent(rmvk(10, as.matrix(dat)))

# this IS an inappropriate input ?
expect_silent(rmvk(10, dat[1,], bw = 1))

expect_silent(rmvk(10, dat, bw = 1))

H <- diag(cov(dat))
expect_silent(rmvk(10, dat, bw = H))

H <- matrix(1, 2, 2)
expect_error(rmvk(10, dat, bw = H))

H <- matrix(NA, k, k)
expect_error(rmvk(10, dat, bw = H))

H <- matrix(Inf, k, k)
expect_error(rmvk(10, dat, bw = H))

# expect_error(H <- matrix(0, k, k)
# expect_error(rmvk(10, dat, bw = H))

expect_silent(rmvk(10, dat, weights = w))

expect_silent(rmvk(10, dat, weights = 1))

w <- c(1,1,1)
expect_error(rmvk(10, dat, weights = w))

w <- rep(-1, n)
expect_error(rmvk(10, dat, weights = w))

w <- rep(NA, n)
expect_error(rmvk(10, dat, weights = w))

w <- rep(Inf, n)
expect_error(rmvk(10, dat, weights = w))

expect_silent(rmvk(10, dat, adjust = 1))

expect_silent(rmvk(10, dat, adjust = 1:10))

expect_error(rmvk(10, dat, adjust = NA))

expect_error(rmvk(10, dat, adjust = 0))

expect_error(rmvk(10, dat, adjust = Inf))

