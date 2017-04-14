


test_that("multivariate Gaussian kernel", {

  dat <- mtcars
  n <- nrow(dat)
  m <- ncol(dat)
  w <- rep(1/n, n)
  H <- cov(dat)*0.5

  expect_silent(rmvg(10, as.matrix(dat)))

  expect_silent(rmvg(10, dat[1, ], bw = 1))
  expect_silent(rmvg(10, dat[1, , drop = FALSE], bw = 1))

  expect_silent(rmvg(10, dat, bw = 1))

  H <- diag(cov(dat))
  expect_silent(rmvg(10, dat, bw = H))

  H <- matrix(1, 2, 2)
  expect_error(rmvg(10, dat, bw = H))

  H <- matrix(NA, m, m)
  expect_error(rmvg(10, dat, bw = H))

  H <- matrix(Inf, m, m)
  expect_error(rmvg(10, dat, bw = H))

  expect_silent(rmvg(10, dat, weights = w))

  expect_silent(rmvg(10, dat, weights = 1))

  w <- c(1,1,1)
  expect_error(rmvg(10, dat, weights = w))

  w <- rep(-1, n)
  expect_error(rmvg(10, dat, weights = w))

  w <- rep(NA, n)
  expect_error(rmvg(10, dat, weights = w))

  w <- rep(Inf, n)
  expect_error(rmvg(10, dat, weights = w))

  expect_silent(rmvg(10, dat, adjust = 1))

  expect_silent(rmvg(10, dat, adjust = 1:10))

  expect_error(rmvg(10, dat, adjust = NA))

  expect_error(rmvg(10, dat, adjust = 0))

  expect_error(rmvg(10, dat, adjust = Inf))


  ## Not run:

  if ( requireNamespace("cramer", quietly = TRUE) ) {

    library(cramer)

    set.seed(0xBEEF)

    bw <- diag(bw.silv(mtcars))
    x <- rmvk(500, mtcars, kernel = "gaussian", bw = sqrt(bw))
    y <- rmvg(500, mtcars, bw = bw)

    stopifnot(dim(x) == dim(y))

    stopifnot(cramer.test(x, y)$result == 0)

  }

  ## End(Not run)

})

