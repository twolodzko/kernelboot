
test_that("multivariate kernels", {

  dat <- mtcars
  n <- nrow(dat)
  m <- ncol(dat)
  w <- rep(1/n, n)
  H <- cov(dat)*0.5

  kernels <- c("gaussian", "epanechnikov", "rectangular",
               "triangular", "biweight", "cosine", "optcosine")

  for (k in kernels) {

    expect_silent(rmvk(10, dat, kernel = k))
    expect_warning(expect_true(all(is.na(rmvk(10, matrix(0, 0L, 1L), kernel = k, bw = 1)))))
    expect_warning(expect_true(all(is.na(rmvk(10, matrix(0, 1L, 0L), kernel = k, bw = 1)))))

  }

  expect_silent(rmvk(10, as.matrix(dat)))

  expect_silent(rmvk(10, dat[1,], bw = 1))
  expect_silent(rmvk(10, dat[1, , drop = FALSE], bw = 1))

  expect_silent(rmvk(10, dat, bw = 1))

  H <- diag(cov(dat))
  expect_silent(rmvk(10, dat, bw = H))

  H <- matrix(1, 2, 2)
  expect_error(rmvk(10, dat, bw = H))

  H <- matrix(NA, m, m)
  expect_error(rmvk(10, dat, bw = H))

  H <- matrix(Inf, m, m)
  expect_error(rmvk(10, dat, bw = H))

  # expect_error(H <- matrix(0, m, m)
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

})
