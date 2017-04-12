
test_that("single column multivariate kernels = univariate kernels", {

  set.seed(1)

  kernels <- c("gaussian", "epanechnikov", "rectangular",
               "triangular", "biweight", "cosine", "optcosine")

  N <- 1e5
  dat <- mtcars[, 1, drop = FALSE]

  for (k in kernels) {

    x <- ruvk(N, drop(dat[, 1]), kernel = k, bw = 2, shrinked = FALSE)
    y <- drop(rmvk(N, dat,       kernel = k, bw = 2, shrinked = FALSE))
    expect_true(suppressWarnings(ks.test(x, y)$p.value > 0.05))

    x <- ruvk(N, drop(dat[, 1]), kernel = k, bw = 2, shrinked = TRUE)
    y <- drop(rmvk(N, dat,       kernel = k, bw = 2, shrinked = TRUE))
    expect_true(suppressWarnings(ks.test(x, y)$p.value > 0.05))

  }

})




