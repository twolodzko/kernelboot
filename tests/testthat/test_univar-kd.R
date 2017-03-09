

y <- rnorm(10)
bw <- bw.nrd0(y)
n <- length(y)
w <- rep(1/n, n)
kernel <- "optcosine"
a <- 1
lp <- FALSE

lp <- TRUE
expect_silent(ruvk(10, y, bw, kernel, w, a, lp))

lp <- FALSE
bw <- 1:10
expect_message(ruvk(10, y, bw, kernel, w, a, lp))

bw <- -1
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

bw <- NA
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

bw <- Inf
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

bw <- 1
kernel <- "fluffy cat"
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

w <- 1
kernel <- "gaussian"
expect_silent(ruvk(10, y, bw, kernel, w, a, lp))

w <- rep(-1, n)
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

w <- rep(Inf, n)
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

w <- rep(NA, n)
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

w <- rep(1/n, n-1)
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

w <- rep(1/n, n+1)
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

w <- rep(1/n, n)
a <- -1
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

w <- rep(1/n, n)
a <- NA
expect_error(ruvk(10, y, bw, kernel, w, a, lp))

w <- rep(1/n, n)
a <- Inf
expect_error(ruvk(10, y, bw, kernel, w, a, lp))


