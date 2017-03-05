

x <- -10:10
y <- rnorm(10)
bw <- bw.nrd0(y)
n <- length(y)
w <- rep(1/n, n)
kernel <- "optcosine"
a <- 1
lp <- FALSE

expect_silent(duvk(x, y))
expect_silent(duvk(x, y, bw, kernel, w, a, lp))

lp <- TRUE
expect_silent(duvk(x, y, bw, kernel, w, a, lp))

lp <- FALSE
bw <- 1:10
expect_message(duvk(x, y, bw, kernel, w, a, lp))

bw <- -1
expect_error(duvk(x, y, bw, kernel, w, a, lp))

bw <- NA
expect_error(duvk(x, y, bw, kernel, w, a, lp))

bw <- Inf
expect_error(duvk(x, y, bw, kernel, w, a, lp))

bw <- 1
kernel <- "fluffy cat"
expect_error(duvk(x, y, bw, kernel, w, a, lp))

w <- 1
kernel <- "gaussian"
expect_silent(duvk(x, y, bw, kernel, w, a, lp))

w <- rep(-1, n)
expect_error(duvk(x, y, bw, kernel, w, a, lp))

w <- rep(Inf, n)
expect_error(duvk(x, y, bw, kernel, w, a, lp))

w <- rep(NA, n)
expect_error(duvk(x, y, bw, kernel, w, a, lp))

w <- rep(1/n, n-1)
expect_error(duvk(x, y, bw, kernel, w, a, lp))

w <- rep(1/n, n+1)
expect_error(duvk(x, y, bw, kernel, w, a, lp))

w <- rep(1/n, n)
a <- -1
expect_error(duvk(x, y, bw, kernel, w, a, lp))

w <- rep(1/n, n)
a <- NA
expect_error(duvk(x, y, bw, kernel, w, a, lp))

w <- rep(1/n, n)
a <- Inf
expect_error(duvk(x, y, bw, kernel, w, a, lp))


