

fun1 <- function(data) coef(lm(mpg ~ ., data = data))
fun2 <- function(data) mean(data)
fun3 <- function(data) data
fun_err <- function(data) stop()


dat <- mtcars

expect_message(expect_error(kernelboot(dat, fun_err, R = 10)))

expect_silent(kernelboot(dat, fun1, R = 10))
expect_silent(kernelboot(dat$mpg, fun2, R = 10))
expect_silent(kernelboot(dat$mpg, fun3, R = 10))

expect_silent(kernelboot(dat, fun1, R = 10, parallel = TRUE))


expect_error(kernelboot(as.list(1:10), fun2, R = 10))

dat[5, 1] <- NA

expect_error(kernelboot(dat, fun1, R = 10))
expect_error(kernelboot(dat$mpg, fun2, R = 10))
expect_error(kernelboot(dat$mpg, fun3, R = 10))


dat <- mtcars
dat[5, 1] <- Inf

expect_error(kernelboot(dat, fun1, R = 10))
expect_error(kernelboot(dat$mpg, fun2, R = 10))
expect_error(kernelboot(dat$mpg, fun3, R = 10))


dat <- mtcars
dat[] <- 0

expect_error(kernelboot(dat, fun1, R = 10))
expect_silent(kernelboot(dat$mpg, fun2, R = 10))
expect_silent(kernelboot(dat$mpg, fun3, R = 10))

dat <- mtcars
n <- nrow(dat)
k <- ncol(dat)
w <- rep(1, n)

expect_silent(kernelboot(dat, fun1, weights = w, R = 10))
expect_silent(kernelboot(dat$mpg, weights = w, fun2, R = 10))
expect_silent(kernelboot(dat$mpg, weights = w, fun3, R = 10))


w <- 1

expect_silent(kernelboot(dat, fun1, weights = w, R = 10))
expect_silent(kernelboot(dat$mpg, weights = w, fun2, R = 10))
expect_silent(kernelboot(dat$mpg, weights = w, fun3, R = 10))


w <- rep(NA, n)

expect_error(kernelboot(dat, fun1, weights = w, R = 10))
expect_error(kernelboot(dat$mpg, weights = w, fun2, R = 10))
expect_error(kernelboot(dat$mpg, weights = w, fun3, R = 10))


w <- rep(-Inf, n)

expect_error(kernelboot(dat, fun1, weights = w, R = 10))
expect_error(kernelboot(dat$mpg, weights = w, fun2, R = 10))
expect_error(kernelboot(dat$mpg, weights = w, fun3, R = 10))


w <- rep(1/n, n-1)

expect_error(kernelboot(dat, fun1, weights = w, R = 10))
expect_error(kernelboot(dat$mpg, weights = w, fun2, R = 10))
expect_error(kernelboot(dat$mpg, weights = w, fun3, R = 10))


w <- rep(1/n, n+1)

expect_error(kernelboot(dat, fun1, weights = w, R = 10))
expect_error(kernelboot(dat$mpg, weights = w, fun2, R = 10))
expect_error(kernelboot(dat$mpg, weights = w, fun3, R = 10))


expect_silent(kernelboot(dat, fun1, bw = diag(k), R = 10))
expect_silent(kernelboot(dat, fun1, bw = rep(1, k), R = 10))
expect_silent(kernelboot(dat$mpg, bw = 1, fun2, R = 10))
expect_silent(kernelboot(dat$mpg, bw = 1, fun3, R = 10))


expect_error(kernelboot(dat, fun1, bw = rep(-1, k), R = 10))
expect_error(kernelboot(dat$mpg, bw = -1, fun2, R = 10))
expect_error(kernelboot(dat$mpg, bw = -1, fun3, R = 10))


expect_error(kernelboot(dat, fun1, bw = matrix(NA, k, k), R = 10))
expect_error(kernelboot(dat, fun1, bw = rep(NA, k), R = 10))
expect_error(kernelboot(dat$mpg, bw = NA, fun2, R = 10))
expect_error(kernelboot(dat$mpg, bw = NA, fun3, R = 10))


expect_error(kernelboot(dat, fun1, bw = matrix(Inf, k, k), R = 10))
expect_error(kernelboot(dat, fun1, bw = rep(Inf, k), R = 10))
expect_error(kernelboot(dat$mpg, bw = Inf, fun2, R = 10))
expect_error(kernelboot(dat$mpg, bw = Inf, fun3, R = 10))


a <- 1

expect_silent(kernelboot(dat, fun1, bw = diag(k), adjust = a, R = 10))
expect_silent(kernelboot(dat, fun1, bw = rep(1, k), adjust = a, R = 10))
expect_silent(kernelboot(dat$mpg, bw = 1, fun2, adjust = a, R = 10))
expect_silent(kernelboot(dat$mpg, bw = 1, fun3, adjust = a, R = 10))



a <- -1

expect_error(kernelboot(dat, fun1, bw = diag(k), adjust = a, R = 10))
expect_error(kernelboot(dat, fun1, bw = rep(1, k), adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun2, adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun3, adjust = a, R = 10))



a <- NA

expect_error(kernelboot(dat, fun1, bw = diag(k), adjust = a, R = 10))
expect_error(kernelboot(dat, fun1, bw = rep(1, k), adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun2, adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun3, adjust = a, R = 10))



a <- Inf

expect_error(kernelboot(dat, fun1, bw = diag(k), adjust = a, R = 10))
expect_error(kernelboot(dat, fun1, bw = rep(1, k), adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun2, adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun3, adjust = a, R = 10))


a <- 1:10

expect_silent(kernelboot(dat, fun1, bw = diag(k), adjust = a, R = 10))
expect_silent(kernelboot(dat, fun1, bw = rep(1, k), adjust = a, R = 10))
expect_silent(kernelboot(dat$mpg, bw = 1, fun2, adjust = a, R = 10))
expect_silent(kernelboot(dat$mpg, bw = 1, fun3, adjust = a, R = 10))


a <- matrix(1, 10, 10)

expect_error(kernelboot(dat, fun1, bw = diag(k), adjust = a, R = 10))
expect_error(kernelboot(dat, fun1, bw = rep(1, k), adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun2, adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun3, adjust = a, R = 10))


a <- as.data.frame(matrix(1, 10, 10))

expect_error(kernelboot(dat, fun1, bw = diag(k), adjust = a, R = 10))
expect_error(kernelboot(dat, fun1, bw = rep(1, k), adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun2, adjust = a, R = 10))
expect_error(kernelboot(dat$mpg, bw = 1, fun3, adjust = a, R = 10))


dat <- mtcars

expect_silent(kernelboot(dat, fun1, R = 10, kernel = "cosine"))
expect_silent(kernelboot(dat$mpg, fun2, R = 10, kernel = "cosine"))
expect_silent(kernelboot(dat$mpg, fun3, R = 10, kernel = "cosine"))


dat <- mtcars

# expect_message(kernelboot(dat, fun1, R = 10, ignore = colnames(dat)))
expect_silent(kernelboot(dat$mpg, fun2, R = 10, ignore = colnames(dat)))
expect_silent(kernelboot(dat$mpg, fun3, R = 10, ignore = colnames(dat)))

