
set.seed(123)

X <- matrix(rnorm(3*20), 20, 3)
mu <- c(0, 0, 0)
Sigma <- matrix(c(1, 0, 0,
                  0, 1, 0,
                  0, 0, 1), 3, 3)

expect_equal(dmvn(X, mu, Sigma), dmvk(X, matrix(mu, 1), Sigma), tolerance = 1e-6)
expect_equal(dmvn(X, mu, Sigma), dnorm(X[,1])*dnorm(X[,2])*dnorm(X[,3]), tolerance = 1e-6)
expect_equal(dmvpk(X, mu, Sigma, kernel = "gaussian"), dnorm(X[,1])*dnorm(X[,2])*dnorm(X[,3]), tolerance = 1e-6)
