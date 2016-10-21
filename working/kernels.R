


sd(rrect(1e6))
sd(rtriang(1e6))
sd(rbiweight(1e6))
sd(rtriweight(1e6))
sd(rempan(1e6))
sd(roptcos(1e6))

mean(abs(rrect(1e6) < 1))
mean(abs(rtriang(1e6) < 1))
mean(abs(rbiweight(1e6) < 1))
mean(abs(rtriweight(1e6) < 1))
mean(abs(rempan(1e6) < 1))
mean(abs(roptcos(1e6) < 1))

hist(rrect(1e6), 100)
hist(rtriang(1e6), 100)
hist(rbiweight(1e6), 100)
hist(rtriweight(1e6), 100)
hist(rempan(1e6), 100)
hist(roptcos(1e6), 100)




# rempan <- function(n) {
#   replicate(n, {
#     u <- runif(3) * 2 - 1
#     au <- abs(u)
#     if (au[3] >= au[2] && au[3] >= au[1]) u[2] else u[3]
#   }) * sqrt(5)
# }
#
# rtriang <- function(n) {
#   (runif(n) - runif(n)) * sqrt(6)
# }

# rrect2 <- function(n) {
#   runif(n, -1, 1) * sqrt(3)
# }

# rbiweight <- function(n) {
#   (rbeta(n, 3, 3) * 2 - 1) * sqrt(7)
# }
#
# rtriweight <- function(n) {
#   (rbeta(n, 4, 4) *2 - 1) * 3
# }

# pi/4 * cos(pi/2 * x)
# roptcosine <- function(n) {
#   u <- runif(n)
#   (2*acos(2*u-1)/pi - 1)/sqrt(1 - 8/pi^2)
# }

# rsmvnorm2 <- function(n, k = length(sd), sd = rep(1, k)) {
#   k <- length(sd)
#   matrix(rnorm(n*k, sd = sd), n, k, byrow = TRUE)
# }
