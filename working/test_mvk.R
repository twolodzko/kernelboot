

rm(list = ls())
library(Rcpp)
library(ggplot2)

sourceCpp("src/multivar-kde.cpp")


vars <- c("mpg", "disp")
X <- as.matrix(mtcars[, vars])

ggplot(as.data.frame(X), aes_string(x = vars[1], y = vars[2])) +
  geom_density2d() + geom_jitter(shape = 2) + theme_minimal()


cov(X)
kernelboot::bw.scott(X)


samp <- cpp_rmvkde(1000, X, kernelboot::bw.scott(X), 1)
Y <- as.data.frame(samp$sample)
colnames(Y) <- vars

ggplot(as.data.frame(X), aes_string(x = vars[1], y = vars[2])) +
  geom_density2d(color = "lightgray") +
  geom_jitter(data = Y, aes_string(x = vars[1], y = vars[2]),
              color = "lightblue2", alpha = 0.5) +
  geom_jitter(shape = 2) +
  theme_minimal()


Y$p <- drop(cpp_dmvkde(as.matrix(Y[,1:2]), X, kernelboot::bw.scott(X), 1)$density)

ggplot(as.data.frame(X), aes_string(x = vars[1], y = vars[2])) +
  geom_density2d(color = "lightgray") +
  geom_jitter(data = Y, aes_string(x = vars[1], y = vars[2], colour = "p"),
              alpha = 0.5) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  geom_density2d(color = "lightgray") +
  geom_jitter(shape = 2) +
  # geom_point(data = as.data.frame(t(colMeans(X))),
  #            aes_string(x = vars[1], y = vars[2]),
  #            color = "red") +
  theme_minimal()



#
#
#
# samp <- cpp_rmvkde(1000, X, cov(X)*0.1, 1)
# Y <- as.data.frame(samp$sample)
# colnames(Y) <- vars
#
# ggplot(as.data.frame(X), aes_string(x = vars[1], y = vars[2])) +
#   geom_density2d(color = "lightgray") +
#   geom_jitter(data = Y, aes_string(x = vars[1], y = vars[2]),
#               color = "lightblue2", alpha = 0.5) +
#   geom_jitter(shape = 2) +
#   theme_minimal()
