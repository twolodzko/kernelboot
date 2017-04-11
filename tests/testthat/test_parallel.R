
# this test behaves strangely on Windows
# expect_silent(kernelboot(mtcars, function(data) coef(lm(mpg ~ ., data = data)), R = 10,
#                          parallel = TRUE))
