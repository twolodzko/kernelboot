

expect_silent(bw.scott(mtcars))
expect_silent(bw.silv(mtcars))

expect_error(bw.scott(1:10))
expect_error(bw.silv(1:10))
