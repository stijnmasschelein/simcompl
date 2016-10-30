c_test <- c(-1, 1)
test_that("check_numeric", {
  expect_error(check_numeric(1, 2), "1 should have 2 items")
  expect_error(check_numeric("hi", 1), '"hi" is not numerical')
  expect_error(check_numeric(c_test, 2, c(0,1)),
               "c_test does not lie in the interval")
})
test_that("create_sample",{
  expect_equal(any(dim(create_sample() == c(200, 7))), TRUE)
})
