test_sample <- create_sample()
test_that("match_reg", {
  expect_error(interaction_reg(test_sample, method = "oops"),
               "method has no valid input. It should be one of the following: traditional, augmented.")
  expect_equal(length(format_reg(match_reg(test_sample))), 5)
  expect_equal(length(format_reg(interaction_reg(test_sample))), 5)
  expect_equal(length(format_reg(interaction_reg(test_sample,
                                                 method = "augmented"))), 5)
})
