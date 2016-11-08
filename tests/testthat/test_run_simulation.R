nsim <- 10
test_that("run simulation", {
  expect_equal(nrow(run_sim(nsim = nsim)), 5 * nsim)
})
