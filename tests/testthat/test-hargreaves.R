test_that("Harvreaves works as expected.", {
  expect_equal(etr_hargreaves(5, 10, 2, 0.5, 150), 1.417284, tolerance = 1e-3)
})
