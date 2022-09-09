test_that("Conversion to radians works properly.", {
  expect_equal(lat_to_radians(90), pi / 2)
})

test_that("Radiation conversion works properly", {
  expect_equal(wm2_to_mj(10), 0.864, tolerance = 1e-3)
})
