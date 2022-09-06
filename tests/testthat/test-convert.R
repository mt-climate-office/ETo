test_that("Conversion to radians works properly.", {
  expect_equal(lat_to_radians(90), pi / 2)
})

test_that("Out of bounds latitude is not allowed.", {
  expect_error(lat_to_radians(361))
  expect_error(lat_to_radians(-361))
})

test_that("Only numeric inputs are allowed.", {
  expect_error(lat_to_radians("a"))
})
