context("Test geometry")

## .quadratic is correct

test_that(".quadratic works", {
  res <- .quadratic(-1, 2, 5)

  expect_true(length(res) == 2)
  expect_true(res[1] < res[2])
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(all(!is.na(res)))
})

test_that(".quadratic can return only one number", {
  res <- .quadratic(1, 2, 1)

  expect_true(length(res) == 1)
  expect_true(!is.na(res))
})

test_that(".quadratic can return only NA", {
  res <- .quadratic(1, 2, 5)

  expect_true(length(res) == 1)
  expect_true(is.na(res))
})
