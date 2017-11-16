context("Test test statistics")

## next_jump.bsFs is correct

test_that("next_jump.bsFs works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)

  res <- next_jump(fit, y)

  expect_true(class(res) == class(fit))
})

test_that("next_jump.bsFs gives a model with one more jump", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)

  res <- next_jump(fit, y)

  jump_before <- binSegInf::jumps(fit)
  jump_after <- binSegInf::jumps(res)

  expect_true(sum(!jump_after %in% jump_before) == 1)
})

######################

## .jump_contrast is correct

test_that(".jump_contrast works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit1 <- binSegInf::binSeg_fixedSteps(y, 1)
  fit2 <- next_jump(fit1, y)

  res <- .jump_contrast(y, binSegInf::jumps(fit1), binSegInf::jumps(fit2))

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

test_that(".jump_contrast computes the correct value", {
  set.seed(10)
  y <- c(rep(0, 5), rep(1, 5), rep(10, 10))
  fit1 <- binSegInf::binSeg_fixedSteps(y, 1)
  fit2 <- next_jump(fit1, y)

  res <- .jump_contrast(y, binSegInf::jumps(fit1), binSegInf::jumps(fit2))

  expect_true(res == 1)
})

#################

## next_jump_statistic is correct

test_that("next_jump_statistic works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)

  res <- next_jump_statistic(y, fit)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

##################

## changepoint_variance is correct

test_that("changepoint_variance works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + rnorm(20)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)

  res <- changepoint_variance(y, fit)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res >= 0)
})
