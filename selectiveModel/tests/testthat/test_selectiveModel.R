context("Test selective model")

## .segments is correct

test_that(".segments works", {
  res <- .segments(10, c(3, 7))

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(3, 10)))
})

test_that(".segments forms the proper matrix", {
  res <- .segments(10, c(3, 7))

  expect_true(all(rowSums(res) == 1))
  for(i in 1:3){
    expect_true(length(unique(res[i,])) == 2)
  }
  expect_true(all(res[1,1:3] != 0))
  expect_true(all(res[1,4:10] == 0))
  expect_true(all(res[2,4:7] != 0))
  expect_true(all(res[2,c(1:3,8:10)] == 0))
  expect_true(all(res[3,8:10] != 0))
  expect_true(all(res[3,1:7] == 0))
})

test_that(".segments can ignore properly", {
  res <- .segments(10, c(3, 7), ignore_jump = 1)

  expect_true(all(dim(res) == c(2, 10)))
  expect_true(all(res[1,] == c(rep(1/7,7), rep(0,3))))
})

########################

## .segment_means is correct

test_that(".segment_means works", {
  set.seed(10)
  y <- rnorm(20)
  segments <- .segments(length(y), c(5, 10, 15))
  res <- .segment_means(y, segments)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 4)
})

test_that(".segment_means returns the proper means", {
  set.seed(10)
  y <- rnorm(20)
  segments <- .segments(length(y), c(5, 10, 15))
  res <- .segment_means(y, segments)

  target <- c(mean(y[1:5]), mean(y[6:10]), mean(y[11:15]), mean(y[16:20]))

  expect_true(sum(abs(res - target)) < 1e-6)
})

############################

## selected_model_inference is correct

test_that("selected_model_inference works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}
  res <- selected_model_inference(y, fit_method, num_samp = 10, verbose = F,
                                  param = list(burn_in = 10))

  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("pval", "test_stat", "null_stat")))
  expect_true(is.numeric(res$pval))
  expect_true(!is.matrix(res$pval))
  expect_true(length(res$pval) == 1)
  expect_true(res$pval >= 0)
  expect_true(res$pval <= 1)
})

test_that("selected_model_inference works for known sigma", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}
  res <- selected_model_inference(y, fit_method, sigma = 1, num_samp = 10, verbose = F,
                                  param = list(burn_in = 10))

  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("pval", "test_stat", "null_stat")))
  expect_true(is.numeric(res$pval))
  expect_true(!is.matrix(res$pval))
  expect_true(length(res$pval) == 1)
  expect_true(res$pval >= 0)
  expect_true(res$pval <= 1)
})

test_that("selected_model_inference works reasonably relatively for hit run", {
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  set.seed(10)
  y1 <- c(rep(0, 10), rep(5, 10)) + 0.01*rnorm(20)
  res1 <- selected_model_inference(y1, fit_method, num_samp = 50, verbose = F,
                                   param = list(burn_in = 10))

  set.seed(10)
  y2 <- c(rep(0, 5), rep(1, 5), rep(5, 10)) + 0.01*rnorm(20)
  res2 <- selected_model_inference(y2, fit_method, num_samp = 50, verbose = F,
                                   param = list(burn_in = 10))

  expect_true(res1$pval > res2$pval)
})

test_that("selected_model_inference works reasonably relatively for rejection", {
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  set.seed(10)
  y1 <- c(rep(0, 5), rep(5, 5)) + 0.01*rnorm(20)
  res1 <- selected_model_inference(y1, fit_method, sample_method = "rejection",
                                   num_samp = 15, verbose = F,
                                   param = list(burn_in = 10))

  set.seed(10)
  y2 <- c(rep(0, 2), rep(1, 3), rep(5, 5)) + 0.01*rnorm(10)
  res2 <- selected_model_inference(y2, fit_method, sample_method = "rejection",
                                   num_samp = 15, verbose = F,
                                   param = list(burn_in = 10))

  expect_true(res1$pval > res2$pval)
})

test_that("selected_model_inference works with the variance test statistic", {
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  set.seed(10)
  y1 <- c(rep(0, 5), rep(5, 5)) + 0.01*rnorm(20)
  res1 <- selected_model_inference(y1, fit_method, test_func = changepoint_variance,
                                   num_samp = 15, verbose = F,
                                   param = list(burn_in = 10))

  set.seed(10)
  y2 <- c(rep(0, 2), rep(1, 3), rep(5, 5)) + 0.01*rnorm(10)
  res2 <- selected_model_inference(y2, fit_method, test_func = changepoint_variance,
                                   num_samp = 15, verbose = F,
                                   param = list(burn_in = 10))

  expect_true(res1$pval > res2$pval)
})

test_that("selected_model_inference works for one problem case", {
  set.seed(77)
  y <- rnorm(20)
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  res <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                  num_samp = 50,
                                  param = list(burn_in = 10, time_limit = 600))

  expect_true(length(res) == 3)
})

test_that("selected_model_inference works for one problem case", {
  set.seed(238)
  y <- rnorm(20)
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  res <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                  num_samp = 50,
                                  param = list(burn_in = 10, time_limit = 600))

  expect_true(length(res) == 3)
})

