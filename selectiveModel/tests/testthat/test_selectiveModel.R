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
  res <- selected_model_inference(y, fit_method, num_samp = 10, verbose = F)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res >= 0)
  expect_true(res <= 1)
})

test_that("selected_model_inference works reasonably relatively for hit run", {
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  set.seed(10)
  y1 <- c(rep(0, 10), rep(5, 10)) + 0.01*rnorm(20)
  res1 <- selected_model_inference(y1, fit_method, num_samp = 50, verbose = F)

  set.seed(10)
  y2 <- c(rep(0, 5), rep(1, 5), rep(5, 10)) + 0.01*rnorm(20)
  res2 <- selected_model_inference(y2, fit_method, num_samp = 50, verbose = F)

  expect_true(res1 > res2)
})

test_that("selected_model_inference works reasonably relatively for rejection", {
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  set.seed(10)
  y1 <- c(rep(0, 5), rep(5, 5)) + 0.01*rnorm(20)
  res1 <- selected_model_inference(y1, fit_method, sample_method = "rejection",
                                   num_samp = 15, verbose = F)

  set.seed(10)
  y2 <- c(rep(0, 2), rep(1, 3), rep(5, 5)) + 0.01*rnorm(10)
  res2 <- selected_model_inference(y2, fit_method, sample_method = "rejection",
                                   num_samp = 15, verbose = F)

  expect_true(res1 > res2)
})
