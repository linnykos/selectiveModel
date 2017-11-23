context("Test Gaussian")

## .transform_gaussian is correct

test_that(".transform_gaussian works", {
  set.seed(10)
  cov_mat <- diag(10)
  cov_mat[1:5,1:5] <- 0.9
  cov_mat[6:10,6:10] <- 0.9
  diag(cov_mat) <- 1
  gaussian <- .gaussian(rep(0, 10), cov_mat)

  vec1 <- rnorm(10); vec2 <- rnorm(10)
  rotation <- .rotation_matrix(vec1, vec2)
  shift <- rnorm(10)

  res <- .transform_gaussian(gaussian, shift = shift, rotation = rotation)

  expect_true(class(res) == "gaussian")
  expect_true(length(res$mean) == length(gaussian$mean))
  expect_true(all(dim(res$covariance) == dim(gaussian$covariance)))
})


############################

## .conditional_gaussian is correct

test_that(".conditional_gaussian works", {
  set.seed(10)
  cov_mat <- diag(10)
  cov_mat[1:5,1:5] <- 0.9
  cov_mat[6:10,6:10] <- 0.9
  diag(cov_mat) <- 1
  gaussian <- .gaussian(rep(0, 10), cov_mat)

  res <- .conditional_gaussian(gaussian, c(5,10))

  expect_true(class(res) == "gaussian")
  expect_true(length(res$mean) == 8)
  expect_true(all(dim(res$covariance) == 8))
})
