context("Test hit and run, Line")

## .remove_nullspace_gaussian is correct

test_that(".remove_nullspace_gaussian works", {
  set.seed(10)
  n <- 10
  gaussian <- .gaussian(rep(0, n), diag(n))
  y <- rnorm(n)
  segments <- .segments(n, 5)

  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  res <- .remove_nullspace_gaussian(gaussian, segments_full, mean_val)

  expect_true(class(res) == "gaussian")
  expect_true(length(res$mean) + nrow(segments) == length(gaussian$mean))
  expect_true(all(dim(res$covariance) + nrow(segments) == dim(gaussian$covariance)))
})

#######

## .remove_nullspace_polyhedra is correct

test_that(".remove_nullspace_polyhedra works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) +rnorm(20)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)
  polyhedra <- binSegInf::polyhedra(fit)

  segments <- .segments(20, binSegInf::jumps(fit))

  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  res <- .remove_nullspace_polyhedra(polyhedra, segments_full, mean_val)

  expect_true(class(res) == "polyhedra")
  expect_true(all(length(res$u) == length(polyhedra$u)))
  expect_true(all(dim(res$gamma) + c(0, nrow(segments)) == dim(polyhedra$gamma)))
})

########

## .remove_nullspace is correct

test_that(".remove_nullspace works", {
  n <- 10
  set.seed(10)
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)
  polyhedra <- binSegInf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  segments <- .segments(n, binSegInf::jumps(fit))
  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  res <- .remove_nullspace(gaussian, polyhedra, segments_full, mean_val)

  expect_true(is.list(res))
  expect_true("forward_translation" %in% names(res))
  expect_true("backward_translation" %in% names(res))
  expect_true(class(res$gaussian) == "gaussian")
  expect_true(class(res$polyhedra) == "polyhedra")
})

########

## .whiten_polyhedra is correct

test_that(".whiten_polyhedra works", {
  n <- 10
  set.seed(10)
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)
  polyhedra <- binSegInf::polyhedra(fit)

  sqrt_cov <- diag(n)
  mean_vec <- rep(0, n)

  res <- .whiten_polyhedra(polyhedra, sqrt_cov, mean_vec)

  expect_true(class(res) == "polyhedra")
  expect_true(all(length(res$u) == length(polyhedra$u)))
  expect_true(all(dim(res$gamma) == dim(polyhedra$gamma)))
})

#######

## .factor_covariance is correct

test_that(".factor_covariance works", {
  set.seed(10)
  mat <- stats::cov(matrix(sample(200), 20, 10))
  res <- .factor_covariance(mat)

  expect_true(all(dim(res$sqrt_cov) == dim(mat)))
  expect_true(all(dim(res$sqrt_inv) == dim(mat)))
})

test_that(".factor_covariance gives the correct matrix", {
  set.seed(10)
  n <- 10
  x <- rnorm(n)
  dat <- cbind(x, 2*x + rnorm(n), .5*x+0.1*rnorm(n))
  mat <- stats::cov(dat)

  res <- .factor_covariance(mat)

  expect_true(sum(abs(res$sqrt_cov %*% t(res$sqrt_cov) - mat)) < 1e-6)
  expect_true(sum(abs(res$sqrt_inv %*% t(res$sqrt_inv) - solve(mat))) < 1e-6)
})

########

## .whiten is correct

test_that(".whiten works", {
  n <- 10
  set.seed(10)
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)
  polyhedra <- binSegInf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  res <- .whiten(gaussian, polyhedra)

  expect_true(is.list(res))
  expect_true("forward_translation" %in% names(res))
  expect_true("backward_translation" %in% names(res))
  expect_true(class(res$gaussian) == "gaussian")
  expect_true(class(res$polyhedra) == "polyhedra")
})

############

## .generate_directions is correct

test_that(".generate_directions works", {
  res <- .generate_directions(10)

  expect_true(all(dim(res) == c(20,10)))
})

############

## .sampler_hit_run_line is correct

test_that(".sampler_hit_run_line works", {
  set.seed(10)
  n <- 10
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)
  polyhedra <- binSegInf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  segments <- .segments(n, binSegInf::jumps(fit), 1)

  res <- .sampler_hit_run_line(y, gaussian, segments, polyhedra, num_samp = 100)

  expect_true(all(dim(res) == c(10, 100)))
})

test_that(".sampler_hit_run_line gives the correct mean for no jumps", {
  set.seed(20)
  n <- 10
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binSegInf::binSeg_fixedSteps(y, 1)
  polyhedra <- binSegInf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  segments <- .segments(n, binSegInf::jumps(fit), 1)

  res <- .sampler_hit_run_line(y, gaussian, segments, polyhedra, num_samp = 100)

  bool_vec <- apply(res, 2, function(x){
    abs(mean(x) - mean(y)) < 1e-6
  })

  expect_true(all(bool_vec))
})
