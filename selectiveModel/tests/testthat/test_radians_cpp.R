context("Test radians C++")

## intersect_intervals_simult is correct

.intersect_intervals <- function(lis){
  Reduce(.intersect_two_intervals, lis)
}

.intersect_two_intervals <- function(mat1, mat2){
  vec <- sort(unique(c(as.numeric(mat1), as.numeric(mat2))))
  vec <- sort(c(vec, .interpolate(vec)))

  bool_vec <- sapply(vec, function(x){
    all(.theta_in_interval(x, mat1), .theta_in_interval(x, mat2))
  })

  idx_mat <- .consecutive_true(bool_vec)

  mat <- matrix(vec[idx_mat], ncol = 2)
  mat <- mat[order(mat[,1]),,drop = F]

  stopifnot(all(as.numeric(t(mat)) == sort(as.numeric(t(mat)))))
  mat
}

.interpolate <- function(vec){
  n <- length(vec)
  sapply(2:n, function(x){
    mean(vec[(x-1):x])
  })
}

.theta_in_interval <- function(theta, mat){
  stopifnot(abs(theta) <= pi/2)

  vec <- any(apply(mat, 1, function(x){
    x[1] <= theta & theta <= x[2]
  }))
}

.consecutive_true <- function(vec){
  idx <- which(vec)
  if(length(idx) == 0) stop("No intersection")
  breakpoint <- which(sapply(2:length(idx), function(x){idx[x]-idx[x-1] != 1}))

  breakpoint <- c(0, breakpoint, length(idx))
  mat <- t(sapply(2:length(breakpoint), function(x){
    c(idx[breakpoint[x-1]+1], idx[breakpoint[x]])
  }))

  #remove singletons
  idx <- which(mat[,1] != mat[,2])
  if(length(idx) == 0) stop("No intersection")
  mat[idx,]
}

.generate_random_list <- function(i){
  set.seed(10*i)
  num_elements <- ceiling(runif(1, 1, 200))
  common_interval <- sort(runif(2, -pi/2, pi/2))
  mid_point <- mean(common_interval)
  lapply(1:num_elements, function(x){
    .interval(range(c(runif(2, -pi/2 + 0.1, pi/2), common_interval)),
              mid_point)
  })
}

test_that("intersect_intervals_simult is correct wrt R implementation", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    lis <- .generate_random_list(x)

    res1 <- .intersect_intervals(lis)
    res2 <- intersect_intervals_simult(lis)

    all.equal(res1, res2)
  })

  expect_true(all(bool))
})


test_that("intersect_intervals_simult works", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/2, 0))
  mat3 <- .partition_interval(c(-3*pi/4, -pi/4))

  lis <- list(mat1, mat2, mat3)
  res <- intersect_intervals_simult(lis)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(ncol(res) == 2)
})

test_that(".intersect_intervals gives the correct output", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/2, 0))
  mat3 <- .partition_interval(c(-3*pi/4, -pi/4))

  res <- .intersect_two_intervals(mat1, mat2)
  res <- .intersect_two_intervals(res, mat3)

  lis <- list(mat1, mat2, mat3)
  res2 <- intersect_intervals_simult(lis)

  expect_true(sum(abs(as.numeric(res) - as.numeric(res2))) < 1e-6)
})

test_that("intersect_intervals_simult is correct", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/3, pi/3))

  res <- intersect_intervals_simult(list(mat1, mat2))

  expect_true(sum(abs(as.numeric(res) - c(-pi/6, pi/3))) < 1e-6)
})

test_that("intersect_intervals_simult can output two intervals", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/2, 0))

  res <- intersect_intervals_simult(list(mat1, mat2))

  expect_true(nrow(res) == 2)
  expect_true(sum(abs(res[1,] - c(-pi/2, -pi/3))) < 1e-6)
  expect_true(sum(abs(res[2,] - c(-pi/6, 0))) < 1e-6)
})

test_that("intersect_intervals_simult can be chained together", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/2, 0))
  mat3 <- .partition_interval(c(-3*pi/4, -pi/4))

  res <- intersect_intervals_simult(list(mat1, mat2))
  res <- intersect_intervals_simult(list(res, mat3))

  expect_true(nrow(res) == 1)
  expect_true(sum(abs(res[1,] - c(-pi/2, -pi/3))) < 1e-6)
})

test_that("intersect_intervals_simult returns empty matrix when there is no intersection", {
  mat1 <- .partition_interval(c(-pi/4, 0))
  mat2 <- .partition_interval(c(pi/4, pi/2))

  res <- intersect_intervals_simult(list(mat1, mat2))
  expect_true(nrow(res) == 0)
})

test_that("intersect_intervals_simult can intersect the entire space", {
  mat1 <- .partition_interval(c(-pi/4, 0))
  mat2 <- .partition_interval(c(-pi/2, pi/2))

  res <- intersect_intervals_simult(list(mat1, mat2))

  expect_true(sum(abs(as.numeric(mat1) - as.numeric(res))) < 1e-6)
})

test_that("intersect_intervals_simult works for many test cases", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    val1 <- runif(1, -3*pi/2, 3*pi/2)
    val2 <- runif(1, max(-3*pi/2, val1 - pi), min(3*pi/2, val1 + pi))
    vec1 <- sort(c(val1, val2))

    val3 <- runif(1, vec1[1], vec1[2])
    val4 <- runif(1, max(-3*pi/2, val3 - pi), min(3*pi/2, val3 + pi))

    mat1 <- .partition_interval(vec1)
    mat2 <- .partition_interval(sort(c(val3, val4)))

    res <- intersect_intervals_simult(list(mat1, mat2))

    seq_vec <- seq(-pi/2, pi/2, length.out = 20)
    bool1 <- sapply(seq_vec, function(x){
      all(.theta_in_interval(x, mat1), .theta_in_interval(x, mat2))
    })
    bool2 <- sapply(seq_vec, function(x){
      .theta_in_interval(x, res)
    })

    all(bool1 == bool2)
  })

  expect_true(all(bool_vec))
})


