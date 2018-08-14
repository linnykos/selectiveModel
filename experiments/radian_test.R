
#######################

## .consecutive_true is correct

test_that(".consecutive_true works", {
  vec <- c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
  res <- .consecutive_true(vec)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(ncol(res) == 2)
})

test_that(".consecutive_true is correct", {
  vec <- c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
  res <- .consecutive_true(vec)

  expect_true(all(dim(res) == c(3,2)))
  expect_true(all(t(res) == c(1,3,6,7,9,10)))
})

test_that(".consecutive_true also works when it does not start with TRUE", {
  vec <- c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
  res <- .consecutive_true(vec)

  expect_true(all(dim(res) == c(3,2)))
  expect_true(all(t(res) == c(2,3,6,7,9,10)))
})

###############################

## .intersect_two_intervals is correct

test_that(".intersect_two_intervals works", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/3, pi/3))

  res <- .intersect_two_intervals(mat1, mat2)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(ncol(res) == 2)
})

test_that(".intersect_two_intervals is correct", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/3, pi/3))

  res <- .intersect_two_intervals(mat1, mat2)

  expect_true(sum(abs(as.numeric(res) - c(-pi/6, pi/3))) < 1e-6)
})

test_that(".intersect_two_intervals can output two intervals", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/2, 0))

  res <- .intersect_two_intervals(mat1, mat2)

  expect_true(nrow(res) == 2)
  expect_true(sum(abs(res[1,] - c(-pi/2, -pi/3))) < 1e-6)
  expect_true(sum(abs(res[2,] - c(-pi/6, 0))) < 1e-6)
})

test_that(".intersect_two_intervals can be chained together", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/2, 0))
  mat3 <- .partition_interval(c(-3*pi/4, -pi/4))

  res <- .intersect_two_intervals(mat1, mat2)
  res <- .intersect_two_intervals(res, mat3)

  expect_true(nrow(res) == 1)
  expect_true(sum(abs(res[1,] - c(-pi/2, -pi/3))) < 1e-6)
})

test_that(".intersect_two_intervals errors when there is no intersection", {
  mat1 <- .partition_interval(c(-pi/4, 0))
  mat2 <- .partition_interval(c(pi/4, pi/2))

  expect_error(.intersect_two_intervals(mat1, mat2))
})

test_that(".intersect_two_intervals can intersect the entire space", {
  mat1 <- .partition_interval(c(-pi/4, 0))
  mat2 <- .partition_interval(c(-pi/2, pi/2))

  res <- .intersect_two_intervals(mat1, mat2)

  expect_true(sum(abs(as.numeric(mat1) - as.numeric(res))) < 1e-6)
})

test_that(".intersect_two_intervals works for many test cases", {
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

    res <- .intersect_two_intervals(mat1, mat2)

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

#########################

## .intersect_intervals is correct

test_that(".intersect_intervals works", {
  mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
  mat2 <- .partition_interval(c(-pi/2, 0))
  mat3 <- .partition_interval(c(-3*pi/4, -pi/4))

  lis <- list(mat1, mat2, mat3)
  res <- .intersect_intervals(lis)

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
  res2 <- .intersect_intervals(lis)

  expect_true(sum(abs(as.numeric(res) - as.numeric(res2))) < 1e-6)
})
