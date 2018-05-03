context("Test radian functions")

## .radius is correct

test_that(".radius works", {
  set.seed(10)
  y <- rnorm(10)
  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)

  res <- .radius(0, y, v, w)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

test_that(".radius is cyclical with period pi, as negatives", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    y <- rnorm(10)
    v <- rnorm(10); w <- rnorm(10)
    v <- v/.l2norm(v)
    w <- .projection(w, v); w <- w/.l2norm(w)
    theta <- runif(1, min = -pi/2, max = pi/2)

    res1 <- .radius(theta, y, v, w)
    res2 <- .radius(theta+pi, y, v, w)

    ifelse(abs(res1 + res2) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

################################

## .initial_theta is correct

test_that(".initial_theta works", {
  set.seed(10)
  y <- rnorm(10)
  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)

  res <- .initial_theta(y, v, w)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res <= pi/2)
  expect_true(res >= -pi/2)
})

test_that(".initial_theta gives a radius of 0", {
  set.seed(10)
  y <- rnorm(10)
  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)
  theta <- .initial_theta(y, v, w)

  radius <- .radius(theta, y, v, w)

  expect_true(abs(radius) < 1e-6)
})

test_that(".initial_theta gives the proper theta", {
  set.seed(10)
  y <- rnorm(10)
  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)
  theta <- .initial_theta(y, v, w)

  res <- .radians_to_data(theta, y, v, w)

  expect_true(sum(abs(y - res)) < 1e-6)
})

###########################

## .radians_to_data is correct

test_that(".radians_to_data works", {
  set.seed(5)
  y <- rnorm(10)
  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)

  res <- .radians_to_data(0, y, v, w)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".radians_to_data preserves the l2 norm", {
  trials <- 100
  set.seed(5)
  y <- rnorm(10)
  target <- .l2norm(y)

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    v <- rnorm(10); w <- rnorm(10)
    v <- v/.l2norm(v)
    w <- .projection(w, v); w <- w/.l2norm(w)
    theta <- runif(1, -pi/2, pi/2)
    ynew <- .radians_to_data(theta, y, v, w)

    res <- .l2norm(ynew)
    ifelse(abs(target - res) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

##########################

## .try_polyhedra is correct

test_that(".try_polyhedra works", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)

  y_mat <- matrix(rnorm(50), ncol = 5)
  res <- .try_polyhedra(y_mat, poly)

  expect_true(is.logical(res))
  expect_true(length(res) == ncol(y_mat))
})

test_that(".try_polyhedra corrects assess if y is in the polyhedra", {
  set.seed(10)
  y <- rnorm(5)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)

  y_mat <- matrix(rnorm(500), ncol = 100)
  res <- .try_polyhedra(y_mat, poly)

  bool_vec <- sapply(1:100, function(x){
    all(poly$gamma %*% y_mat[,x] >= poly$u)
  })

  bool_vec2 <- sapply(1:100, function(x){
    obj2 <- binSegInf::binSeg_fixedSteps(y_mat[,x], 1)
    all(all(binSegInf::jumps(obj) == binSegInf::jumps(obj2)),
        all(sign(binSegInf::jump_cusum(obj)) == sign(binSegInf::jump_cusum(obj2))))
  })

  expect_true(all(res == bool_vec))
  expect_true(all(res == bool_vec2))
})

test_that(".try_polyhedra works for a single vector", {
  set.seed(10)
  y <- rnorm(5)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)
  y_new <- rnorm(5)

  res <- .try_polyhedra(y_new, poly)
  expect_true(is.logical(res))
  expect_true(length(res) == 1)
})

test_that(".try_polyhedra can FALSE when fixing l2 norm", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)

  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)
  theta <- .initial_theta(y, v, w)

  y_new <- .radians_to_data(theta + pi/2, y, v, w)

  res <- .try_polyhedra(y_new, poly)
  expect_true(!res)
})

########################

## .range_theta_polyhedra is correct

test_that(".range_theta_polyhedra works", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)

  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)

  res <- .range_theta_polyhedra(y, v, w, poly)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
})

##########################

## .euclidean_to_radian is correct

test_that(".euclidean_to_radian works", {
  circle <- .circle(center = c(1,1), radius = sqrt(2))
  point <- c(1,1+sqrt(2))
  res <- .euclidean_to_radian(circle, point)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res <= pi/2)
  expect_true(res >= -pi/2)
})

test_that(".euclidean_to_radian works with the origin", {
  circle <- .circle(center = c(1,1), radius = sqrt(2))
  point <- c(0,0)
  res <- .euclidean_to_radian(circle, point)

  expect_true(abs(2*(sin(res)*1+ cos(res)*1)) < 1e-6)
})

test_that(".euclidean_to_radian returns the correct theta", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    center <- rnorm(2)
    circle <- .circle(center = center, radius = sqrt(sum(center^2)))
    rad <- runif(1, -2*pi, 2*pi)
    point <- 2*(sin(rad)*center[1] + cos(rad)*center[2])*c(sin(rad), cos(rad))

    res <- .euclidean_to_radian(circle, point)
    point2 <- 2*(sin(res)*center[1] + cos(res)*center[2])*c(sin(res), cos(res))

    ifelse(.l2norm(point - point2) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".euclidean_to_radian works with .intersect_circle_line", {
  set.seed(10)
  plane <- .plane(rnorm(2), b = 0)
  center <- rnorm(2)
  circle <- .circle(center = center, radius = sqrt(sum(center^2)))
  points <- .intersect_circle_line(plane, circle)

  theta <- apply(points, 2, function(x){
    .euclidean_to_radian(circle, x)
  })

  bool_vec <- apply(points, 2, function(x){ #check to be on plane.
    abs(plane$a %*%x - plane$b) < 1e-16
  })
  expect_true(all(bool_vec))

  bool_vec <- apply(points, 2, function(x){ #check to be on plane
    abs(.l2norm(x - circle$center) - circle$radius) < 1e-6
  })
  expect_true(all(bool_vec))

  expect_true(all(theta <= pi/2))
  expect_true(all(theta >= -pi/2))
  expect_true(is.numeric(theta))
  expect_true(length(theta) == 2)
})

#######################

## .partition_interval is correct

test_that(".partition_interval works", {
  res <- .partition_interval(c(pi/4, pi/3))

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(ncol(res) == 2)
})

test_that(".partition_interval works when endpoints are different signs, no wrap-around", {
  res <- .partition_interval(c(-pi/3, pi/4))

  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(res[1,] == c(-pi/3, 0)))
  expect_true(all(res[2,] == c(0, pi/4)))
})

test_that(".partition_interval works when endpoints are different signs, yes wrap-around", {
  res <- .partition_interval(c(pi/3, 3*pi/4))

  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(res[1,] == c(-pi/2, -pi/4)))
  expect_true(all(res[2,] == c(pi/3, pi/2)))

  res <- .partition_interval(c(-3*pi/4, -pi/3))

  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(res[1,] == c(-pi/2, -pi/3)))
  expect_true(all(res[2,] == c(pi/4, pi/2)))
})

test_that(".partition_interval works when endpoints are same signs, yes wrap-around", {
  res <- .partition_interval(c(-7*pi/6, -pi/3))

  expect_true(all(dim(res) == c(3,2)))
  expect_true(sum(abs(res[1,] - c(-pi/2, -pi/3))) < 1e-6)
  expect_true(sum(abs(res[2,] - c(-pi/6, 0))) < 1e-6)
  expect_true(sum(abs(res[3,] - c(0, pi/2))) < 1e-6)

  res <- .partition_interval(c(pi/3, 7*pi/6))
  expect_true(all(dim(res) == c(3,2)))
  expect_true(sum(abs(res[1,] - c(-pi/2, 0))) < 1e-6)
  expect_true(sum(abs(res[2,] - c(0, pi/6))) < 1e-6)
  expect_true(sum(abs(res[3,] - c(pi/3, pi/2))) < 1e-6)
})

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
