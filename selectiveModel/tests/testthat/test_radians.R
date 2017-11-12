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
  expect_true(!is.matrix(res))
})

test_that(".range_theta_polyhedra returns correct theta's", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)

  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)

  res <- .range_theta_polyhedra(y, v, w, poly)

  bool_vec <- sapply(res, function(x){
    .try_polyhedra(.radians_to_data(x, y, v, w), poly)
  })

  expect_true(all(bool_vec))
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
    rad <- runif(1,-pi/2, pi/2)
    point <- 2*(sin(rad)*center[1] + cos(rad)*center[2])*c(sin(rad), cos(rad))

    res <- .euclidean_to_radian(circle, point)

    ifelse(abs(res - rad) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})
