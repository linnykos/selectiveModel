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

##########################

## .theta_seq is correct

test_that(".theta_seq works", {
  res <- .theta_seq(0, 3)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 2^3)
})

test_that(".theta_seq forms the correct filtration-like sequence", {
  res1 <- .theta_seq(0, 1)
  res2 <- .theta_seq(0, 2)
  res3 <- .theta_seq(0, 3)

  expect_true(length(res1) + length(res2) + length(res3) ==
                length(unique(c(res1, res2, res3))))

  target <- seq(-pi/2, pi/2, length.out = 2^4+1)
  target <- target[-c(1, length(target))]
  target <- target[-which(target == 0)]
  expect_true(sum(abs(sort(c(res1, res2, res3)) - target)) < 1e-6)
})

##########################

## .binary_search is correct

