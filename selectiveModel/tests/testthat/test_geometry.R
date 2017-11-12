context("Test geometry")

## .quadratic is correct

test_that(".quadratic works", {
  res <- .quadratic(-1, 2, 5)

  expect_true(length(res) == 2)
  expect_true(res[1] < res[2])
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(all(!is.na(res)))
})

test_that(".quadratic can return only one number", {
  res <- .quadratic(1, 2, 1)

  expect_true(length(res) == 1)
  expect_true(!is.na(res))
})

test_that(".quadratic can return only NA", {
  res <- .quadratic(1, 2, 5)

  expect_true(length(res) == 1)
  expect_true(is.na(res))
})

#################################

## .point_on_plane is correct

test_that(".point_on_plane works", {
  obj <- .plane(c(1,2,3,4), 4)
  res <- .point_on_plane(obj)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 4)
})

test_that(".point_on_plane returns a point on the plane", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    obj <- .plane(rnorm(10), rnorm(1))
    res <- .point_on_plane(obj)

    ifelse(abs(res%*%obj$a - obj$b) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".point_on_plane handles cases parallel to axes", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- rep(0, 10)
    idx <- sample(1:10, 5)
    vec[idx] <- rnorm(5)
    obj <- .plane(vec, rnorm(1))
    res <- .point_on_plane(obj)

    ifelse(abs(res%*%obj$a - obj$b) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".point_on_plane handles cases where parallel to only but one axes", {
  obj <- .plane(c(0,0,1), 0)
  res <- .point_on_plane(obj)

  expect_true(abs(res%*%obj$a - obj$b) < 1e-6)
})

##########################

## .distance_point_to_plane is correct

test_that(".distance_point_to_plane works", {
  plane <- .plane(c(0,0,1), 0)
  point <- c(1,0,1)
  res <- .distance_point_to_plane(point, plane)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res >= 0)
})

test_that(".distance_point_to_plane is always non-negative", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    plane <- .plane(rnorm(10), rnorm(1))
    point <- rnorm(10)
    res <- .distance_point_to_plane(point, plane)

    ifelse(res >= 0, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".distance_point_to_plane satisifies triangle inequality", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    plane <- .plane(rnorm(10), rnorm(1))
    point1 <- rnorm(10)
    point2 <- rnorm(10)
    res1 <- .distance_point_to_plane(point1, plane)
    res2 <- .distance_point_to_plane(point2, plane)
    dist <- .l2norm(point1 - point2)

    ifelse(abs(res2-res1) <= dist, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

####################

## .intersect_circle_line is correct

test_that(".intersect_circle_line works", {
  plane <- .plane(c(0,1), 0)
  circle <- .circle(c(0,0), 1)
  res <- .intersect_circle_line(plane, circle)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(2,2)))
})

test_that(".intersect_circle_line gives a proper point", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    plane <- .plane(rnorm(2), b = 0)
    circle <- .circle(c(0,0), 1)
    points <- .intersect_circle_line(plane, circle)

    #check to see points are in circle
    res1 <- apply(points, 1, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - 1)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 1, function(x){abs(plane$a%*%x) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line can give NA", {
  plane <- .plane(c(1,-1), b = 10)
  circle <- .circle(c(0,0), 1)
  res <- .intersect_circle_line(plane, circle)

  expect_true(length(res) == 1)
  expect_true(is.na(res))
})

test_that(".intersect_circle_line can give one point", {
  plane <- .plane(c(1,0), b = 1)
  circle <- .circle(c(0,0), 1)
  res <- .intersect_circle_line(plane, circle)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(1,2)))
})
