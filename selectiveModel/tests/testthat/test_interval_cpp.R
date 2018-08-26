context("Test interval c++")

.euclidean_to_radian <- function(circle, point, tol = 1e-3){
  stopifnot(length(point) == 2, length(circle$center) == 2)
  stopifnot(abs(circle$radius^2 - sum(circle$center^2)) < tol)
  stopifnot(abs(sum((point-circle$center)^2) - circle$radius^2) < tol)

  if(point[2] != 0){
    atan(point[1]/point[2])
  } else {
    atan(-circle$center[2]/circle$center[1])
  }
}

.circle <- function(center, radius){
  structure(list(center = as.numeric(center), radius = as.numeric(radius)), class = "circle")
}

####################################

## c_euclidean_to_radian is correct

test_that("c_euclidean_to_radian gives the correct answer", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    center <- rnorm(2)
    circle <- .circle(center = center, radius = sqrt(sum(center^2)))
    rad <- runif(1, -2*pi, 2*pi)
    point <- 2*(sin(rad)*center[1] + cos(rad)*center[2])*c(sin(rad), cos(rad))

    res1 <- .euclidean_to_radian(circle, point)
    res2 <- .c_euclidean_to_radian_tester(center, sqrt(sum(center^2)),
                                          point)

    ifelse(abs(res1 - res2) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that("c_euclidean_to_radian works", {
  point <- c(1,1+sqrt(2))
  res <- .c_euclidean_to_radian_tester(c(1,1), sqrt(2), point)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res <= pi/2)
  expect_true(res >= -pi/2)
})

test_that("c_euclidean_to_radian works with the origin", {
  point <- c(0,0)
  res <- .c_euclidean_to_radian_tester(c(1,1), sqrt(2), point)

  expect_true(abs(2*(sin(res)*1+ cos(res)*1)) < 1e-6)
})

test_that("c_euclidean_to_radian returns the correct theta", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    center <- rnorm(2)
    rad <- runif(1, -2*pi, 2*pi)
    point <- 2*(sin(rad)*center[1] + cos(rad)*center[2])*c(sin(rad), cos(rad))

    res <- .c_euclidean_to_radian_tester(center, sqrt(sum(center^2)), point)
    point2 <- 2*(sin(res)*center[1] + cos(res)*center[2])*c(sin(res), cos(res))

    ifelse(.l2norm(point - point2) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that("c_euclidean_to_radian works with c_intersect_circle_line", {
  set.seed(10)
  a <- rnorm(2); b <- 0
  center <- rnorm(2)
  radius <- sqrt(sum(center^2))
  points <- .c_intersect_circle_tester(a, b, center, radius)

  theta <- apply(points, 2, function(x){
    .c_euclidean_to_radian_tester(center, radius, x)
  })

  bool_vec <- apply(points, 2, function(x){ #check to be on plane.
    abs(a %*%x - b) < 1e-16
  })
  expect_true(all(bool_vec))

  bool_vec <- apply(points, 2, function(x){ #check to be on plane
    abs(.l2norm(x - center) - radius) < 1e-6
  })
  expect_true(all(bool_vec))

  expect_true(all(theta <= pi/2))
  expect_true(all(theta >= -pi/2))
  expect_true(is.numeric(theta))
  expect_true(length(theta) == 2)
})
