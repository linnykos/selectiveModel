context("Test plane Rcpp")

.plane <- function(a, b = 0){
  if(!is.matrix(a)){a <- matrix(a, nrow = 1)}
  stopifnot(length(b) == nrow(a))
  stopifnot(length(which(a != 0)) > 0)

  l2_vec <- apply(a, 1, .l2norm)
  b <- as.numeric(b/l2_vec)
  if(length(l2_vec) > 1){
    a <- diag(1/l2_vec)%*%a
  } else {
    a <- a/l2_vec
  }

  structure(list(a = a, b = b), class = "plane")
}

.intersect_plane_basis <- function(plane, y, v, w, tol = 1e-6){
  a <- plane$a%*%cbind(v, w)
  if(length(which(abs(a) > tol)) > 0){
    .plane(a, plane$b - plane$a%*%y)
  } else{
    NA
  }
}

.point_on_plane <- function(plane){
  if(nrow(plane$a) == 1){
    d <- length(plane$a)
    vec <- rep(0, d)
    idx <- which(plane$a != 0)
    stopifnot(length(idx) >= 1)
    idx <- idx[1]
    vec[-idx] <- 1
    vec[idx] <- as.numeric(plane$b - plane$a[-idx]%*%vec[-idx])/plane$a[idx]

    vec
  } else {
    k <- nrow(plane$a); n <- ncol(plane$a)
    mat <- matrix(0, ncol = 3*n, nrow = k+2*n)

    mat[1:k,1:n] <- plane$a
    mat[1:k,(n+1):(2*n)] <- -plane$a

    diag(mat[(k+1):nrow(mat),(2*n+1):(3*n)]) <- 1
    diag(mat[(k+1):(k+n),1:n]) <- -1

    diag(mat[(k+1+n):nrow(mat),(2*n+1):(3*n)]) <- 1
    diag(mat[(k+1+n):nrow(mat),(n+1):(2*n)]) <- -1

    vec <- c(plane$b, rep(0, 2*n))
    res <- lpSolve::lp(objective.in = c(rep(0, 2*n), rep(1, n)), const.mat = mat,
                       const.dir = c(rep("=", k), rep(">=", 2*n)),
                       const.rhs = vec)

    if(res$status == 2) {
      stop("LP to find point on plane failed")
    }
    res$solution[1:n] - res$solution[(n+1):(2*n)]
  }
}

.distance_point_to_plane <- function(point, plane){
  stopifnot(length(point) == length(plane$a))
  stopifnot(nrow(plane$a) == 1)

  x <- .point_on_plane(plane)
  .l2norm(plane$a%*%(point - x))/.l2norm(plane$a)
}

#########################

test_that("Plane can be generated properly", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1)

    res1 <- .plane(a = a, b = b)
    res2 <- new(Plane, a, b)

    sum(abs(res1$a - res2$a)) <= 1e-6 & abs(res1$b - res2$b) <= 1e-6
  })

  expect_true(all(bool))
})

########

test_that("c_intersect_plane_basis works properly", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1)
    y <- rnorm(10); v <- rnorm(10); w <- rnorm(10)

    res1 <- .plane(a = a, b = b)
    res1 <- .intersect_plane_basis(res1, y, v, w)

    res2 <- new(Plane, a, b)
    res2$c_intersect_basis(y, v, w)

    sum(abs(res1$a - res2$a)) <= 1e-6 & abs(res1$b - res2$b) <= 1e-6
  })

  expect_true(all(bool))
})

############

test_that("c_point_on_plane works properly", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1)

    plane <- new(Plane, a, b)
    res <- plane$c_point_on_plane()

    sum(abs(res %*% plane$a - plane$b)) <= 1e-6
  })

  expect_true(all(bool))
})

###########

test_that("c_distance_point_to_plane works properly", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1)
    point <- rnorm(10)

    plane1 <- .plane(a = a, b = b)
    res1 <- .distance_point_to_plane(point, plane1)

    plane2 <- new(Plane, a, b)
    res2 <- plane2$c_distance_point_to_plane(point)

    abs(res1 - res2) <= 1e-6
  })

  expect_true(all(bool))
})


