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

.quadratic <- function(a, b, c, tol = 1e-8){
  stopifnot(is.numeric(a), is.numeric(b), is.numeric(c))
  stopifnot(length(c(a,b,c)) == 3)

  term <- b^2 - 4*a*c
  if(term < 0) return(NA)
  if(abs(term) < tol) return(-b/(2*a))

  sort(c((-b-sqrt(term))/(2*a), (-b+sqrt(term))/(2*a)))
}

.intersect_circle_line <- function(plane, circle, tol = 1e-6, tol2 = 1e-9,
                                   full = F){
  stopifnot(length(plane$a) == 2, length(which(plane$a != 0)) > 0)
  stopifnot(nrow(plane$a) == 1)

  dis <- .distance_point_to_plane(circle$center, plane)
  if(dis > circle$radius + tol) {
    if(full) {
      return(matrix(NA, 2, 2))
    } else {
      return(NA)
    }
  }

  if(abs(plane$a[2]) < tol){
    #treat plane$a[2] as zero
    x <- plane$b/plane$a[1]
    y <- circle$center[2] + c(-1, 1) * sqrt(circle$radius^2 - (x - circle$center[1])^2)
    mat <- rbind(x, y)

  } else if(abs(plane$a[1]) < tol) {
    #treat plane$a[1] as zero
    y <- plane$b/plane$a[2]
    x <- circle$center[1] + c(-1, 1) * sqrt(circle$radius^2 - (y - circle$center[2])^2)
    mat <- rbind(x, y)

  } else {
    a1 <- plane$a[1]; a2 <- plane$a[2]
    c1 <- circle$center[1]; c2 <- circle$center[2]

    a <- 1 + (a1/a2)^2
    b <- -2*(a1/a2)*(plane$b/a2 - c2) -2*c1
    c <- -circle$radius^2 + (plane$b/a2 - c2)^2 + c1^2

    x <- .quadratic(a, b, c)
    stopifnot(all(!is.na(x)))
    y <- (plane$b - a1*x)/a2

    if((length(x) == 1 || abs(x[1]-x[2]) < tol2) && (length(y) == 1 || abs(y[1]-y[2]) < tol2)){
      mat <- matrix(c(x[1], y[1]), nrow = 2)
      if(full) mat <- rbind(mat, NA)
    } else {
      mat <- rbind(x, y)
    }
  }

  colnames(mat) <- NULL
  rownames(mat) <- NULL

  mat
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

#################

test_that("c_intersect_circle works properly", {
  trials <- 100

  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(2); b <- 1
    plane <- .plane(a, b)
    center <- rnorm(2); radius <- 5
    circle <- .circle(center, radius)
    dis <- .distance_point_to_plane(center, plane)

    if(dis >= radius){
      return(TRUE)
    } else {
      res1 <- .intersect_circle_line(plane, circle, full = T)

      res2 <- .c_intersect_circle_tester(a, b, center, radius)

      sum(abs(res1 - res2)) < 1e-6
    }
  })

  expect_true(all(bool))
})

test_that("c_intersect_circle works", {
  res <- .c_intersect_circle_tester(c(0,1), 0, c(0,0), 1)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(2,2)))
})

test_that("c_intersect_circle gives a proper point", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- rnorm(2)
    points <- .c_intersect_circle_tester(a, 0, c(0,0), 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - 1)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that("c_intersect_circle gives a proper point with b", {
  trials <- 100
  radius <- 100
  b <- 3
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- rnorm(2)
    points <- .c_intersect_circle_tester(a, b, c(0,0), radius)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - radius)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x - b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that("c_intersect_circle gives a proper point with small values of a[1]", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- c(0, rnorm(1))
    points <- .c_intersect_circle_tester(a, b, c(0,0), radius*b)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - b*radius)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x - b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives two points with small values of a[1] always gives 2 points", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a_vec <- c(0, rnorm(1))
    points <- .c_intersect_circle_tester(a_vec, b, c(0,0), radius*b)

    ncol(points) == 2
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with small values of a[1] and offcenter circle", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- c(0, rnorm(1))
    center <- c(1,2)
    points <- .c_intersect_circle_tester(a, b, center, radius*b)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt((x[1]-center[1])^2+(x[2]-center[2])^2) - radius*b)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x - b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with small values of a[2] and offcenter circle", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- c(rnorm(1), 0)
    center <- c(1,2)
    points <- .c_intersect_circle_tester(a, b, center, radius*b)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt((x[1]-center[1])^2+(x[2]-center[2])^2) - radius*b)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x - b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line can give NA", {
  res <- .c_intersect_circle_tester(c(1, -1), 10, c(0,0), 1)

  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(is.na(res)))
})

test_that(".intersect_circle_line can give one point", {
  res <- .c_intersect_circle_tester(c(1, 0), 1, c(0,0), 1)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(is.na(res[,2])))
  expect_true(all(res[,1] == c(1,0)))
})


