#' Construct a 2-dimensional circle class
#'
#' @param center vector of length 2
#' @param radius numeric
#'
#' @return a \code{circle} object
.circle <- function(center, radius){
  structure(list(center = as.numeric(center), radius = as.numeric(radius)), class = "circle")
}

#' Construct a line object
#'
#' @param point vector
#' @param direction vector
#'
#' @return a \code{line} object
.line <- function(point, direction){
  stopifnot(length(point) == length(direction))

  direction <- direction/.l2norm(direction)
  structure(list(point = point, direction = direction), class = "line")
}

#' Construct a hyperplane
#'
#' Hyperplane is of form \code{a%*%x = b}.
#' If \code{a} is a matrix, this object represents an intersection of hyperplanes.
#'
#' @param a matrix
#' @param b vector
#'
#' @return a \code{plane} object
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

#' Convert a higher-dimensional plane into the two-dimensional coordinate system
#'
#' If the coordinate system set up by \code{v} and \code{w} is parallel to \code{plane},
#' an \code{NA} is returned.
#'
#' @param plane \code{plane} object
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#' @param tol small positive number
#'
#' @return either a \code{plane} object or \code{NA}
.intersect_plane_basis <- function(plane, y, v, w, tol = 1e-6){
  a <- plane$a%*%cbind(v, w)
  if(length(which(abs(a) > tol)) > 0){
    .plane(a, plane$b - plane$a%*%y)
  } else{
    NA
  }
}

#' Returns a point on the plane
#'
#' If \code{plane} is only one plane, this function uses a deterministic
#' method to find a point that lines along a hyperplane.
#'
#' If \code{plane} is an intersection of many hyperplanes, this function
#' uses a linear program to find a point on the plane. The linear program
#' is minimize sum t_1,...t_n such that Ay = b, t_i >= y_i, t_i >= -y_i for all
#' i in 1 to n.
#'
#' @param plane \code{plane} object
#'
#' @return a vector
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

#' Compute euclidean distance from point to plane
#'
#' @param point vector
#' @param plane \code{plane} object
#'
#' @return numeric
#' @source \url{https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_plane}
.distance_point_to_plane <- function(point, plane){
  stopifnot(length(point) == length(plane$a))
  stopifnot(nrow(plane$a) == 1)

  x <- .point_on_plane(plane)
  .l2norm(plane$a%*%(point - x))/.l2norm(plane$a)
}

#' Intersect circle with plane
#'
#' The \code{plane} object must be in 2-dimensions
#'
#' The output is a 2x2 matrix, where the first row represents the first
#' point and the second row represents the second point.
#'
#' If there is only one point, then the output is a 2x1 matrix.
#' If no points are found, return \code{NA}.
#'
#' \code{tol} is used to differentiate zero from non-zero. \code{tol2} is used
#' to determine if there are one or two roots.
#'
#' @param plane \code{plane} object
#' @param circle \code{circle} object
#' @param tol small positive number
#'
#' @return matrix of size 2x2 or \code{NA}
.intersect_circle_line <- function(plane, circle, tol = 1e-6){
  stopifnot(length(plane$a) == 2, length(which(plane$a != 0)) > 0)
  stopifnot(nrow(plane$a) == 1)

  dis <- .distance_point_to_plane(circle$center, plane)
  if(dis > circle$radius + tol) return(NA)

  if(abs(plane$a[1]) > tol){
    a1 <- plane$a[1]; a2 <- plane$a[2]
  } else {
    a1 <- plane$a[2]; a2 <- plane$a[1]
  }

  a <- 1 + (a2/a1)^2
  b <- -2*(a2/a1)*(plane$b/a1 - circle$center[1]) -
    2*circle$center[2]
  c <- -circle$radius^2 + (plane$b/a1 - circle$center[1])^2 +
    circle$center[2]^2

  y <- .quadratic(a, b, c)
  stopifnot(all(!is.na(y)))
  x <- (plane$b -a2*y)/a1

  if(length(y) == 1 || abs(y[1]-y[2]) < tol){
    mat <- matrix(c(x[1], y[1]), ncol = 2)
  } else {
    mat <- cbind(x, y)
  }
  colnames(mat) <- NULL

  if(abs(plane$a[1]) < tol) mat <- mat[,c(2,1)]
  mat
}

#' Determine a point on plane closest to origin
#'
#' @param plane \code{plane} object
#'
#' @return vector
#' @source \url{https://math.stackexchange.com/questions/832279/closest-point-to-a-vector-in-a-subspace}
.closest_point_to_origin <- function(plane){
  rowspace <- .sample_matrix_space(plane$a, null = F)
  point <- .point_on_plane(plane)

  #project against each vector
  mat <- sapply(1:ncol(rowspace), function(x){
    rowspace[,x]%*%t(rowspace[,x])%*%point
  })

  rowSums(mat)
}

#' Sample from a n-1 dimensional unit sphere uniformally
#'
#' @param n positive integer larger than 1
#'
#' @return vector of length n
.sample_sphere <- function(n){
  vec <- stats::rnorm(n)
  vec/.l2norm(vec)
}

#' Quadratic formula
#'
#' The function returns a vector of length two if there are two roots,
#' the first smaller than the second. If there is only one root, then only
#' one numeric is returned. If no roots are found, then \code{NA} is returned.
#'
#' @param a numeric
#' @param b numeric
#' @param c numeric
#' @param tol small positive number
#'
#' @return vector of length 2
.quadratic <- function(a, b, c, tol = 1e-4){
  stopifnot(is.numeric(a), is.numeric(b), is.numeric(c))
  stopifnot(length(c(a,b,c)) == 3)

 term <- b^2 - 4*a*c
 if(term < 0) return(NA)
 if(abs(term) < tol) return(-b/(2*a))

 sort(c((-b-sqrt(term))/(2*a), (-b+sqrt(term))/(2*a)))
}

#' Intersect a polyhedron with a line
#'
#' Computes one endpoint by solving the following optimization problem:
#' max a s.t. polyhedra$gamma * (point + a*direction) >= polyhedra$u.
#' Since the linear program forces all variables to be non-negative, additional
#' transformations are needed. The other endpoint is solved by a minimization problem.
#'
#' This function assumes that \code{line$point} is a valid point inside the
#' polyhedra.
#'
#' @param polyhedra a \code{polyhedra} object
#' @param line a \code{line} object
#'
#' @return a vector of length 2
.intersect_polyhedron_line <- function(polyhedra, line){
  stopifnot(ncol(polyhedra$gamma) == length(line$direction))
  stopifnot(all(polyhedra$gamma %*% line$point >= polyhedra$u))
  d <- length(polyhedra$u)

  vec <- polyhedra$gamma %*% line$direction
  rhs_vec <- polyhedra$u - polyhedra$gamma %*% line$point
  constr_mat <- cbind(vec, -vec)

  res1 <- lpSolve::lp(direction = "min",
                      objective.in = c(1, -1), const.mat = constr_mat,
                      const.dir = rep(">=", d), const.rhs = rhs_vec)

  res2 <- lpSolve::lp(direction = "max",
                      objective.in = c(1, -1), const.mat = constr_mat,
                      const.dir = rep(">=", d), const.rhs = rhs_vec)

  val1 <- ifelse(res1$status == 0, res1$solution[1] - res1$solution[2],
                 -Inf)
  val2 <- ifelse(res2$status == 0, res2$solution[1] - res2$solution[2],
                 Inf)

  stopifnot(0 >= val1, 0 <= val2)

  c(val1, val2)
}

.transform_point <- function(val, shift, rotation){

}
