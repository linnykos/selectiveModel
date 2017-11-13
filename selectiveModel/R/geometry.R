#' Construct a 2-dimensional circle class
#'
#' @param center vector of length 2
#' @param radius numeric
#'
#' @return a \code{circle} object
.circle <- function(center, radius){

  structure(list(center = as.numeric(center), radius = as.numeric(radius)), class = "circle")
}

#' Construct a hyperplane
#'
#' Hyperplane is of form \code{t(a)%*%x = b}.
#'
#' @param a vector
#' @param b vector
#'
#' @return a \code{plane} object
.plane <- function(a, b = 0){
  stopifnot(length(b) == 1)
  stopifnot(length(which(a != 0)) > 0)

  b <- as.numeric(b/.l2norm(a))
  a <- as.numeric(a/.l2norm(a))

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
#' Uses a deterministic method to find a point that lines along a hyperplane
#'
#' @param plane \code{plane} object
#'
#' @return a vector
.point_on_plane <- function(plane){
  d <- length(plane$a)
  vec <- rep(0, d)
  idx <- which(plane$a != 0)
  stopifnot(length(idx) >= 1)
  idx <- idx[1]
  vec[-idx] <- 1
  vec[idx] <- as.numeric(plane$b - plane$a[-idx]%*%vec[-idx])/plane$a[idx]

  vec
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

  x <- .point_on_plane(plane)
  .l2norm((point - x)%*%plane$a)/.l2norm(plane$a)
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
#' @param plane \code{plane} object
#' @param circle \code{circle} object
#' @param tol small positive number
#'
#' @return matrix of size 2x2 or \code{NA}
.intersect_circle_line <- function(plane, circle, tol = 1e-6){
  stopifnot(length(plane$a) == 2, length(which(plane$a != 0)) > 0)

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
