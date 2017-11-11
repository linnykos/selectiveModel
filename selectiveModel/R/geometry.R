#' Construct a 2-dimensional circle class
#'
#' @param center vector of length 2
#' @param radius numeric
#'
#' @return a \code{circle} object
.circle <- function(center, radius){

  structure(list(center = center, radius = radius), class = "circle")
}

#' Construct a plane
#'
#' Plane is of form \code{t(a)%*%x = b}.
#'
#' @param a vector
#' @param b vector
#'
#' @return a \code{line} object
.plane <- function(a, b = 0){
  stopifnot(length(b) == 1)
  structure(list(a = a, b = b), class = "plane")
}

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
#'
#' @return numeric
#' @source \url{https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_plane}
.distance_point_to_plane <- function(point, plane){
  x <- .point_on_plane(plane)
  .l2norm((point - x)%*%plane$a)/.l2norm(plane$a)
}

#' Intersect circle with line
#'
#' The output is a 2x2 matrix, where the first row represents the first
#' point and the second row represents the second point.
#' If no points are found, return \code{NA}.
#'
#' @param line \code{line} object
#' @param circle \code{circle} object
#'
#' @return matrix of size 2x2 or \code{NA}
.intersect_circle_line <- function(line, circle){

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

 sort(c((-b-term)/(2*a), (-b+term)/(2*a)))
}
