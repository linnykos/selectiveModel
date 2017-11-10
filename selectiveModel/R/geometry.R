#' Construct a 2-dimensional circle class
#'
#' @param center vector of length 2
#' @param radius numeric
#'
#' @return a \code{circle} object
.circle <- function(center, radius){

}

#' Construct a 2-d line class
#'
#' Line is of form \code{t(a)%*%x = b}.
#'
#' @param a vector of length 2
#' @param b vector of length 2
#'
#' @return a \code{line} object
.line <- function(a, b){

}

#' Compute euclidean distance from point to line
#'
#' @param line \code{line} object
#' @param point vector of length 2
#'
#' @return numeric
.distance_point_to_line <- function(line, point){

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
