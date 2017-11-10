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
#' The function always return a vector of length two. If no roots
#' are found, both values are \code{NA}. If only one root is found, the value
#' is repeated twice.
#'
#' @param a numeric
#' @param b numeric
#' @param c numeric
#'
#' @return vector of length 2
.quadratic <- function(a, b, c){

}
