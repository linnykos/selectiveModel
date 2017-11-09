#' Find the range of theta's that are in the polyhedra
#'
#' Find an interval (vector of length 2) of radians such that for any
#' theta inside this interval,
#' \code{.radians_to_data(theta_start, y, v, w)} is inside the polyhedra
#' dictated by \code{polyhedra}.
#'
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#' @param polyhedra \code{polyhedra} object
#' @param tol small positive number
#'
#' @return vector of 2 radians, the first being smaller than the other
.range_theta_polyhedra <- function(y, v, w, polyhedra, tol = 1e-4){

}

#' Radius function
#'
#' @param theta radians
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#'
#' @return numeric
.radius <- function(theta, y, v, w){

}

#' Compute the arctan to find the initial theta
#'
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#'
#' @return radians between -pi/2 and pi/2.
.initial_theta <- function(y, v, w){

}

#' Convert radians into a sample
#'
#' @param theta radians
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#'
#' @return vector with length of \code{length(y)}
.radians_to_data <- function(theta, y, v, w){

}

#' Determine if each y is in polyhedra
#'
#' @param y_mat matrix where each column represents a different sample
#' @param polyhedra \code{polyhedra} object
#'
#' @return vector of booleans, one for each column of \code{y_mat}
.try_polyhedra <- function(y_mat, polyhedra){

}

#' Perform binary search to determine cutoff point
#'
#' \code{theta_start} is defined to be the endpoint where
#' \code{.radians_to_data(theta_start, y, v, w)} is inside the polyhedra
#' dictated by \code{polyhedra}. Likewise, \code{theta_end} is defined to be
#' endpoint where \code{.radians_to_data(theta_end, y, v, w)} is not.
#'
#' @param theta_start radians endpoint of search space
#' @param theta_end radians endpoint of search space
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#' @param polyhedra \code{polyhedra} object
#' @param tol small positive number
#'
#' @return radians
.binary_search <- function(theta_start, theta_end, y, v, w, polyhedra, tol = 1e-4){

}

#' Creates a vector of theta's
#'
#' Creates a vector of thetas to try starting at \code{theta_start}
#' (includsive) to \code{theta_end} (exclusive) based on \code{generation}.
#'
#' @param theta_start radians
#' @param theta_end radians
#' @param generation positive integer
#'
#' @return vector of radians
.theta_seq <- function(theta_start, theta_end, generation){

}


