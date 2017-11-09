.radius <- function(theta, y, v, w){

}

#' Compute the arctan to find the initial theta
#'
#' @param y vector
#' @param v vector
#' @param w vector
#'
#' @return radians between -pi/2 and pi/2.
.initial_theta <- function(y, v, w){

}

.radians_to_data <- function(y, v, w, rad){

}

.try_polyhedra <- function(y_mat, polyhedra){

}

.binary_search <- function(theta_start, theta_end, y, v, w, polyhedra, tol = 1e-4){

}

#' Creates a vector of theta's
#'
#' Creates a vector of thetas to try starting at \code{theta_start}
#' (includsive) to \code{theta_end} (exclusive) based on \code{generation}.
#'
#' @param theta_start numeric
#' @param theta_end numeric
#' @param generation numeric
#'
#' @return
.theta_seq <- function(theta_start, theta_end, generation){

}

.range_theta_polyhedra <- function(y, v, w, polyhedra, tol = 1e-4){

}
