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
.range_theta_polyhedra <- function(y, v, w, polyhedra, tol = 1e-2){

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
  -2*as.numeric(t(y)%*%(v*sin(theta)+w*cos(theta)))
}

#' Compute the arctan to find the initial theta
#'
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#'
#' @return radians between -pi/2 and pi/2.
.initial_theta <- function(y, v, w){
  atan(-as.numeric(t(y)%*%w)/as.numeric(t(y)%*%v))
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
  y + .radius(theta, y, v, w)*(sin(theta)*v + cos(theta)*w)
}

#' Determine if each y is in polyhedra
#'
#' @param y_mat matrix where each column represents a different sample, or
#' a vector
#' @param polyhedra \code{polyhedra} object
#'
#' @return vector of booleans, one for each column of \code{y_mat} (or just one)
.try_polyhedra <- function(y_mat, polyhedra){
  if(is.matrix(y_mat)){
    stopifnot(ncol(polyhedra$gamma) == nrow(y_mat))
  } else {
    stopifnot(ncol(polyhedra$gamma) == length(y_mat))
  }

  res <- polyhedra$gamma %*% y_mat

  if(is.matrix(y_mat)){
    sapply(1:ncol(res), function(x){all(res[,x] >= polyhedra$u)})
  } else {
    all(res >= polyhedra$u)
  }
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
.binary_search <- function(theta_start, theta_end, y, v, w, polyhedra, tol = 1e-2){

}

#' Creates a vector of theta's
#'
#' Creates a vector of thetas to try starting at \code{theta-pi/2}
#' (includsive) to \code{theta+pi/2} (exclusive) based on \code{generation}.
#'
#' @param theta radians
#' @param generation positive integer
#'
#' @return vector of radians
.theta_seq <- function(theta = 0, generation = 1){
  vec <- seq(theta-pi/2, theta+pi/2, length.out = 2^(generation+1)+1)
  n <- length(vec)
  vec[seq(2, n, by = 2)]
}


