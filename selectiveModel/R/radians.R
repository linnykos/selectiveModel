#' Find the range of theta's that are in the polyhedra
#'
#' Brute-force finds theta's such that for each theta,
#' \code{.radians_to_data(theta, y, v, w)} is inside the polyhedra.
#'
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#' @param polyhedra \code{polyhedra} object
#' @param attempts positive integer
#'
#' @return vector of possible theta's
.range_theta_polyhedra <- function(y, v, w, polyhedra, attempts = 10){
  stopifnot(.try_polyhedra(y, polyhedra))
  stopifnot(abs(.l2norm(v) - 1) < 1e-4,  abs(.l2norm(w) - 1) < 1e-4, abs(v%*%w) < 1e-4)

  #create possible theta's
  theta <- .initial_theta(y, v, w)
  seq_vec <- seq(theta-pi/2, theta+pi/2, length.out = attempts + 1)
  seq_vec <- sort(unique(c(seq_vec[-length(seq_vec)], theta)))

  #try each theta
  y_mat <- sapply(seq_vec, function(x){
    .radians_to_data(x, y, v, w)
  })
  bool_vec <- .try_polyhedra(y_mat, polyhedra)

  stopifnot(sum(bool_vec) >= 1)
  seq_vec[which(bool_vec)]
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


