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
  theta_vec <- seq(theta-pi/2, theta+pi/2, length.out = attempts + 1)
  theta_vec <- sort(unique(c(theta_vec[-length(theta_vec)], theta)))

  #try each theta
  y_mat <- sapply(theta_vec, function(x){
    .radians_to_data(x, y, v, w)
  })
  bool_vec <- .try_polyhedra(y_mat, polyhedra)

  stopifnot(sum(bool_vec) >= 1)
  theta_vec[which(bool_vec)]
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
  stopifnot(as.numeric(t(y)%*%v) != 0)
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

#' Construct the interval of radians
#'
#' Output is always a vector of length four. We are representing the
#' interval as two distinct intervals. The first four numbers represent
#' the interval as a subset of \code{c(-pi/2,0)}. The last four numbers
#' represent the interval as a subset of \code{c(0, pi/2)}.
#'
#' Each of the four numbers are two pairs, and all four numbers are
#' in increasing order. For example, if the interval \code{c(-pi/2, -pi/3)}
#' and \code{c(-pi/4, 0)} were to be expressed, the four numbers would be
#' \code{c(-pi/2, -pi/3, -pi/4, 0)}. If the interval \code{c(-pi/2, -pi/3)}
#' were to be expressed, the four numbers would be \code{c(-pi/2, -pi/3, NA, NA)}.
#' If the interval \code{c(-pi/2, 0)} were to be expressed, the interval would be
#' \code{c(-pi/2, -pi/4, -pi/4, 0)} by default.
#'
#' In short, for every four numbers, the first is always -pi/2, 0, or \code{NA},
#' while the last is always 0, pi/2, or \code{NA}.
#'
#' If an interval does not intersect with (say)  \code{c(-pi/2,0)}, four \code{NA}s
#' are placed there instead.
#'
#' \code{theta} is a radian that is guaranteed to be "in" the interval.
#'
#' @param endpoints vector of length 2
#' @param theta numeric
#'
#' @return vector of length 8
.interval <- function(endpoints, theta){
  stopifnot(all(c(endpoints, theta) <= pi/2), all(c(endpoints, theta) >= -pi/2))

  interval <- .basic_interval(endpoints, theta)
  interval <- .partition_interval(interval)
}

.basic_interval <- function(endpoints, theta){
  endpoints <- sort(endpoints)
  if(endpoints[1] <= theta & theta <= endpoints[2]){
    endpoints
  } else {
    c(endpoints[2], endpoints[1]+pi)
  }
}

#' Partition an interval of radians
#'
#' This function depends heavily on the representation used in
#' \code{.basic_interval}. Given an interval, this outputs all the
#' segments that are equivalent partitioned between -pi/2, 0, and pi/2.
#'
#' @param interval vector of length 2
#'
#' @return 2-column matrix
.partition_interval <- function(interval){
  stopifnot(interval[1] < interval[2])
  a <- interval[1]; b <- interval[2]
  case <- abs(floor(a/(pi/2)) - floor(b/(pi/2)))
  if(case == 0) {
    mat <- matrix(c(a,b), ncol = 2)
  } else if (case == 1){
    if(a < 0) {
      mat <- rbind(c(a, 0), c(0, b))
    } else {
      mat <- rbind(c(-pi/2, b-pi), c(a, pi/2))
    }
  } else {
    stopifnot(b-a >= pi/2)
    if(a < 0){
      mat <- rbind(c(-pi/2, a), c(0, pi/2), c(b-pi, a))
    } else {
      mat <- rbind(c(-pi/2, 0), c(0, a), c(b-pi, pi/2))
    }
  }

  stopifnot(all(as.numeric(t(mat)) == sort(as.numeric(t(mat)))))
  mat
}

#' Intersect intervals
#'
#' \code{mat} is a 8-column matrix formed by concatenating
#' \code{.construct_interval} row-wise.
#'
#' Returns a matrix with 2 columns, where each row represents a
#' closed, connected interval of radians. The union of the rows represents
#' the intersection of all the intervals.
#'
#' @param mat matrix with 8 columns
#'
#' @return matrix
.intersect_intervals <- function(mat){

}

#' Convert euclidean points on circle into radians
#'
#' Only works if the circle is of the form: (x-a)^2 + (y-b)^2 = a^2+b^2.
#'
#' @param circle \code{circle} object
#' @param point point on circle, as a vector of length 2
#'
#' @return radian
.euclidean_to_radian <- function(circle, point, tol = 1e-6){
  stopifnot(length(point) == 2, length(circle$center) == 2)
  stopifnot(abs(circle$radius^2 - sum(circle$center^2)) < tol)
  stopifnot(abs(sum((point-circle$center)^2) - circle$radius^2) < tol)

  if(point[2] != 0){
    atan(point[1]/point[2])
  } else {
    atan(-circle$center[2]/circle$center[1])
  }
}
