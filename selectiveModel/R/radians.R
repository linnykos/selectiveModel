#' Find the range of theta's that are in the polyhedra
#'
#' Brute-force finds theta's such that for each theta,
#' \code{.radians_to_data(theta, y, v, w)} is inside the polyhedra.
#'
#' @param y vector
#' @param v unit vector
#' @param w unit vector orthogonal to \code{v}
#' @param polyhedra \code{polyhedra} object
#'
#' @return vector of possible theta's
.range_theta_polyhedra <- function(y, v, w, polyhedra){
  stopifnot(.try_polyhedra(y, polyhedra))
  stopifnot(abs(.l2norm(v) - 1) < 1e-4,  abs(.l2norm(w) - 1) < 1e-4, abs(v%*%w) < 1e-4)

  interval_list <- lapply(1:nrow(polyhedra$gamma), function(x){
    plane <- .plane(polyhedra$gamma[x,], polyhedra$u[x])
    plane <- .intersect_plane_basis(plane, y, v, w)
    if(any(is.na(plane))) return(matrix(c(-pi/2, pi/2), ncol = 2))
    center <- c(-y%*%v, -y%*%w)
    radius <- sqrt(sum(center^2))
    circle <- .circle(center, radius)
    dis <- .distance_point_to_plane(center, plane)

    if(dis >= radius){
      matrix(c(-pi/2, pi/2), ncol = 2)
    } else {
      mat <- .intersect_circle_line(plane, circle)
      vec <- apply(mat, 1, .euclidean_to_radian, circle = circle)
      init_theta <- .initial_theta(y, v, w)
      .interval(vec, init_theta)
    }
  })

  interval <- .intersect_intervals(interval_list)
  stopifnot(all(interval[,1] < interval[,2]))

  interval
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
  .partition_interval(interval)
}

#' Convert endpoints into an interval
#'
#' Interval is of length 2. \code{theta} dictates which side of the interval
#' is included.
#'
#' @param endpoints vector of length 2
#' @param theta initial theta
#'
#' @return vector of length 2
.basic_interval <- function(endpoints, theta){
  stopifnot(all(abs(endpoints) <= pi/2))

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
#' @param tol small positive number
#'
#' @return 2-column matrix
.partition_interval <- function(interval, tol = 1e-6){
  stopifnot(interval[1] < interval[2])
  a <- interval[1]; b <- interval[2]

  vec <- seq(-2*pi, 2*pi, by = pi/2)
  vec <- c(a, vec[intersect(which(vec >= a), which(vec <= b))], b)
  lis <- lapply(1:(length(vec)-1), function(x){c(vec[c(x,x+1)])})

  idx <- sapply(lis, function(x){mid <- mean(x); sign(mid)*ceiling(abs(mid)/(pi/2))})

  lis <- lapply(1:length(lis), function(x){
    if(abs(idx[x]) <= 1) return(lis[[x]])
    lis[[x]] - sign(idx[x])*pi
  })

  if(length(lis) == 1){
    mat <- matrix(lis[[1]], ncol = 2)
  } else {
    mat <- do.call(rbind, lis)
    mat <- mat[order(mat[,1]),]
  }

  idx <- which(mat[,2]-mat[,1] > tol)
  mat <- mat[idx,,drop = F]

  stopifnot(all(as.numeric(t(mat)) == sort(as.numeric(t(mat)))))
  stopifnot(all(abs(mat) <= pi/2))
  mat
}

#' Intersect intervals
#'
#' \code{lis} is a list of matrices created by \code{.interval}.
#'
#' Returns a matrix with 2 columns, where each row represents a
#' closed, connected interval of radians. The union of the rows represents
#' the intersection of all the intervals.
#'
#' This function errors if there is no intersection.
#'
#' @param lis list of matrices
#'
#' @return matrix
.intersect_intervals <- function(lis){
  Reduce(.intersect_two_intervals, lis)
}

#' Intersect two intervals
#'
#' Both \code{mat1} and \code{mat2} are outputs of \code{.interval}. This
#' finds all the radian intervals that lie in their intersection, and returns
#' it as another matrix that is similar in layout to those outputted by
#' \code{.interval}.
#'
#' This function errors if there is no intersection.
#'
#' @param mat1 matrix
#' @param mat2 matrix
#'
#' @return matrix
.intersect_two_intervals <- function(mat1, mat2){
  vec <- sort(unique(c(as.numeric(mat1), as.numeric(mat2))))
  vec <- sort(c(vec, .interpolate(vec)))

  bool_vec <- sapply(vec, function(x){
    all(.theta_in_interval(x, mat1), .theta_in_interval(x, mat2))
  })

  idx_mat <- .consecutive_true(bool_vec)

  mat <- matrix(vec[idx_mat], ncol = 2)
  mat <- mat[order(mat[,1]),,drop = F]

  stopifnot(all(as.numeric(t(mat)) == sort(as.numeric(t(mat)))))
  mat
}

#' Produce a vector with midpoints
#'
#' @param vec vector
#'
#' @return
.interpolate <- function(vec){
  n <- length(vec)
  sapply(2:n, function(x){
    mean(vec[(x-1):x])
  })
}

#' Is theta in an interval
#'
#' The \code{mat} is an output of \code{.interval}.
#'
#' @param theta radian
#' @param mat matrix
#'
#' @return boolean
.theta_in_interval <- function(theta, mat){
  stopifnot(abs(theta) <= pi/2)

  vec <- any(apply(mat, 1, function(x){
    x[1] <= theta & theta <= x[2]
  }))
}

#' Finding consecutive sets of TRUE's in a vector
#'
#' Returns a 2-column matrix where each row represents a separate pair (start
#' and end index) of consecutive TRUE's. This function ignores indices
#' for TRUE that are singleton.
#'
#' @param vec vector
#'
#' @return 2-column matrix
.consecutive_true <- function(vec){
  idx <- which(vec)
  if(length(idx) == 0) stop("No intersection")
  breakpoint <- which(sapply(2:length(idx), function(x){idx[x]-idx[x-1] != 1}))

  breakpoint <- c(0, breakpoint, length(idx))
  mat <- t(sapply(2:length(breakpoint), function(x){
    c(idx[breakpoint[x-1]+1], idx[breakpoint[x]])
  }))

  #remove singletons
  idx <- which(mat[,1] != mat[,2])
  if(length(idx) == 0) stop("No intersection")
  mat[idx,]
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
