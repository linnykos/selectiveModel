#' Hit and run sampler, Line
#'
#' For known sigma
#'
#' @param start_y initial y to draw from
#' @param gaussian a \code{gaussian} object where the mean represents the data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param num_samp number of desired samples from null distribution
#' @param burn_in positive integer of the first few samples to throw out per core
#' @param lapse positive integer, where we sample \code{num_samp*burn_in}
#' samples from the null distribution and return every \code{burn_in}th sample
#' @param verbose boolean
#'
#' @return matrix with \code{num_samp} columns and \code{length(gaussian$mean)} rows
.sampler_hit_run_line <- function(start_y, gaussian, segments, polyhedra, num_samp = 100,
                                  burn_in = 500, lapse = 2, verbose = F){

  stopifnot(all(polyhedra$gamma %*% start_y >= polyhedra$u))
  num_col <- burn_in + num_samp*lapse

  prev_y <- start_y
  y_mat <- matrix(NA, nrow = length(start_y), ncol = num_samp)
  seq_idx <- burn_in + (1:num_samp)*lapse

  for(i in 1:num_col){
    next_y <- .hit_run_next_point_line(prev_y, segments, polyhedra, gaussian)
    if(i %in% seq_idx){
      y_mat[,which(seq_idx == i)] <- next_y
    }

    prev_y <- next_y
  }

  stopifnot(all(.try_polyhedra(y_mat, polyhedra)))
  y_mat
}

#' Output the next point for hit-and-run sampler, Line
#'
#' For known sigma
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param gaussian \code{gaussian} object
#'
#' @return vector
.hit_run_next_point_line <- function(y, segments, polyhedra, gaussian){
  stopifnot(length(y) == length(gaussian$mean))

  n <- length(y)
  v <- as.numeric(.sample_matrix_space(segments, 1, null = T))
  stopifnot(abs(.l2norm(v)-1) < 1e-6)
  line <- .line(y, v)
  interval <- .intersect_polyhedron_line(polyhedra, line)

  rotation <- .rotation_matrix(v, c(1, rep(0, n-1)))
  gaussian <- .transform_gaussian(gaussian, y, rotation)

  univariate <- .conditional_gaussian(gaussian, rep(0, n-1))

  alpha <- .sampler_truncated_gaussian(univariate, interval[1], interval[2])
  y + alpha*line$direction
}
