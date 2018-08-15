#' Hit and run sampler, Radial
#'
#' For unknown sigma
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param num_samp number of desired samples from null distribution
#' @param cores umber of cores
#' @param burn_in positive integer of the first few samples to throw out per core
#' @param lapse positive integer, where we sample \code{num_samp*burn_in}
#' samples from the null distribution and return every \code{burn_in}th sample
#' @param verbose boolean
#'
#' @return matrix with \code{num_samp} columns and \code{length(y)} rows
.sampler_hit_run_radial <- function(y, segments, polyhedra, num_samp = 100,
                                    burn_in = 500, lapse = 2, verbose = F){
  num_col <- burn_in + num_samp*lapse
  y_mat <- matrix(NA, nrow = length(y), ncol = num_samp)
  seq_idx <- burn_in + (1:num_samp)*lapse
  polyhedra <- .remove_rows(polyhedra, .l2norm(y))

  prev_y <- y

  for(i in 1:num_col){
    next_y <- .hit_run_next_point_radial(prev_y, segments, polyhedra)
    if(i %in% seq_idx){
      y_mat[,which(seq_idx == i)] <- next_y
    }

    prev_y <- next_y
  }

  stopifnot(all(.try_polyhedra(y_mat, polyhedra)))
  y_mat
}

#' Output the next point for hit-and-run sampler, Radial
#'
#' For unknown sigma
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#'
#' @return vector
.hit_run_next_point_radial <- function(y, segments, polyhedra){
  stopifnot(ncol(segments) == length(y))

  tmp <- .sample_matrix_space(segments, 2, null = T)
  v <- tmp[,1]; w <- tmp[,2]
  interval <- .range_theta_polyhedra(y, v, w, polyhedra)

  if(nrow(interval) == 1) {
    theta <- stats::runif(1, interval[1], interval[2])
  } else {
    len <- interval[,2] - interval[,1]
    row_idx <- sample(1:nrow(interval), 1, prob = len)
    theta <- stats::runif(1, interval[row_idx,1], interval[row_idx,2])
  }

  y_new <- .radians_to_data(theta, y, v, w)

  stopifnot(.try_polyhedra(y_new, polyhedra))

  y_new
}

.remove_rows <- function(polyhedra, radius, tol = 1e-6){
  bool <- sapply(1:length(polyhedra$u), function(x){
    plane <- .plane(polyhedra$gamma[x,], polyhedra$u[x])
    n <- ncol(polyhedra$gamma)
    dis <- .distance_point_to_plane(rep(0, n), plane)
    dis <= radius + tol
  })

  idx <- which(bool)
  binSegInf::polyhedra(polyhedra$gamma[idx,,drop = F], polyhedra$u[idx])
}
