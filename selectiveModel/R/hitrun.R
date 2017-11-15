#' Hit and run sampler
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param num_samp number of desired samples from null distribution
#' @param cores umber of cores
#' @param burn_in positive integer, where we sample \code{num_samp*burn_in}
#' samples from the null distribution and return every \code{burn_in}th sample
#' @param seed numeric to adjust the seed of each core
#' @param verbose boolean
#'
#' @return matrix with \code{num_samp} columns and \code{length(y)} rows
.sampler_hit_run <- function(y, segments, polyhedra, num_samp = 100,
                             cores = 1, burn_in = 2,
                             seed = 1, verbose = F){
  if(!is.na(cores)) {
    doMC::registerDoMC(cores = cores)
    num_col <- ceiling(num_samp*burn_in/cores)
  } else {
    num_col <- num_samp*burn_in
  }

  n <- length(y)

  func <- function(i){
    mat <- matrix(0, ncol = num_col, nrow = n)
    prev_y <- y

    for(j in 1:num_col){
      set.seed((j^i)*seed)
      if(verbose & i == 1 & j %% floor(num_col/10) == 0) cat('*')
      mat[,j] <- .hit_run_next_point(prev_y, segments, polyhedra)
      prev_y <- mat[,j]
    }

    mat
  }

  i <- 0 #debugging reasons
  if(!is.na(cores)) {
    y_mat <- do.call(cbind, foreach::"%dopar%"(foreach::foreach(i = 1:cores),
                     func(i)))
  } else {
    y_mat <- func(1)
  }

  seq_vec <- seq(burn_in, ncol(y_mat), by = burn_in)
  stopifnot(length(seq_vec) >= num_samp)
  if(length(seq_vec) > num_samp) seq_vec <- seq_vec[1:num_samp]

  y_mat <- y_mat[,seq_vec]
  stopifnot(all(.try_polyhedra(y_mat, polyhedra)))
  y_mat
}

#' Output the next point for hit-and-run sampler
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#'
#' @return vector
.hit_run_next_point <- function(y, segments, polyhedra){
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
