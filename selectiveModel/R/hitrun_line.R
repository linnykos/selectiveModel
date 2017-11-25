#' Hit and run sampler, Line
#'
#' For known sigma
#'
#' @param gaussian a \code{gaussian} object where the mean represents the data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param num_samp number of desired samples from null distribution
#' @param cores umber of cores
#' @param burn_in positive integer of the first few samples to throw out per core
#' @param lapse positive integer, where we sample \code{num_samp*burn_in}
#' samples from the null distribution and return every \code{burn_in}th sample
#' @param verbose boolean
#'
#' @return matrix with \code{num_samp} columns and \code{length(gaussian$mean)} rows
.sampler_hit_run_line <- function(gaussian, segments, polyhedra, num_samp = 100,
                                  cores = 1, burn_in = 500, lapse = 2, verbose = F){
  if(!is.na(cores)) {
    doMC::registerDoMC(cores = cores)
    num_col <- ceiling(burn_in + num_samp*lapse/cores)
  } else {
    num_col <- burn_in + num_samp*lapse
  }

  y <- gaussian$mean
  n <- length(y)

  func <- function(i){
    set.seed(i)
    mat <- matrix(0, ncol = num_col, nrow = n)
    prev_y <- y

    for(j in 1:num_col){
      if(verbose & i == 1 & j %% floor(num_col/10) == 0) cat('*')
      mat[,j] <- .hit_run_next_point_line(prev_y, segments, polyhedra, gaussian)
      prev_y <- mat[,j]
    }

    mat[,burn_in:ncol(mat)]
  }

  i <- 0 #debugging reasons
  if(!is.na(cores)) {
    y_mat <- do.call(cbind, foreach::"%dopar%"(foreach::foreach(i = 1:cores),
                                               func(i)))
  } else {
    y_mat <- func(1)
  }

  seq_vec <- seq(1, ncol(y_mat), by = lapse)
  stopifnot(length(seq_vec) >= num_samp)
  if(length(seq_vec) > num_samp) seq_vec <- seq_vec[1:num_samp]

  y_mat <- y_mat[,seq_vec]
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
  c <- y + interval[1]*v
  d <- y + interval[2]*v

  gaussian <- .transform_gaussian(gaussian, c, rotation)

  # checking
  direction <- d-c
  vec <- rotation%*%direction
  vec <- vec/.l2norm(vec)
  stopifnot(sum(abs(vec - c(1, rep(0, n-1)))) < 1e-6)

  univariate <- .conditional_gaussian(gaussian, rep(0, n-1))

  alpha <- .sampler_truncated_gaussian(univariate, 0, interval[2]-interval[1])
  y + (alpha+interval[1])*v
}
