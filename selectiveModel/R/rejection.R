#' Rejection sampler
#'
#' Function returns with a warning if \code{time_limit} is reached.
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param num_samp number of desired samples from null distribution
#' @param cores umber of cores
#' @param time_limit number of seconds as the time limit for function
#' @param verbose boolean
#'
#' @return matrix with at most \code{num_samp} columns and \code{length(y)} rows
.sampler_rejection <- function(y, segments, polyhedra, num_samp = 100,
                               cores = 1, time_limit = 3600, verbose = F){
  seg_means <- .segment_means(y, segments)
  plane <- .plane(segments, seg_means)
  center <- .closest_point_to_origin(plane)
  major_radius <- .l2norm(y)
  minor_radius <- sqrt(major_radius^2 - .l2norm(center)^2)
  nullspace <- .sample_matrix_space(segments)

  if(!is.na(cores)) {
    doMC::registerDoMC(cores = cores)
    num_col <- ceiling(num_samp/cores)
  } else {
    num_col <- num_samp
  }
  n <- length(y)
  start_time <- proc.time()["elapsed"]

  func <- function(i){
    set.seed(i)
    j <- 1
    mat <- matrix(0, ncol = num_col, nrow = n)

    while(TRUE){
      if(verbose & i == 1 & j %% floor(num_col/10) == 0) cat('*')
      if(proc.time()["elapsed"] - start_time > time_limit) break()
      if(j > num_col) break()

      x <- .sample_sphere(ncol(nullspace))
      y_new <- .change_basis(x, center, nullspace, minor_radius)

      stopifnot(sum(abs(seg_means - .segment_means(y_new, segments))) < 1e-6)
      stopifnot(abs(.l2norm(y_new) - major_radius) < 1e-6)

      obj <- binSegInf::binSeg_fixedSteps(y_new, 1)
      #print(binSegInf::jump_sign(obj))

      if(all(polyhedra$gamma %*% y_new >= polyhedra$u)) {
        #print("score")
        mat[,j] <- y_new
        j <- j+1
      }
    }

    mat
  }

  if(!is.na(cores)) {
    i <- 0 #debugging reasons
    y_mat <- do.call(cbind, foreach::"%dopar%"(foreach::foreach(i = 1:cores),
                                               func(i)))
  } else {
    y_mat <- func(1)
  }

  #cleanup
  if(any(is.na(y_mat[1,]))){
    idx <- which(is.na(y_mat[1,]))
    y_mat <- y_mat[,-idx,drop = F]
    warning("Rejection sampler ran out of time")
  }

  y_mat
}
