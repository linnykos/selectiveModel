#' Hit and run sampler
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param num_samp number of desired samples from null distribution
#' @param cores umber of cores
#' @param burn_in positive integer, where we sample \code{num_samp*burn_in}
#' samples from the null distribution and return every \code{burn_in}th sample
#' @param verbose boolean
#' @param tol small positive number
#' @param max_iter maximum number of iterations
#'
#' @return matrix with \code{num_samp} columns and \code{length(y)} rows
.sampler_hit_run <- function(y, segments, polyhedra, num_samp,
                             cores = 1, burn_in = 3, tol = 1e-2, max_iter = 10){
  doMC::registerDoMC(cores = cores)
  n <- length(y)
  num_col <- ceiling(num_samp*burn_in/cores)

  func <- function(i){
    mat <- matrix(0, ncol = num_col, nrow = n)
    prev_y <- y

    for(j in 1:num_col){
      if(verbose & i == 1 & j %% floor(num_col/10) == 0) cat('*')
      mat[,j] <- .hit_run_next_point(prev_y, segments, polyhedra, tol, max_iter)
      prev_y <- mat[,j]
    }

    mat
  }

  i <- 0 #debugging reasons
  y_new_mat <- do.call(cbind, foreach::"%dopar%"(foreach::foreach(i = 1:cores),
                     func(i)))

  seq_vec <- seq(burn_in, ncol(y_new_mat), by = burn_in)
  stopifnot(length(seq_vec) < num_samp)
  if(length(seq_vec) > num_samp) seq_vec <- seq_vec[1:num_samp]

  y_new_mat[,seq_vec]
}

#' Output the next point for hit-and-run sampler
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param tol small positive number
#' @param max_iter maximum number of iterations
#'
#' @return vector
.hit_run_next_point <- function(y, segments, polyhedra, tol = 1e-2, max_iter = 10){
  stopifnot(ncol(segments) == length(y))

  tmp <- .sample_nullspace(segments, 2)
  v <- tmp[,1]; w <- tmp[,2]
  interval <- .range_theta_polyhedra(y, v, w, polyhedra, tol, max_iter)
  theta <- runif(1, interval[1], interval[2])
  .radians_to_data(theta, y, v, w)
}

#' Projection of vector onto another vector
#'
#' Returns the component of \code{vec1} that is orthogonal to \code{vec2}
#'
#' @param vec1 vector
#' @param vec2 vector
#' @param tol small positive number
#'
#' @return vector
.projection <- function(vec1, vec2, tol = 1e-6){
  stopifnot(length(vec1) == length(vec2))

  d <- length(vec1)
  vec2 <- vec2/.l2norm(vec2)
  as.numeric((diag(d) - vec2%*%t(vec2))%*%vec1)
}

#' Projection of vector onto rows of a matrix
#'
#' Returns the component of \code{vec} that is orthogonal to all the rows
#' of \code{mat}
#'
#' @param vec vector
#' @param mat matrix
#'
#' @return vector
.projection_matrix <- function(vec, mat){
  stopifnot(length(vec) == ncol(mat), ncol(mat) >= nrow(mat))
  n <- length(vec)
  proj_mat <- diag(n) - t(mat) %*% solve(mat %*% t(mat)) %*% mat

  as.numeric(proj_mat %*% vec)
}

#' Sample unit vectors from the null space of a matrix
#'
#' This function assumes the row space of \code{mat} is the number of
#' rows in \code{mat}.
#'
#' @param mat matrix
#' @param num_vec positive integer
#'
#' @return a matrix of vectors, with number of rows equal to \code{ncol(mat)}
#' and number of columns equal to \code{num_vec}
.sample_nullspace <- function(mat, num_vec = 2){
  stopifnot(num_vec > 0, num_vec %% 1 == 0)
  stopifnot(nrow(mat) + num_vec <= ncol(mat))

  n <- ncol(mat)
  vec_mat <- sapply(1:num_vec, function(x){
    vec <- stats::rnorm(n)
    .projection_matrix(vec, mat)
  })

  if(num_vec > 1){
    for(i in 2:num_vec){
      vec_mat[,i] <- .projection_matrix(vec_mat[,i], t(vec_mat[,1:(i-1),drop = F]))
    }
  }

  sapply(1:ncol(vec_mat), function(x){vec_mat[,x]/.l2norm(vec_mat[,x])})
}
