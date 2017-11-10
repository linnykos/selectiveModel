#' Hit and run sampler
#'
#' @param y data
#' @param seg_mean vector of segment means
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param num_samp number of desired samples from null distribution
#' @param cores umber of cores
#' @param burn_in positive integer, where we sample \code{num_samp*burn_in}
#' samples from the null distribution and return every \code{burn_in}th sample
#' @param verbose boolean
#'
#' @return matrix with \code{num_samp} columns and \code{length(y)} rows
.sampler_hit_run <- function(y, seg_mean, segments, polyhedra, num_samp,
                             cores = NA, burn_in = 3, verbose = F){

  #define a function that we'll be looping

  #initialize y

  #compute v and w

  #compute radius function

  #find range of radius function

  #sample next y

  #account for burn-in
}

#' Output the next point for hit-and-run sampler
#'
#' @param seg_mean vector of segment means
#' @param segments matrix created by \code{.segments}
#' @param initial initial data vector
#'
#' @return vector
.hit_run_next_point <- function(seg_mean, segments, initial){

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
    vec <- stat::rnorm(n)
    .projection_matrix(vec, mat)
  })

  if(num_vec > 1){
    for(i in 2:num_vec){
      vec_mat[,i] <- .projection_matrix(vec_mat[,i], t(vec_mat[,1:(i-1),drop = F]))
    }
  }

  sapply(1:ncol(vec_mat), function(x){vec_mat[,x]/.l2norm(vec_mat[,x])})
}
