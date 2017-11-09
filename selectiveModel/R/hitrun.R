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
#' @param segments 2-row matrix created by \code{.segments}
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
#'
#' @return vector
.projection <- function(vec1, vec2){

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

}

#' Sample vectors from the null space of a matrix
#'
#' This function assumes the row space of \code{mat} is the number of
#' rows in \code{mat}.
#'
#' @param mat matrix
#' @param num_vec positive integer
#'
#' @return a list of vectors
.sample_nullspace <- function(mat, num_vec = 2){

}
