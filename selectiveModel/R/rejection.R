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

}
