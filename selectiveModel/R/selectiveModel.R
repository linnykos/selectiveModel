selective_model_inference <- function(y, fit_method, test_func, num_samp, sample_method,
                                      param = list(), cores = NA, verbose = T){

  #fit the model

  #create the polyhedra

  #run the test function on the observed model

  #pass to sampler

  #for each sample, run the test function

  #compute the quantile
}

#' Compute the matrix of segments
#'
#' Return a matrix \code{mat} such that
#' for a data vector \code{y}, we have
#' \code{mat %*% y} would output the empirical means of each segment of \code{y}
#' according to \code{jump_vec}.
#'
#' \code{jump_vec} is a vector of indices (not necessarily in any particular
#' order) less than or equal to \code{n}. A value of \code{i} in \code{jump_vec}
#' denotes that there is jump between index \code{i} and \code{i+1}.
#'
#' @param n number of samples
#' @param jump_vec location of jumps.
#'
#' @return 2-column matrix
.segments <- function(n, jump_vec){

}

.segment_means <- function(y, segments){

}

.l2norm <- function(vec){
  sqrt(sum(vec^2))
}
