#' Selected model inference
#'
#' Fits a changepoint model using \code{fit_method} to \code{y}. The returned
#' model needs to have a function defined via the generic function
#' \code{binSeg::polyhedra}. Inference is performed by sampling data from the
#' null distribution matching the estimated changepoint model
#' where each sample is drawn from a Gaussian and matches
#' in empirical means in each segment and in L2-norm (the sufficient statistics).
#'
#' Currently the test is a 2-sided test.
#'
#' The \code{test_func} should be a function that takes in \code{y} as the first
#' input and the output of \code{fit_method} as the second input, and \code{...}
#' as the third input.
#'
#' @param y data of length n
#' @param fit_method function to fit changepoint model
#' @param test_func function to apply to each \code{y} and fitted model
#' @param num_samp number of desired samples from null distribution
#' @param sample_method either \code{hit_run} or \code{rejection}
#' @param param additional parameters for \code{sample_method} passed in as a list
#' @param cores number of cores
#' @param verbose boolean
#' @param ... optional inputs for \code{test_func}
#'
#' @return quantile
#' @export
selected_model_inference <- function(y, fit_method,
                                     test_func = next_jump_statistic,
                                     num_samp = 100,
                                     sample_method = "hitrun",
                                     param = list(burn_in = 2, time_limit = 3600),
                                     cores = 1, verbose = T, ...){

  #fit the model, fit polyhedra, and compute test statistic on the observed model
  n <- length(y)
  fit <- fit_method(y)
  polyhedra <- binSegInf::polyhedra(fit)
  test_stat <- test_func(y, fit, ...)

  #prepare sampler
  segments <- .segments(n, binSegInf::jumps(fit))

  #pass to sampler
  if(sample_method == "hitrun") {
    stopifnot(c("burn_in") %in% names(param))
    samples <- .sampler_hit_run(y, segments, polyhedra, num_samp = num_samp,
                                cores = cores, burn_in = param$burn_in,
                                verbose = verbose)
  } else if(sample_method == "rejection") {
    stopifnot(c("time_limit") %in% names(param))
    samples <- .sampler_rejection(y, segments, polyhedra, num_samp = num_samp,
                                cores = cores, time_limit = param$time_limit,
                                verbose = verbose)
  } else {
    stop("sample_method not appropriate")
  }

  #for each sample, run the test function
  null_stat <- apply(samples, 2, function(x){
    test_func(x, fit, ...)
  })

  #compute the quantile
  length(which(abs(null_stat) >= abs(test_stat)))/length(null_stat)
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
#' \code{jump_vec} does not need to include \code{n}.
#'
#' @param n number of samples
#' @param jump_vec vector of indices
#'
#' @return s-row matrix, where s is equal to \code{length(jump_vec+1)}
.segments <- function(n, jump_vec){
  stopifnot(all(jump_vec >= 1), all(jump_vec <= n), all(jump_vec %% 1 == 0))
  jump_vec <- unique(sort(c(0,jump_vec, n)))

  len <- length(jump_vec)
  mat <- sapply(1:(len-1), function(x){
    vec <- rep(0, n)
    vec[(jump_vec[x]+1):jump_vec[(x+1)]] <- 1/(jump_vec[x+1]-jump_vec[x])
    vec
  })

  t(mat)
}

#' Computes the empirical mean of each segment
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#'
#' @return vector
.segment_means <- function(y, segments){
  stopifnot(!is.matrix(y), length(y) == ncol(segments))
  as.numeric(segments %*% y)
}

.l2norm <- function(vec){
  as.numeric(sqrt(sum(vec^2)))
}
