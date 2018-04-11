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
#' If \code{sigma} is \code{NA}, then the function conditions on \code{.l2norm(y)}.
#' Otherwise, \code{sigma} represents the standard deviation of the marginal distribution
#' of \code{y}, asserting that each coordinate is independent and Gaussian.
#'
#' @param y data of length n
#' @param fit_method function to fit changepoint model
#' @param test_func function to apply to each \code{y} and fitted model
#' @param num_samp number of desired samples from null distribution
#' @param sample_method either \code{hit_run} or \code{rejection}
#' @param sigma postive number or \code{NA}.
#' @param ignore_jump which jump to ignore
#' @param param additional parameters for \code{sample_method} passed in as a list
#' @param cores number of cores, only for rejection sampling
#' @param verbose boolean
#' @param ... optional inputs for \code{test_func}
#'
#' @return quantile
#' @export
selected_model_inference <- function(y, fit_method,
                                     test_func = next_jump_statistic,
                                     num_samp = 100,
                                     sample_method = "hitrun",
                                     sigma = NA,
                                     ignore_jump = NA,
                                     param = list(burn_in = default_burn_in(),
                                                  lapse = default_lapse()),
                                     cores = 1, verbose = T, ...){

  #fit the model, fit polyhedra, and compute test statistic on the observed model
  n <- length(y)
  fit <- fit_method(y)
  polyhedra <- binSegInf::polyhedra(fit)
  test_stat <- test_func(y, fit, jump = ignore_jump, ...)

  #prepare sampler
  segments <- .segments(n, binSegInf::jumps(fit), ignore_jump = ignore_jump)
  param <- .fill_in_arguments(param)

  #pass to sampler
  if(is.na(sigma)) {
    if(sample_method == "hitrun") {
      samples <- .sampler_hit_run_radial(y, segments, polyhedra, num_samp = num_samp,
                                         burn_in = param$burn_in, lapse = param$lapse,
                                         verbose = verbose)
    } else if(sample_method == "rejection") {
      samples <- .sampler_rejection_radial(y, segments, polyhedra, num_samp = num_samp,
                                           cores = cores, verbose = verbose)
    } else {
      stop("sample_method not appropriate")
    }
  } else {
    if(sample_method == "hitrun") {
      gaussian <- .gaussian(rep(0, n), sigma^2*diag(n))
      samples <- .sampler_hit_run_line(y, gaussian, segments, polyhedra, num_samp = num_samp,
                                       burn_in = param$burn_in)
    } else {
      stop("sample_method not appropriate")
    }
  }

  #for each sample, run the test function
  null_stat <- apply(samples, 2, function(x){
    test_func(x, fit, jump = ignore_jump, ...)
  })

  #compute the quantile
  pval <- length(which(abs(null_stat) >= abs(test_stat)))/length(null_stat)

  list(pval = pval, test_stat = test_stat, null_stat = null_stat)
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
#' @param ignore_jump which jump to ignore
#'
#' @return s-row matrix, where s is equal to \code{length(jump_vec+1)}
.segments <- function(n, jump_vec, ignore_jump = NA){
  if(all(is.na(jump_vec))) return(matrix(rep(1/n,n), nrow = 1))

  stopifnot(all(jump_vec >= 1), all(jump_vec <= n), all(jump_vec %% 1 == 0))
  jump_vec <- sort(jump_vec)
  if(!is.na(ignore_jump)){
    stopifnot(ignore_jump %% 1 == 0, ignore_jump > 0, ignore_jump <= length(jump_vec))
    jump_vec <- jump_vec[-ignore_jump]
  }
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

.fill_in_arguments <- function(param){
  burn_in <- ifelse(!is.null(param$burn_in), param$burn_in, default_burn_in())
  lapse <- ifelse(!is.null(param$lapse), param$lapse, default_lapse())

  list(burn_in = burn_in, lapse = lapse)
}

#' Default burn in
#'
#' @return numeric
#' @export
default_burn_in <- function(){2000}

#' Default lapse
#'
#' @return numeric
#' @export
default_lapse <- function(){2}
