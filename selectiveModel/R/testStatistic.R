#' Compute the next jump for a \code{bsFs} object
#'
#' @param obj \code{bsFs} object
#' @param y data
#' @param ... additional parameters
#'
#' @return \code{bsFs} object with one more jump
#' @export
next_jump.bsFs <- function(obj, y, ...){

}

#' Variance of data with respect to changepoint model
#'
#' The changepoint model is the Gaussian model where the mean is the same
#' in segments dictated by \code{jump_vec}.
#'
#' \code{jump_vec} is a vector of indices (not necessarily in any particular
#' order) less than or equal to \code{n}. A value of \code{i} in \code{jump_vec}
#' denotes that there is jump between index \code{i} and \code{i+1}.
#' \code{jump_vec} does not need to include \code{n}.
#'
#' @param y data
#' @param jump_vec vector of indices
#'
#' @return positive numeric
.model_variance <- function(y, jump_vec){

}

#' Compute distance of next jump
#'
#' @param y data
#' @param jump_vec_before vector of indices
#' @param jump_vec_after vector of indices with exactly one additional element
#' compared to \code{jump_vec_before}
#'
#' @return positive numeric
.jump_contrast <- function(y, jump_vec_before, jump_vec_after){

}
