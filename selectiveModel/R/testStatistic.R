#' Compute the next jump for a \code{bsFs} object
#'
#' @param obj \code{bsFs} object
#' @param y data
#' @param ... additional parameters
#'
#' @return \code{bsFs} object with one more jump
#' @export
next_jump.bsFs <- function(obj, y, ...){
  num_steps <- obj$numSteps
  binSegInf::binSeg_fixedSteps(y, num_steps+1)
}

#' Compute distance of next jump
#'
#' \code{jump_vec} is a vector of indices (not necessarily in any particular
#' order) less than or equal to \code{n}. A value of \code{i} in \code{jump_vec}
#' denotes that there is jump between index \code{i} and \code{i+1}.
#' \code{jump_vec} does not need to include \code{n}.
#'
#' @param y data
#' @param jump_vec_before vector of indices
#' @param jump_vec_after vector of indices with exactly one additional element
#' compared to \code{jump_vec_before}
#'
#' @return positive numeric
.jump_contrast <- function(y, jump_vec_before, jump_vec_after){
  new_jump <- which(!jump_vec_after %in% jump_vec_before)
  stopifnot(length(idx) == 1)

  n <- length(y)
  jump_sorted <- c(0, sort(jump_vec_after), n)
  idx <- which(jump_sorted == new_jump)
  stopifnot(idx != n)
  interval_left <- (jump_sorted[idx-1]+1):jump_sorted[idx]
  interval_right <- (jump_sorted[idx]+1):jump_sorted[idx+1]

  mean(y[interval_right]) - mean(y[interval_left])
}

#' Test statistic for testing empirical size of next jump
#'
#' @param y data
#' @param fit changepoint fitted model
#' @param ... nothing
#'
#' @return numeric
.next_jump_statistic <- function(y, fit, ...){
  next_fit <- next_jump(fit, y)
  .jump_contrast(y, binSegInf::jumps(fit), binSegInf::jumps(next_fit))
}
