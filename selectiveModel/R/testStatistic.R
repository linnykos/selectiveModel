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
  binseginf::binSeg_fixedSteps(y, num_steps+1)
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
  new_jump <- jump_vec_after[which(!jump_vec_after %in% jump_vec_before)]
  stopifnot(length(new_jump) == 1)

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
#' @export
next_jump_statistic <- function(y, fit, ...){
  next_fit <- next_jump(fit, y)
  .jump_contrast(y, binseginf::jumps(fit), binseginf::jumps(next_fit))
}

#' Compute the height of a specific jump
#'
#' \code{jump_number} refers to the \code{jump_number}th jump starting from
#' the left. This function outputs the difference between the mean of \code{y}
#' to the left and to the right of that jump (right minus left)
#'
#' @param y data vector
#' @param fit fitted changepoint object
#' @param jump positive integer
#' @param ... not used
#'
#' @return numeric
#' @export
segment_difference <- function(y, fit, jump = 1, ...){
  stopifnot(jump %% 1 == 0, jump > 0)

  jumps <- sort(binseginf::jumps(fit))
  stopifnot(jump <= length(jumps))
  jump_idx <- jumps[jump]
  stopifnot(jump_idx < length(y))

  jumps <- sort(unique(c(0, jumps, length(y))))
  idx <- which(jumps == jump_idx)
  stopifnot(idx != 1)

  mean(y[(jumps[idx]+1):jumps[idx+1]]) -
    mean(y[(jumps[idx-1]+1):jumps[idx]])
}
