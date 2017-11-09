#' Compute next jump
#'
#' @param obj object
#' @param y data
#' @param ... additional parameters
#'
#' @return object with one more jump
#' @export
next_jumps <- function(obj, y, ...) {UseMethod("next_jumps")}
