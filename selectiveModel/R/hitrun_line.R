.sampler_hit_run_line <- function(){

}

#' Output the next point for hit-and-run sampler, Line
#'
#' For known sigma
#'
#' @param y data
#' @param segments matrix created by \code{.segments}
#' @param polyhedra \code{polyhedra} object
#' @param gaussian \code{gaussian} object
#'
#' @return vector
.hit_run_next_point_line <- function(y, segments, polyhedra, gaussian){
  stopifnot(length(y) == length(gaussian$mean))

  n <- length(y)
  v <- as.numeric(.sample_matrix_space(segments, 1, null = T))
  line <- .line(y, v)
  interval <- .intersect_polyhedron_line(polyhedra, line)

  rotation <- .rotation_matrix(v, c(1, rep(0, n-1)))
  gaussian <- .transform_gaussian(gaussian, y, rotation)

  univariate <- .conditional_gaussian(gaussian, rep(0, n-1))

  alpha <- .sampler_truncated_gaussian(univariate, interval[1], interval[2])
  y + alpha*v
}
