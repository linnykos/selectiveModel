.gaussian <- function(mean, covariance){
  structure(list(mean = mean, covariance = covariance), class = "gaussian")
}

#' Compute the density of a Gaussian transformed
#'
#' Shifts the mean of a Gaussian and rotations the density
#'
#' @param gaussian a \code{gaussian} object
#' @param shift vector
#' @param rotation square matrix
#'
#' @return a \code{gaussian} object
#' @source \url{https://en.wikipedia.org/wiki/Covariance_matrix#Basic_properties}
.transform_gaussian <- function(gaussian, shift, rotation){
  stopifnot(length(gaussian$mean) == length(shift))
  stopifnot(all(dim(gaussian$covariance) == dim(gaussian$rotation)))

  gaussian$mean <- rotation%*%(gaussian$mean - shift)
  gaussian$covariance <- rotation%*%gaussian$covariance%*%t(rotation)

  gaussian
}

#' Compute the conditional density of a Gaussian distribution
#'
#' If \code{val} is length \code{k} long, then this function returns the
#' \code{gaussian} object corresponding to density conditioned on the last
#' \code{k} values equal to \code{val}. For example, if \code{gaussian} represents
#' a 5-dimensional distribution, and \code{val = c(5,10)}, then this function
#' returns a \code{gaussian} object corresponding to the Gaussian distribution
#' where the 4th coordinate is set to 5 and the 5th coordinate is set to 10.
#'
#' @param gaussian a \code{gaussian} object
#' @param val the values of the last coordinates
#'
#' @return a \code{gaussian} object
#' @source \url{https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions}
.conditional_gaussian <- function(gaussian, val){
  d <-  length(gaussian$mean)
  len <- d - length(val)
  inv <- solve(gaussian$covariance[(len+1):d,(len+1):d])

  mean <- gaussian$mean[1:len] + gaussian$covariance[1:len, (len+1):d]%*%inv%*%(val-gaussian$mean[(len+1):d])
  cov <- gaussian$covariance[1:len, 1:len] - gaussian$covariance[1:len, (len+1):d] %*% inv %*% t(gaussian$covariance[1:len, (len+1):d])

  .gaussian(mean = mean, covariance = cov)
}

.sampler_truncated_gaussian <- function(gaussian, lower, upper){

}

.segment_gaussian <- function(gaussian, segments, seg_means){

}
