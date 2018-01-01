.gaussian <- function(mean, covariance){
  stopifnot(!is.matrix(mean))
  stopifnot(length(covariance) == 1 | is.matrix(covariance))

  if(length(covariance) == 1) covariance <- matrix(covariance)

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

  mean <- as.numeric(rotation%*%(gaussian$mean - shift))
  covariance <- rotation%*%gaussian$covariance%*%t(rotation)

  .gaussian(mean, covariance)
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
  inv <- solve(gaussian$covariance[(len+1):d,(len+1):d, drop = F])

  mean <- gaussian$mean[1:len] + gaussian$covariance[1:len, (len+1):d]%*%inv%*%(val-gaussian$mean[(len+1):d])
  cov <- gaussian$covariance[1:len, 1:len, drop = F] -
    gaussian$covariance[1:len, (len+1):d, drop = F] %*% inv %*% t(gaussian$covariance[1:len, (len+1):d, drop = F])

  .gaussian(mean = as.numeric(mean), covariance = cov)
}

#' Sample from a truncated gaussian
#'
#' This only works for one-dimensional Gaussian distributions
#'
#' @param gaussian a \code{gaussian} object
#' @param lower numeric
#' @param upper numeric
#'
#' @return numeric
#' @source \url{https://arxiv.org/pdf/0907.4010.pdf}
.sampler_truncated_gaussian <- function(gaussian, lower, upper){
  stopifnot(length(gaussian$mean) == 1, length(gaussian$covariance) == 1)

  lower <- (lower - gaussian$mean)/sqrt(gaussian$covariance)
  upper <- (upper - gaussian$mean)/sqrt(gaussian$covariance)

  # one-sided fix
  if(is.infinite(lower) | is.infinite(upper)){

    if(is.infinite(lower) & is.infinite(upper)){
      return(sqrt(gaussian$covariance)*stats::rnorm(1)+gaussian$mean)
    }

    if(is.infinite(lower)){
      lower <- upper - 10
    } else {
      upper <- lower + 10
    }
  }

  while(TRUE){
    z <- stats::runif(1, lower, upper)

    if(0 >= lower & 0 <= upper){
      thres <- exp(-z^2/2)
    } else if(upper < 0){
      thres <- exp((upper^2 - z^2)/2)
    } else {
      thres <- exp((lower^2 - z^2)/2)
    }

    u <- stats::runif(1)
    if(u <= thres) break()
  }

  as.numeric(sqrt(gaussian$covariance)*z+gaussian$mean)
}
