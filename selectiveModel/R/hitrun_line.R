## code modified from https://github.com/selective-inference/R-software/blob/master/forLater/maxZ/funs.constraints.R

.sampler_hit_run_line <- function(start_y, gaussian, segments, polyhedra, num_samp = 100,
                                  burn_in = 500, lapse = 2, verbose = F){

  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  setting <- .remove_nullspace(gaussian, polyhedra, segments_full, mean_val)

}

.remove_nullspace <- function(gaussian, polyhedra, segments_full, mean_val){
  new_gaussian <- .remove_nullspace_gaussian(gaussian, segments_full, mean_val)
  new_polyhedra <- .remove_nullspace_polyhedra(polyhedra, segments_full, mean_val)

  forward_translation <- function(y){
    as.numeric(segments_full %*% y)[1:length(mean_val)]
  }

  backward_translation <- function(z){
    as.numeric(solve(segments_full) %*% c(z, mean_val))
  }

  list(gaussian = new_gaussian, polyhedra = new_polyhedra,
       forward_translation = forward_translation,
       backward_translation = backward_translation)
}

.remove_nullspace_gaussian <- function(gaussian, segments_full, mean_val){
  new_gaussian <- .gaussian(mean = segments_full %*% gaussian$mean,
                            covariance = segments_full %*% gaussian$covariance %*% t(segments_full))
  conditional_gaussian <- .conditional_gaussian(new_gaussian, val = mean_val)
}

.remove_nullspace_polyhedra <- function(polyhedra, segments_full, mean_val){
  stopifnot(nrow(segments_full) == ncol(segments_full))
  n <- ncol(segments_full)
  k <- length(mean_val)

  gamma <- polyhedra$gamma %*% solve(segments_full)
  u <- polyhedra$u - gamma[,(n-k+1):n]*mean_val
  gamma <- gamma[,1:(n-k)]

  binSegInf::polyhedra(gamma, u)
}

.whiten <- function(gaussian, polyhedra){

  n <- length(gaussian$mean)
  factor_cov <- .factor_covariance(gaussian$covariance)
  sqrt_cov <- factor_cov$sqrt_cov
  sqrt_inv <- factor_cov$sqrt_inv

  new_polyhedra <- .whiten_polyhedra(polyhedra, sqrt_cov, gaussian$mean)

  forward_translation <- function(y){
    as.numeric(sqrt_inv %*% (y - gaussian$mean))
  }

  backward_translation <- function(z){
    as.numeric(sqrt_cov %*% z) + gaussian$mean
  }

  list(gaussian = .gaussian(rep(0, n), diag(n)),
       polyhedra = new_polyhedra,
       forward_translation = forward_translation,
       backward_translation = backward_translation)
}

.whiten_polyhedra <- function(polyhedra, sqrt_cov, mean_vec){
  gamma <- polyhedra$gamma %*% sqrt_cov
  u <- polyhedra$u - polyhedra$gamma %*% mean_vec

  scaling <- sqrt(apply(gamma^2, 1, sum))
  gamma <- gamma/scaling
  u <- u/scaling

  binSegInf::polyhedra(gamma, u)
}

.factor_covariance = function(mat) {
  k <- Matrix::rankMatrix(mat)
  svd_X <- svd(S, nu=k, nv=k)
  sqrt_cov <- t(sqrt(svd_X$d[1:k]) * t(svd_X$u[,1:k]))
  sqrt_inv <- t((1. / sqrt(svd_X$d[1:k])) * t(svd_X$u[,1:k]))

  list(sqrt_cov=sqrt_cov, sqrt_inv=sqrt_inv)
}
