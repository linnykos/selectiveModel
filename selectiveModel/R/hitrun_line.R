## code modified from https://github.com/selective-inference/R-software/blob/master/forLater/maxZ/funs.constraints.R

#' @useDynLib selectiveModel sample_truncnorm_white
.sampler_hit_run_line <- function(start_y, gaussian, segments, polyhedra, num_samp = 100,
                                  burn_in = 500){

  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%start_y)
  segments_full <- rbind(t(nullspace_mat), segments)

  setting_1 <- .remove_nullspace(gaussian, polyhedra, segments_full, mean_val)
  setting_2 <- .whiten(setting_1$gaussian, setting_1$polyhedra)
  new_polyhedra <- setting_2$polyhedra

  start_z <- setting_2$forward_translation(setting_1$forward_translation(start_y))
  n <- length(start_z)
  directions <- .generate_directions(n)
  alphas <- new_polyhedra$gamma %*% t(directions)
  slack <- new_polyhedra$gamma %*% start_z - new_polyhedra$u

  z_sample <- matrix(rep(0, n * num_samp), nrow = n, ncol = num_samp)

  result <- .C("sample_truncnorm_white",
              as.numeric(start_z),
              as.numeric(slack),
              as.numeric(t(directions)),
              as.numeric(alphas),
              output=z_sample,
              as.integer(nrow(new_polyhedra$gamma)),
              as.integer(nrow(directions)),
              as.integer(length(start_z)),
              as.integer(burn_in),
              as.integer(num_samp))
  z_sample <- result$output

  apply(z_sample, 2, function(x){
    setting_1$backward_translation(setting_2$backward_translation(x))
  })
}

.remove_nullspace <- function(gaussian, polyhedra, segments_full, mean_val){
  new_gaussian <- .remove_nullspace_gaussian(gaussian, segments_full, mean_val)
  new_polyhedra <- .remove_nullspace_polyhedra(polyhedra, segments_full, mean_val)

  n <- ncol(polyhedra$gamma)
  k <- length(mean_val)

  forward_translation <- function(y){
    stopifnot(sum(abs(as.numeric(segments_full %*% y)[(n-k+1):n] - mean_val)) < 1e-6)

    as.numeric(segments_full %*% y)[1:(n-k)]
  }

  backward_translation <- function(z){
    as.numeric(solve(segments_full) %*% c(z, mean_val))
  }

  list(gaussian = new_gaussian, polyhedra = new_polyhedra,
       forward_translation = forward_translation,
       backward_translation = backward_translation)
}

.remove_nullspace_gaussian <- function(gaussian, segments_full, mean_val){
  new_gaussian <- .gaussian(mean = as.numeric(segments_full %*% gaussian$mean),
                            covariance = segments_full %*% gaussian$covariance %*% t(segments_full))
  .conditional_gaussian(new_gaussian, val = mean_val)
}

.remove_nullspace_polyhedra <- function(polyhedra, segments_full, mean_val){
  stopifnot(nrow(segments_full) == ncol(segments_full))
  n <- ncol(segments_full)
  k <- length(mean_val)

  gamma <- polyhedra$gamma %*% solve(segments_full)
  u <- polyhedra$u - gamma[,(n-k+1):n,drop = F] %*% mean_val
  gamma <- gamma[,1:(n-k),drop = F]

  binSegInf::polyhedra(gamma, u)
}

.whiten <- function(gaussian, polyhedra){

  n <- length(gaussian$mean)
  sqrt_cov <- t(base::chol(gaussian$covariance))
  sqrt_inv <- solve(chol_mat)

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

  #remove redundant rows
  idx <- sort(unique(c(which(is.finite(u)), which(u >= 0))))
  gamma <- gamma[idx,,drop = F]
  u <- u[idx]

  binSegInf::polyhedra(gamma, u)
}

.generate_directions <- function(n){
  mat <- rbind(diag(rep(1, n)),
               matrix(rnorm(n^2), n, n))

  scaling <- apply(mat, 1, .l2norm)
  mat/scaling
}


