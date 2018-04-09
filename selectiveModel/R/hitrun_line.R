## code modified from https://github.com/selective-inference/R-software/blob/master/forLater/maxZ/funs.constraints.R

.sampler_hit_run_line <- function(start_y, gaussian, segments, polyhedra, num_samp = 100,
                                  burn_in = 500, verbose = F){

  n <- length(start_y)
  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  setting_1 <- .remove_nullspace(gaussian, polyhedra, segments_full, mean_val)
  setting_2 <- .whiten(setting_1$gaussian, setting_1$polyhedra)
  new_polyhedra <- setting_2$polyhedra

  directions <- .generate_directions(n)
  start_z <- setting_2$forward_translation(setting_1$forward_translation(start_y))
  alphas <- new_polyhedra$gamma %*% directions
  slack <- new_polyhedra$gamma %*% start_z - new_polyhedra$u

  num_draw <- burn_in+num_samp
  z_sample <- matrix(rep(0, n * num_draw), nrow = n, ncol = num_draw)
  result <- .C("sample_truncnorm_white",
              as.numeric(start_z),
              as.numeric(slack),
              as.numeric(t(directions)),
              as.numeric(alphas),
              output=z_sample,
              as.integer(nconstraint),
              as.integer(ndirection),
              as.integer(nstate),
              as.integer(burnin),
              as.integer(ndraw),
              package="selectiveModel")
  z_sample <- result$output

  apply(z_sample, 2, function(x){
    setting_1$backward_translation(setting_2$backward_translation(x))
  })
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

  #remove redundant rows
  idx <- sort(unique(c(which(is.finite(u)), which(u >= 0))))
  gamma <- gamma[idx,,drop = F]
  u <- u[idx]

  binSegInf::polyhedra(gamma, u)
}

.factor_covariance = function(mat){
  k <- Matrix::rankMatrix(mat)
  svd_X <- svd(S, nu=k, nv=k)
  sqrt_cov <- t(sqrt(svd_X$d[1:k]) * t(svd_X$u[,1:k]))
  sqrt_inv <- t((1. / sqrt(svd_X$d[1:k])) * t(svd_X$u[,1:k]))

  list(sqrt_cov=sqrt_cov, sqrt_inv=sqrt_inv)
}

.generate_directions <- function(n){
  mat <- rbind(diag(rep(1, n)),
               matrix(rnorm(n^2), n, n))

  scaling <- apply(mat, 1, .l2norm)
  mat / scaling
}


