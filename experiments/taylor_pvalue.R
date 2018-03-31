library("selectiveInference")
source("../experiments/taylor_helper2.R")

taylor_pvalue_hitrun_sampler <- function(y, gaussian, poly, segments, num_samples = 1000, burn_in = 2000){
  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  n <- ncol(segments)
  k <- nrow(segments)
  z <- segments_full %*% y

  cov_mat_full <- segments_full %*% gaussian$covariance %*% t(segments_full)
  mean_cond <- cov_mat_full[1:(n-k),(n-k+1):n]%*%solve(cov_mat_full[(n-k+1):n,(n-k+1):n])%*%mean_val
  cov_cond <- cov_mat_full[1:(n-k),1:(n-k)] - cov_mat_full[1:(n-k),(n-k+1):n]%*%solve(cov_mat_full[(n-k+1):n,(n-k+1):n])%*%t(cov_mat_full[1:(n-k),(n-k+1):n])

  new_poly_gamma <- poly$gamma %*% solve(segments_full)

  new_poly_u <- poly$u - new_poly_gamma[,(n-k+1):n]*mean_val
  new_poly_gamma <- new_poly_gamma[,1:(n-k)]

  samp <- sample_from_constraints(linear_part = -new_poly_gamma,
                                  offset = -new_poly_u,
                                  mean_param = mean_cond,
                                  covariance = cov_cond,
                                  initial_point = z[1:(n-k)],
                                  ndraw = num_samples, burnin = burn_in)

  apply(samp, 1, function(x){
    solve(segments_full) %*% c(x, mean_val)
  })
}
