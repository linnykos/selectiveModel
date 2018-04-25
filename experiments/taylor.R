rm(list=ls())
set.seed(10)
#devtools::install_github("selective-inference/R-software", subdir = "selectiveInference-currentCRAN")
library("selectiveInference")
source("../experiments/taylor_helper.R")
source("../experiments/taylor_helper2.R")

n <- 6
trials <- 5000

# determine polyhedron
y <- rnorm(n)
obj <- binSegInf::binSeg_fixedSteps(y, 1)
poly <- binSegInf::polyhedra(obj)

segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj), ignore_jump = 1)
segments_full <- rbind(c(1, rep(0, n-1)), c(0, 1, rep(0, n-2)),
                       c(0,0,1,0,0,0), c(0,0,0,1,0,0), c(0,0,0,0,1,0),
                       segments)
mean_val <- as.numeric(segments%*%y)
cov_mat <- diag(n)

cov_mat_full <- segments_full %*% cov_mat %*% t(segments_full)
mean_cond <- cov_mat_full[1:5,6]%*%solve(cov_mat_full[6,6])%*%mean_val
cov_cond <- cov_mat_full[1:5,1:5] - cov_mat_full[1:5,6]%*%solve(cov_mat_full[6,6])%*%t(cov_mat_full[1:5,6])

new_poly_gamma <- poly$gamma %*% solve(segments_full)

# check
z <- segments_full %*% y
stopifnot(all(new_poly_gamma %*% z >= poly$u))

new_poly_u <- poly$u - new_poly_gamma[,6]*mean_val
new_poly_gamma <- new_poly_gamma[,1:5]

#check again
stopifnot(all(new_poly_gamma %*% z[1:5] >= new_poly_u))

#flip constraints to match specification
samp <- sample_from_constraints(linear_part = -new_poly_gamma,
                                offset = -new_poly_u,
                                mean_param = mean_cond,
                                covariance = cov_cond,
                                initial_point = z[1:5],
                                ndraw = 500, burnin = 100)

#convert back to y's
samp <- apply(samp, 1, function(x){
  solve(segments_full) %*% c(x, mean_val)
})

#ensure in constraint
bool_vec <- apply(samp, 2, function(x){
  all(poly$gamma %*% x >= poly$u)
})
stopifnot(all(bool_vec))

#now compute the distribution
segment_real <-  .segments(length(y), jump_vec = binSegInf::jumps(obj))
segment_dist <- apply(samp, 2, function(x){
  diff(segment_real %*% x)
})
