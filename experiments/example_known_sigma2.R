rm(list=ls())
library(binseginf)
library(selectiveModel)
set.seed(10)
n <- 5
k <- 1
dat <- rnorm(n, c(rep(0, n/5), rep(1, n/5), rep(0, n/5), rep(-1, n/5), rep(0, n/5)))
fit_method <- function(x){binseginf::bsfs(x, numSteps = k)}
test_func <- selectiveModel::segment_difference
num_samp <- 10
cores <- NA

# set.seed(10)
# # run selective model test for known sigma, testing the first jump (aka: jump with the smallest index)
# system.time(res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
#                          num_samp = num_samp, ignore_jump = 1, sigma = 1,
#                          cores = cores, verbose = F, param = list(burn_in = 200, lapse = 1)))
#
# res$pval

set.seed(10)
system.time(res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                            num_samp = num_samp, ignore_jump = 1,
                                            cores = cores, verbose = F, param = list(burn_in = 10, lapse = 1)))

############################

rm(list=ls())
library(binseginf)
library(selectiveModel)
set.seed(10)
n <- 5
k <- 1
dat <- rnorm(n, c(rep(0, n/5), rep(1, n/5), rep(0, n/5), rep(-1, n/5), rep(0, n/5)))
fit_method <- function(x){binseginf::bsfs(x, numSteps = k)}
test_func <- selectiveModel::segment_difference
num_samp <- 10
cores <- NA

set.seed(10)
y <- dat
ignore_jump = 1
param = list(burn_in = 10, lapse = 1)
n <- length(y)
fit <- fit_method(y)
polyhedra <- binseginf::polyhedra(fit)
test_stat <- test_func(y, fit, jump = ignore_jump)

#prepare sampler
segments <- selectiveModel:::.segments(n, binseginf::jumps(fit), ignore_jump = ignore_jump)
param <- selectiveModel:::.fill_in_arguments(param)
lapse = param$lapse
burn_in = param$burn_in

num_col <- burn_in + num_samp*lapse
y_mat <- matrix(NA, nrow = length(y), ncol = num_samp)
seq_idx <- burn_in + (1:num_samp)*lapse
prev_y <- y
null_mat <- selectiveModel:::.compute_nullspace(segments)

for(i in 1:13){
  print(paste0("i: ", i))
  next_y <- selectiveModel:::.hit_run_next_point_radial(prev_y, null_mat, polyhedra)
  if(i %in% seq_idx){
    y_mat[,which(seq_idx == i)] <- next_y
  }

  prev_y <- next_y
}

y <- prev_y
i = 14
tmp <- selectiveModel:::.sample_two_vectors(null_mat)
v <- tmp[,1]; w <- tmp[,2]
interval <- selectiveModel:::.range_theta_polyhedra(y, v, w, polyhedra)
interval_list <- lapply(1:nrow(polyhedra$gamma), function(x){
  selectiveModel:::.c_form_interval(polyhedra$gamma[x,], polyhedra$u[x], y, v, w)
})
print(interval_list)

interval <- selectiveModel:::.c_intersect_intervals(interval_list)
print(interval)


