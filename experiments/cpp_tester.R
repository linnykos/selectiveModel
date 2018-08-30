rm(list=ls())
n = 6
set.seed(10)
k <- 2
dat <- rnorm(n)
fit_method <- function(x){binseginf::bsfs(x, numSteps = k)}
test_func <- selectiveModel::segment_difference
num_samp <- 2000
cores <- NA

i = 5
set.seed(10*i)
dat <- rnorm(n)

tmp <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                num_samp = num_samp, ignore_jump = 1,
                                cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))

#########################################

rm(list=ls())
i = 99
print(i)
set.seed(i)
y <- rnorm(6)
fit_method <- function(x){binseginf::bsfs(x, numSteps = 1)}
test_func <- selectiveModel::segment_difference
num_samp <- 2000
cores <- NA
res <- selected_model_inference(y, fit_method = fit_method, test_func = test_func,
                                num_samp = num_samp, ignore_jump = 1,
                                cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))

# sigma = NA
# ignore_jump = 1
# verbose = F
# param = list(burn_in = 2000, lapse = 1)
#
# n <- length(y)
# fit <- fit_method(y)
# polyhedra <- binseginf::polyhedra(fit)
# test_stat <- test_func(y, fit, jump = ignore_jump)
#
# #prepare sampler
# segments <- .segments(n, binseginf::jumps(fit), ignore_jump = ignore_jump)
# param <- .fill_in_arguments(param)
#
# burn_in = param$burn_in
# lapse = param$lapse
#
# num_col <- burn_in + num_samp*lapse
# y_mat <- matrix(NA, nrow = length(y), ncol = num_samp)
# seq_idx <- burn_in + (1:num_samp)*lapse
# prev_y <- y
#
# for(i in 1:513){
#   print(paste0(i , " in ", num_col, ": ", round(prev_y[1],3)))
#
#   next_y <- .hit_run_next_point_radial(prev_y, segments, polyhedra)
#   if(i %in% seq_idx){
#     y_mat[,which(seq_idx == i)] <- next_y
#   }
#
#   prev_y <- next_y
# }
#
# # # 514TH ONE
# # .hit_run_next_point_radial(prev_y, segments, polyhedra)
# y = prev_y
# tmp <- .sample_matrix_space(segments, 2, null = T)
# v <- tmp[,1]; w <- tmp[,2]
# x = 1
# .c_form_interval(polyhedra$gamma[x,], polyhedra$u[x], y, v, w)
#
# ##############
#
# load("../experiments/debugging_env.RData")
# c_form_interval(polyhedra$gamma[x,], polyhedra$u[x], y, v, w)
