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
