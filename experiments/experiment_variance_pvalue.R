rm(list=ls())
library(binseginf)
library(selectiveModel)
set.seed(10)
n <- 200
k <- 4
lev <- 0
dat <- rnorm(n, c(rep(0, n/5), rep(lev, n/5), rep(0, n/5), rep(-lev, n/5), rep(0, n/5)))
fit_method <- function(x){binseginf::bsfs(x, numSteps = k)}
test_func <- selectiveModel::segment_difference
num_samp <- 2000
cores <- NA

# trials <- 100
# # run selective model test for known sigma, testing the first jump (aka: jump with the smallest index)
# res_vec_known_0 <- rep(NA, trials)
# for(x in 1:trials){
#   print(x)
#   set.seed(x*10)
#   res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
#                            num_samp = num_samp, ignore_jump = 1, sigma = 1,
#                            cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
#   res_vec_known_0[x] <- res$pval
# }
# save.image("experiment_variance_pvalue.RData")
#
# res_vec_unknown_0 <- rep(NA, trials)
# for(x in 1:trials){
#   print(x)
#   set.seed(x*10)
#   res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
#                                   num_samp = num_samp, ignore_jump = 1,
#                                   cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
#   res_vec_unknown_0[x] <- res$pval
# }
# save.image("experiment_variance_pvalue.RData")


#######################

set.seed(10)
n <- 200
k <- 1
lev <- 0
dat <- rnorm(n, c(rep(0, n/5), rep(lev, n/5), rep(0, n/5), rep(-lev, n/5), rep(0, n/5)))
fit_method <- function(x){binseginf::bsfs(x, numSteps = k)}
test_func <- selectiveModel::segment_difference
num_samp <- 2000
cores <- NA

trials <- 100
# run selective model test for known sigma, testing the first jump (aka: jump with the smallest index)
res_vec_known_4 <- rep(NA, trials)
for(x in 1:trials){
  print(x)
  set.seed(x*10)
  res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                  num_samp = num_samp, ignore_jump = 1, sigma = 1,
                                  cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
  res_vec_known_4[x] <- res$pval
}
save.image("experiment_variance_pvalue2.RData")

res_vec_unknown_4 <- rep(NA, trials)
for(x in 1:trials){
  print(x)
  set.seed(x*10)
  res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                  num_samp = num_samp, ignore_jump = 1,
                                  cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
  res_vec_unknown_4[x] <- res$pval
}

save.image("experiment_variance_pvalue2.RData")
