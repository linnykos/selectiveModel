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
num_samp <- 4000
burn_in <- 4000
cores <- NA
doMC::registerDoMC(cores = 10)


trials <- 100
# run selective model test for known sigma, testing the first jump (aka: jump with the smallest index)
func <- function(x){
  print(x)
  set.seed(x*10)
  res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                  num_samp = num_samp, ignore_jump = 1, sigma = 1,
                                  cores = cores, verbose = F, param = list(burn_in = burn_in, lapse = 1))
  res$pval
}
res_vec_known_0 <- unlist(foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x)))
save.image("experiment_variance_pvalue.RData")

##

func <- function(x){
  print(x)
  set.seed(x*10)
  res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                  num_samp = num_samp, ignore_jump = 1,
                                  cores = cores, verbose = F, param = list(burn_in = burn_in, lapse = 1))
  res$pval
}
res_vec_unknown_0 <- unlist(foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x)))
save.image("experiment_variance_pvalue.RData")


#######################

set.seed(10)
lev <- 1
while(TRUE){
  dat <- rnorm(n, c(rep(0, n/5), rep(lev, n/5), rep(0, n/5), rep(-lev, n/5), rep(0, n/5)))
  res <- fit_method(dat)
  idx <- sort(binseginf::jumps(res))
  if(abs(idx[1] - n/5) <= 2) break()
}

func <- function(x){
  print(x)
  set.seed(x*10)
  res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                  num_samp = num_samp, ignore_jump = 1, sigma = 1,
                                  cores = cores, verbose = F, param = list(burn_in = burn_in, lapse = 1))
  res$pval
}
res_vec_known_signal <- unlist(foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x)))
save.image("experiment_variance_pvalue.RData")

##

func <- function(x){
  print(x)
  set.seed(x*10)
  res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                  num_samp = num_samp, ignore_jump = 1,
                                  cores = cores, verbose = F, param = list(burn_in = burn_in, lapse = 1))
  res$pval
}
res_vec_unknown_signal <- unlist(foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x)))
save.image("experiment_variance_pvalue.RData")
