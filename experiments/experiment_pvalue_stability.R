rm(list=ls())
library(binseginf)
library(selectiveModel)
set.seed(10)
n <- 200
k <- 2
lev <- 1
fit_method <- function(x){binseginf::bsfs(x, numSteps = k)}
test_func <- selectiveModel::segment_difference
num_samp <- 2000
burn_in <- 2000
cores <- NA
doMC::registerDoMC(cores = 15)
middle_mutation <- function(lev, n=200){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+40)] <- lev
  mn
}

trials <- 2000

func <- function(x){
  print(x)
  set.seed(x*10)
  dat <- middle_mutation(lev) + stats::rnorm(n)
  fit <- fit_method(dat)

  jump_idx <- sort(binseginf::jumps(fit))

  res <- rep(NA, 2)
  if(abs(jump_idx[1] - n/2) <= 2){
    set.seed(x*10)
    tmp <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                       num_samp = num_samp, ignore_jump = 1,
                                       cores = cores, verbose = F, param = list(burn_in = burn_in, lapse = 1))
    res[1] <- tmp$pval
  }

  if(abs(jump_idx[2] - (n/2+40)) <= 2){
    set.seed(x*10)
    tmp <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                       num_samp = num_samp, ignore_jump = 2,
                                       cores = cores, verbose = F, param = list(burn_in = burn_in, lapse = 1))
    res[2] <- tmp$pval
  }

  res
}

res_mat <- unlist(foreach::"%dopar%"(foreach::foreach(x = 1:trials), func(x)))
save.image("experiment_pvalue_stability_unknown.RData")
