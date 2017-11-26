rm(list=ls())
library(selectiveModel)
library(simulation)

rule_closure <- function(n){
  function(vec, ...){
    lev <- vec[1]

    n5 <- ceiling(n/5)
    mean_vec <- c(rep(0, n5), rep(lev, n5), rep(0, n5), rep(-2*lev, n5), rep(0, n5))
    mean_vec + stats::rnorm(length(mean_vec))
  }
}

criterion_closure <- function(fit_method,
                              test_func = selectiveModel::segment_difference, num_samp = 250,
                              cores = 14, verbose = F){
  function(dat, vec, ...){
    res <- selectiveModel::selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                             num_samp = num_samp, ignore_jump = vec[2], cores = cores,
                             verbose = verbose)
    res$pval
  }
}

##################

n <- 60
trials <- 250
paramMat <- as.matrix(expand.grid(c(0, 0.5, 1, 1.5, 2), c(1:4)))
fit_method = function(x){binSegInf::binSeg_fixedSteps(x, numSteps = 4)}

rule <- rule_closure(n)
criterion <- criterion_closure(fit_method)

simulation::simulation_generator(rule, criterion, paramMat, trials = trials, cores = 1,
                                 as_list = F, filepath = "../simulation/tmp.RData")
