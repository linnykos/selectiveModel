rm(list=ls())
library(selectiveModel)
library(simulation)
library(binSegInf)

rule_closure <- function(n){
  function(vec, ...){
    lev <- vec[1]

    mean_vec <- c(rep(0, floor(n/2)), rep(lev, floor(n/2)))

    while(TRUE){
      y  <- mean_vec + stats::rnorm(length(mean_vec))
      fit <- fit_method(dat)
      jumps <- binSegInf::jumps(fit)

      if(jumps[1] == floor(n/2) & jumps[2] > jumps[1]) break()
    }

    y
  }
}

criterion_closure <- function(fit_method,
                              test_func = selectiveModel::segment_difference, num_samp = 2000, cores = NA){
  function(dat, vec, y, ...){
    fit <- fit_method(dat)
    ignore_jump <- 2

    known_sigma <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                            num_samp = num_samp, ignore_jump = ignore_jump, cores = cores, verbose = F, sigma = 1, param = list(burn_in = 2000, lapse = 1))

    unknown_sigma <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                              num_samp = num_samp, ignore_jump = ignore_jump, cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))


    c(known_sigma$pval, unknown_sigma$pval, binSegInf::jumps(fit))
  }
}

##################

n <- 20
trials <- 500
paramMat <- as.matrix(1)
fit_method <- function(x){binSegInf::binSeg_fixedSteps(x, numSteps = 2)}

rule <- rule_closure(n)
criterion <- criterion_closure(fit_method)

res <- simulation::simulation_generator(rule, criterion, paramMat, trials = trials, cores = 10,
                                        as_list = F, filepath = "../simulation/tmp.RData")

save.image("../simulation/flat_region_pvalue.RData")
