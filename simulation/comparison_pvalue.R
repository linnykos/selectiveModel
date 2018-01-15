rm(list=ls())
library(selectiveModel)
library(simulation)
library(binSegInf)

rule_closure <- function(n){
  function(vec, ...){
    lev <- vec[1]

    mean_vec <- c(rep(0, floor(n/2)), rep(lev, floor(n/2)))

     while(TRUE){
       y <- mean_vec + stats::rnorm(length(mean_vec))
       fit <- binSegInf::binSeg_fixedSteps(y, 2)
       jumps <- binSegInf::jumps(fit)
       if(all(abs(jumps - floor(n/2)) <= 1)) break()
     }

    y
  }
}

criterion_closure <- function(fit_method,
                              test_func = selectiveModel::segment_difference, num_samp = 1000,
                              cores = NA, verbose = T){
  function(dat, vec, ...){
    fit <- binSegInf::binSeg_fixedSteps(dat, 1)
    poly <- binSegInf::polyhedra(fit)
    contrast <- binSegInf::contrast_vector(fit, vec[2])
    saturated_pval <- binSegInf:::poly.pval2(dat, poly, contrast, sigma = 1, bits = 1000)$pv

    # selected <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
    #                                       num_samp = num_samp, ignore_jump = vec[2], cores = cores,
    #                                       verbose = F, sigma = 1,
    #                                       param = list(burn_in = 7500, lapse = 10))

    # print(paste0("saturated: ", round(saturated_pval,3), "// selected: ", round(selected$pval,3)))
    # c(saturated_pval, selected$pval)
    # c(selected$pval, binSegInf::jumps(fit))
    saturated_pval
  }
}

##################

n <- 10
trials <- 1000
paramMat <- cbind(seq(0, 1.5, length.out = 4), 1)
fit_method <- function(x){binSegInf::binSeg_fixedSteps(x, numSteps = 1)}

rule <- rule_closure(n)
criterion <- criterion_closure(fit_method)

res <- simulation::simulation_generator(rule, criterion, paramMat, trials = trials, cores = 1,
                                 as_list = F, filepath = "../simulation/tmp.RData")

save.image("../simulation/comparison_pvalue_saturated.RData")
