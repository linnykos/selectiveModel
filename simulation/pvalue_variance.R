rm(list=ls())
library(selectiveModel)
library(simulation)
library(binSegInf)

rule_closure <- function(n){
  function(vec, ...){
    y_list[[vec["Level"]]]
  }
}

criterion_closure <- function(fit_method,
                              test_func = selectiveModel::segment_difference,
                              verbose = T){
  function(dat, vec, ...){
    fit <- binSegInf::binSeg_fixedSteps(dat, 1)
    poly <- binSegInf::polyhedra(fit)
    contrast <- binSegInf::contrast_vector(fit, 1)

    selected <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                                          num_samp = vec["Samples"], ignore_jump = 1, cores = NA,
                                          verbose = F, sigma = 1,
                                          param = list(burn_in = vec["Burn.in"], lapse = vec["Lapse"]))

    selected$pval
  }
}


###################

n <- 10
trials <- 100
paramMat <- as.matrix(expand.grid(1:4, c(100, 1000, 7500), c(1, 10, 100), c(100, 500, 1000)))
colnames(paramMat) <- c("Level", "Burn.in", "Lapse", "Samples")
fit_method <- function(x){binSegInf::binSeg_fixedSteps(x, numSteps = 1)}

# create the data
y_list <- lapply(seq(0, 1.5, length.out = 4), function(lev){
  set.seed(10)
  mean_vec <- c(rep(0, floor(n/2)), rep(lev, floor(n/2)))

  while(TRUE){
    y <- mean_vec + stats::rnorm(length(mean_vec))
    fit <- binSegInf::binSeg_fixedSteps(y, 1)
    jumps <- binSegInf::jumps(fit)
    if(all(abs(jumps - floor(n/2)) <= 1)) break()
  }

  y
})


rule <- rule_closure(n)
criterion <- criterion_closure(fit_method)

res <- simulation::simulation_generator(rule, criterion, paramMat, trials = trials, cores = 10,
                                        as_list = F, filepath = "../simulation/tmp.RData")

save.image("../simulation/pvalue_variance.RData")
