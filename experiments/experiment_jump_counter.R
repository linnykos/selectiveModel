rm(list=ls())
library(simulation)
library(binseginf)

trials <- 12000
paramMat <- as.matrix(cbind(expand.grid(c(0,.25,.5,1,2,4), c(2:4)), 200))
colnames(paramMat) <- c("SnR", "k", "n")
middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}

rule <- function(vec){
  middle_mutation(lev = vec["SnR"], n = vec["n"]) + stats::rnorm(vec["n"])
}

criterion_closure <- function(method){
  function(dat, vec, y){
    fit <- method(dat, numSteps = vec["k"])
    binseginf::jumps(fit)
  }
}

criterion_bs <- criterion_closure(function(x){binseginf::bsfs(x, ...)})
criterion_fl <- criterion_closure(function(x){binseginf::fLasso_fixedSteps(x, ...)})

########################

bs_res <- simulation::simulation_generator(rule = rule, criterion = criterion_bs,
                                           paramMat = paramMat, trials = trials,
                                           cores = 15, as_list = F,
                                           filepath = "experiment_jump_counter_tmp.RData",
                                           verbose = T)
save.image("experiment_jump_counter_2.RData")

fl_res <- simulation::simulation_generator(rule = rule, criterion = criterion_fl,
                                           paramMat = paramMat, trials = trials,
                                           cores = 15, as_list = F,
                                           filepath = "experiment_jump_counter_tmp.RData",
                                           verbose = T)
save.image("experiment_jump_counter_2.RData")
