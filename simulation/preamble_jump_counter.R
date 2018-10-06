rm(list=ls())
library(simulation)
library(binseginf)

trials <- 10000
paramMat <- as.matrix(expand.grid(c(0,.25,.5,1,2,4), c(2:4), 200))
colnames(paramMat) <- c("SnR", "k", "n")
middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}

rule <- function(vec){
  middle_mutation(lev = vec["SnR"], n = vec["n"]) + stats::rnorm(vec["n"])
}

criterion_bs <- function(dat, vec, y){
  fit <- binseginf::bsfs(dat, numSteps = vec["k"])
  binseginf::jumps(fit)
}

criterion_fl <- function(dat, vec, y){
  fit <- binseginf::fLasso_fixedSteps(dat, numSteps = vec["k"])
  binseginf::jumps(fit)
}

########################

bs_res <- simulation::simulation_generator(rule = rule, criterion = criterion_bs,
                                           paramMat = paramMat, trials = trials,
                                           cores = 15, as_list = F,
                                           filepath = "preamble_jump_counter_tmp.RData",
                                           verbose = T)
save.image("preamble_jump_counter.RData")

fl_res <- simulation::simulation_generator(rule = rule, criterion = criterion_fl,
                                           paramMat = paramMat, trials = trials,
                                           cores = 15, as_list = F,
                                           filepath = "experiment_jump_counter_tmp.RData",
                                           verbose = T)
save.image("preamble_jump_counter.RData")
