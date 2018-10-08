rm(list=ls())
library(simulation)
library(binseginf)

trials <- 10000
n <- 200
paramMat <- as.matrix(expand.grid(c(0,.25,.5,1,2,4), c(1), 200))
colnames(paramMat) <- c("SnR", "ksteps", "n")

middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}

edge_mutation <- function(lev, n=200){
  mn = rep(0,n)
  mn[seq(from=n-40+1, to=n)] = lev
  return(mn)
}

rule <- function(vec){
  edge_mutation(lev = vec["SnR"], n = n) + stats::rnorm(n)
}

criterion_bs <- function(dat, vec, y){
  fit <- binseginf::bsfs(dat, numSteps = vec["ksteps"])
  binseginf::jumps(fit)
}

criterion_fl <- function(dat, vec, y){
  fit <- binseginf::fLasso_fixedSteps(dat, numSteps = vec["ksteps"])
  binseginf::jumps(fit)
}

########################

bs_res <- simulation::simulation_generator(rule = rule, criterion = criterion_bs,
                                           paramMat = paramMat, trials = trials,
                                           cores = 15, as_list = F,
                                           filepath = "preamble_jump_counter_onejump_tmp.RData",
                                           verbose = T)
save.image("preamble_jump_counter_onejump.RData")

fl_res <- simulation::simulation_generator(rule = rule, criterion = criterion_fl,
                                           paramMat = paramMat, trials = trials,
                                           cores = 15, as_list = F,
                                           filepath = "experiment_jump_counter_onejump_tmp.RData",
                                           verbose = T)
save.image("preamble_jump_counter_onejump.RData")
