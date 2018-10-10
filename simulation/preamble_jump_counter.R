rm(list=ls())
library(simulation)
library(binseginf)

trials <- 10000
n <- 200
paramMat <- rbind(as.matrix(expand.grid(1, c(0,.25,.5,1,2,4), c(1,2), 1)),
                  as.matrix(expand.grid(2, c(0,.25,.5,1,2,4), c(1,2), c(2:4))))
colnames(paramMat) <- c("Type", "SnR", "method", "ksteps")

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
  if(vec["Type"] == 1){
    edge_mutation(lev = vec["SnR"], n = n) + stats::rnorm(n)
  } else {
    middle_mutation(lev = vec["SnR"], n = n) + stats::rnorm(n)
  }
}

criterion <- function(dat, vec, y){
  if(vec["method"] == 1){
    fit <- binseginf::bsfs(dat, numSteps = vec["ksteps"])
  } else {
    fit <- binseginf::fLasso_fixedSteps(dat, numSteps = vec["ksteps"])
  }
  binseginf::jumps(fit)
}

########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                           paramMat = paramMat, trials = trials,
                                           cores = 15, as_list = F,
                                           filepath = "preamble_jump_tmp.RData",
                                           verbose = T)
save.image("preamble_jump.RData")
