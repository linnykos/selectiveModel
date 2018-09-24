rm(list=ls())
library(binseginf)
middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}

lev = 2
n = 200
mn = middle_mutation(lev=lev, n=n)
res = lapply(1:6, function(seed){
  set.seed(seed)
  y = mn + rnorm(n)
  fit = binseginf::fLasso_fixedSteps(y, numStep=2)
  binseginf:::jump_sign(fit)
})

save(res, file = "experiment_fl_check.RData")
