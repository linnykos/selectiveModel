rm(list=ls())
load("../experiments/flasso_pvalue.RData")

sat_pval
sel_pval

library(genlasso)
gen_fit <- genlasso::fusedlasso1d(y)
coef_vec <- coef(gen_fit, df = 3) #check

ylim <- range(c(y, as.numeric(res$null_samples)))
plot(y, pch = 16, ylim = ylim, col = "red")

for(i in 1:100){
  points(res$null_samples[,i], col = rgb(0,0,0,0.1))
}
