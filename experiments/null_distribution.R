rm(list=ls())
#library(selectiveModel)
trials <- 250

func <- function(i){
  write.csv(i, file = "../experiments/tmp.csv")
  print(i)
  set.seed(i)
  y <- rnorm(20)
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  res <- selected_model_inference(y, fit_method,
                                  test_func = selectiveModel::segment_difference,
                                  verbose = F, cores = NA,
                                  num_samp = 250, sigma = 1, ignore_jump = 1,
                                  param = list(time_limit = 600))
  res$pval
}

#doMC::registerDoMC(cores = 14)
#res <- foreach::"%dopar%"(foreach::foreach(trial = 1:trials), func(trial))
#res <- unlist(res)
res <- rep(NA, trials)
for(i in 1:trials){
  if(i %% floor(trials/10) == 0) cat('*')
  res[i] <- func(i)
}

plot(sort(res), seq(0, 1, length.out = length(res)), asp = T, pch = 16)
lines(c(0,1), c(0,1), col = "red", lwd = 1)

jit <- min(diff(sort(unique(res))))
res2 <- res + runif(length(res), 0, jit/10)
ks_pval <- stats::ks.test(res2, punif)

save.image("null_distribution.RData")
