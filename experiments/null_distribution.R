rm(list=ls())
library(selectiveModel)
trials <- 500

func <- function(i){
  print(i)
  set.seed(i)
  y <- rnorm(6)
  fit_method <- function(x){binseginf::bsfs(x, numSteps = 1)}
  test_func <- selectiveModel::segment_difference
  num_samp <- 2000
  cores <- NA

  res <- selected_model_inference(y, fit_method = fit_method, test_func = test_func,
                                  num_samp = num_samp, ignore_jump = 1,
                                  cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
  res$pval
}

# doMC::registerDoMC(cores = 14)
# res <- foreach::"%dopar%"(foreach::foreach(trial = 1:trials), func(trial))
# res <- unlist(res)
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
