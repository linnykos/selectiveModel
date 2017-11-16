library(selectiveModel)
n_seq <- ceiling(seq(5, 100, length.out = 10)/2)*2
time_vec <- rep(0, length(n_seq))

for(i in 1:length(time_vec)){
  cat('*')
  set.seed(i)
  n <- n_seq[i]
  y <- c(rep(0, n/2), rep(5, n/2)) + rnorm(n)
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  start <- proc.time()["elapsed"]
  res <- selected_model_inference(y, fit_method, verbose = F, cores = 10,
                                  num_samp = 100,
                                  param = list(burn_in = 3, seed = 1,
                                               time_limit = 600))
  end <- proc.time()["elapsed"]

  time_vec[i] <- end - start
  save.image("timing_test.RData")
}
