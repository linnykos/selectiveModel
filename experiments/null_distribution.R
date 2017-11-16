trials <- 1000

func <- function(i){
  write.csv(i, file = "../experiments/tmp.csv")
  set.seed(i)
  y <- c(rep(0, 5), rep(5, 5)) + rnorm(10)
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  first_val <- function(y, fit, ...){y[1]}

  selected_model_inference(y, fit_method, test_func = first_val, verbose = F, cores = NA,
                           num_samp = 500,
                           param = list(burn_in = 3, seed = 1, time_limit = 600)
                           #, sample_method = "rejection")
                           )
}

doMC::registerDoMC(cores = 3)
res <- foreach::"%dopar%"(foreach::foreach(trial = 1:trials), func(trial))
res <- unlist(res)

plot(sort(res), seq(0, 1, length.out = length(res)), asp = T, pch = 16)
lines(c(0,1), c(0,1), col = "red", lwd = 1)
