trials <- 100

func <- function(i){
  print(i)
  set.seed(i)
  y <- c(rep(0, 10), rep(5, 10)) + rnorm(20)
  fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

  selected_model_inference(y, fit_method, verbose = F, cores = NA)
}

doMC::registerDoMC(cores = 3)
res <- foreach::"%dopar%"(foreach::foreach(trial = 1:trials), func(trial))
res <- unlist(res)

plot(sort(res), seq(0, 1, length.out = length(res)), asp = T, pch = 16)
lines(c(0,1), c(0,1), col = "red", lwd = 1)
