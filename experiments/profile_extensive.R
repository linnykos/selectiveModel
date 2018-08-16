n_vec <- c(4,6,9,12,16,24)

for(n in n_vec){
  print(paste0("Doing n = ", n))
  set.seed(10)
  k <- 2
  dat <- rnorm(n)
  fit_method <- function(x){binSegInf::binSeg_fixedSteps(x, numSteps = k)}
  test_func <- selectiveModel::segment_difference
  num_samp <- 2000
  cores <- NA

  library(profvis)
  p <- profvis({
    for(i in 1:10){
      set.seed(10*i)
      dat <- rnorm(n)

      tmp <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                               num_samp = num_samp, ignore_jump = 1,
                               cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
    }
  })
  htmlwidgets::saveWidget(p, paste0("../experiments/profile_n_", n, "_k_", k, ".html"))
}
