set.seed(10)
dat <- rnorm(6)
fit_method <- function(x){binSegInf::binSeg_fixedSteps(x, numSteps = 1)}
test_func <- selectiveModel::segment_difference
num_samp <- 2000
cores <- NA

Rprof(tmp <- tempfile())
selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                         num_samp = num_samp, ignore_jump = ignore_jump, cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
Rprof()
summaryRprof(tmp)

##########

library(profvis)
p <- profvis({
  set.seed(10)
  selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                           num_samp = num_samp, ignore_jump = ignore_jump, cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
})
htmlwidgets::saveWidget(p, "../experiments/profile.html")
browseURL("../experiments/profile.html")
