set.seed(10)
n <- 24
k <- 2
dat <- rnorm(n)
fit_method <- function(x){binSegInf::binSeg_fixedSteps(x, numSteps = k)}
test_func <- selectiveModel::segment_difference
num_samp <- 2000
cores <- NA

# Rprof(tmp <- tempfile())
# selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
#                          num_samp = num_samp, ignore_jump = ignore_jump, cores =
#                            cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
# Rprof()
# summaryRprof(tmp)

##########

library(profvis)
p <- profvis({
  set.seed(10)
  selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                           num_samp = num_samp, ignore_jump = 1,
                           cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))
})
htmlwidgets::saveWidget(p, paste0("../experiments/profile_n_", n, "_k_", k, ".html"))
# browseURL("../experiments/profile.html")
