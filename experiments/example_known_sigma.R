set.seed(10)
n <- 200
k <- 4
dat <- rnorm(n, c(rep(0, n/5), rep(1, n/5), rep(0, n/5), rep(-1, n/5), rep(0, n/5)))
fit_method <- function(x){binSegInf::binSeg_fixedSteps(x, numSteps = k)}
test_func <- selectiveModel::segment_difference
num_samp <- 2000
cores <- NA

set.seed(10)
# run selective model test for known sigma, testing the first jump (aka: jump with the smallest index)
res <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
                         num_samp = num_samp, ignore_jump = 1, sigma = 1,
                         cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))

res$pval
