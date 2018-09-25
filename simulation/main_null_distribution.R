rm(list=ls())
library(simulation)
library(binseginf)
library(selectiveModel)

paramMat <- cbind(0, 200)
colnames(paramMat) <- c("SnR", "n")

paramMat_bs <- cbind(paramMat, 3800)
colnames(paramMat_bs)[3] <- "trials"
paramMat_fl <- cbind(paramMat, 2300)
colnames(paramMat_fl)[3] <- "trials"

numSteps <- 4

middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}
true_jumps <- c(100, 140)
test_func <- selectiveModel::segment_difference
num_samp <- 4000
burn_in <- 1000

rule <- function(vec){
  middle_mutation(lev = vec["SnR"], n = vec["n"]) + stats::rnorm(vec["n"])
}

criterion_closure <- function(fit_method){
  function(dat, vec, y){
    fit <- fit_method(dat)
    jump_vec <- sort(binseginf::jumps(fit))

    res <- c(jump_vec, rep(NA, numSteps))
    names(res) <- c(paste0("Jump ", 1:numSteps), paste0("Pvalue ", 1:numSteps))

    for(i in 1:numSteps){
      if(any(abs(jump_vec[i] - true_jumps) <= 2)){
        set.seed(10*y)
        tmp <- selected_model_inference(dat, fit_method = fit_method,
                                        test_func = test_func, num_samp = num_samp,
                                        ignore_jump = i, sigma = 1, cores = NA,
                                        verbose = F, param = list(burn_in = burn_in,
                                                                  lapse = 1))
        res[i+numSteps] <- tmp$pval
      }
    }

    res
  }
}

fit_method_bs <- function(x){binseginf::bsfs(x, numSteps = numSteps)}
fit_method_fl <- function(x){binseginf::fLasso_fixedSteps(x, numSteps = numSteps)}

criterion_bs <- criterion_closure(fit_method_bs)
criterion_fl <- criterion_closure(fit_method_fl)

###########################

bs_res <- simulation::simulation_generator(rule = rule_bs, criterion = criterion_bs,
                                           paramMat = paramMat_bs, trials = paramMat_bs[,"trials"],
                                           cores = 15, as_list = F,
                                           filepath = "main_null_distribution_tmp.RData",
                                           verbose = T)
save.image("main_null_distribution.RData")

fl_res <- simulation::simulation_generator(rule = rule_fl, criterion = criterion_fl,
                                           paramMat = paramMat_fl, trials = paramMat_fl[,"trials"],
                                           cores = 15, as_list = F,
                                           filepath = "main_null_distribution_tmp.RData",
                                           verbose = T)
save.image("main_null_distribution.RData")
