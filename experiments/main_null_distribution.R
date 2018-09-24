rm(list=ls())
library(simulation)
library(binseginf)
library(selectiveModel)

paramMat <- cbind(0, 200, 500)
colnames(paramMat) <- c("SnR", "n", "trials")
numSteps <- 4

middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}
true_jumps <- c(100, 140)
test_func <- selectiveModel::segment_difference
num_samp <- 2000
burn_in <- 2000

rule_closure <- function(fit_method){
  function(vec){
    while(TRUE){
      dat <- middle_mutation(lev = vec["SnR"], n = vec["n"]) + stats::rnorm(vec["n"])
      fit <- fit_method(dat)
      bool_vec <- sapply(true_jumps, function(true_jump){
        any(which(abs(binseginf::jumps(fit) - true_jump) <= 2))
      })
      if(any(bool_vec)) break()
    }

    dat
  }
}

criterion_closure <- function(fit_method){
  function(dat, vec, y){
    fit <- fit_method(dat)
    jump_vec <- sort(binseginf::jumps(fit))

    res <- c(jump_vec, rep(NA, numSteps))
    names(res) <- c(paste0("Jump ", 1:numSteps), paste0("Pvalue ", 1:numSteps))

    for(i in 1:numSteps){
      if(abs(jump_vec[i] - true_jumps[i]) <= 2){
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

rule_bs <- rule_closure(fit_method_bs)
rule_fl <- rule_closure(fit_method_fl)
criterion_bs <- criterion_closure(fit_method_bs)
criterion_fl <- criterion_closure(fit_method_fl)

###########################

bs_res <- simulation::simulation_generator(rule = rule_bs, criterion = criterion_bs,
                                           paramMat = paramMat, trials = paramMat[,"trials"],
                                           cores = 15, as_list = F,
                                           filepath = "main_powercurve_null_tmp.RData",
                                           verbose = T)
save.image("main_powercurve_null.RData")

fl_res <- simulation::simulation_generator(rule = rule_fl, criterion = criterion_fl,
                                           paramMat = paramMat, trials = paramMat[,"trials"],
                                           cores = 15, as_list = F,
                                           filepath = "main_powercurve_null_tmp.RData",
                                           verbose = T)
save.image("main_powercurve_null.RData")
