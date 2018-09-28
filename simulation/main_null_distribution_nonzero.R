rm(list=ls())
library(simulation)
library(binseginf)
library(selectiveModel)

contrast_vector <- function(obj, jump.idx, n, ...){
  jump.vec <- jumps(obj, sorted = T)
  jump <- jump.vec[jump.idx]

  jump_mat <- jump_sign(obj)
  jump_sign <- jump_mat[which(jump_mat[,"Jump"] == jump), "Sign"]

  v <- binseginf:::.contrast_vector_segment(obj, jump, n)
  v * jump_sign
}


paramMat <- cbind(4, 200)
colnames(paramMat) <- c("SnR", "n")

paramMat_bs <- cbind(paramMat, c(1200))
colnames(paramMat_bs)[3] <- "trials"
paramMat_fl <- cbind(paramMat, c(1200))
colnames(paramMat_fl)[3] <- "trials"

middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}
true_jumps <- c(100, 140)

test_func_closure <- function(contrast){
  function(y, fit = NA, jump = NA){
    as.numeric(contrast %*% y)
  }
}

num_samp <- 6000
burn_in <- 4000
numSteps <- 4

rule <- function(vec){
  middle_mutation(lev = vec["SnR"], n = vec["n"]) + stats::rnorm(vec["n"])
}

criterion_closure <- function(fit_method){
  function(dat, vec, y){
    fit <- fit_method(dat)
    jump_vec <- binseginf::jumps(fit)

    res <- rep(NA, 2*numSteps)
    len <- length(jump_vec)
    names(res) <- c(paste0("Jump ", 1:numSteps),
                    paste0("Pvalue ", 1:numSteps))
    res[1:len] <- jump_vec

    if(all(true_jumps %in% jump_vec)){
      idx <- which(true_jumps %in% true_jumps)
      for(i in 1:length(idx)){
        contrast <- contrast_vector(fit, i, vec["n"])
        val <- contrast %*% middle_mutation(lev = vec["SnR"], n = vec["n"])

        null_mean <- rep(0, vec["n"])
        null_mean[which(contrast > 0)] <- val
        test_func <- test_func_closure(contrast)

        tmp <- selectiveModel::selected_model_inference(dat, fit_method = fit_method,
                                                        test_func = test_func,
                                                        null_mean = null_mean,
                                                        num_samp = num_samp,
                                                        direction = NA,
                                                        ignore_jump = i, sigma = 1,
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

# set.seed(1); criterion_bs(rule(paramMat[1,]), paramMat[1,], 1)

###########################

fl_res <- simulation::simulation_generator(rule = rule, criterion = criterion_fl,
                                           paramMat = paramMat_fl, trials = paramMat_fl[,"trials"],
                                           cores = 15, as_list = F,
                                           filepath = "main_null_distribution_nonzero_tmp.RData",
                                           verbose = T)
save.image("main_null_distribution_nonzero.RData")

bs_res <- simulation::simulation_generator(rule = rule, criterion = criterion_bs,
                                           paramMat = paramMat_bs, trials = paramMat_bs[,"trials"],
                                           cores = 15, as_list = F,
                                           filepath = "main_null_distribution_nonzero_tmp.RData",
                                           verbose = T)
save.image("main_null_distribution_nonzero.RData")

