rm(list=ls())
library(simulation)
library(binseginf)
library(selectiveModel)

paramMat <- cbind(c(.25,.5,1,2,4), 200)
colnames(paramMat) <- c("SnR", "n")

paramMat_bs <- cbind(paramMat, c(4400, 1900, 800, 500, 400))
colnames(paramMat_bs)[3] <- "trials"
paramMat_fl <- cbind(paramMat, c(2500, 1250, 650, 450, 350))
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
num_samp <- 4000
burn_in <- 1000
numSteps <- 4

rule <- function(vec){
  middle_mutation(lev = vec["SnR"], n = vec["n"]) + stats::rnorm(vec["n"])
}

criterion_closure <- function(fit_method){
  function(dat, vec, y){
    dat <- rule(paramMat[4,])
    fit <- fit_method(dat)
    sign_mat <- binseginf::jump_sign(fit)
    cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                              how_close = 2,
                                              desired_jumps = true_jumps)
    cluster_list

    res <- rep(NA, 3*numSteps)
    len <- length(cluster_list$jump_vec)
    res[1:len] <- cluster_list$jump_vec
    names(res) <- c(paste0("Jump ", 1:numSteps), paste0("Direction ", 1:numSteps),
                    paste0("Pvalue ", 1:numSteps))

    for(i in 1:len){
      if(cluster_list$target_bool[i]){
        set.seed(10*y)
        contrast <- contrast_from_cluster(cluster_list, vec["n"], i)
        test_func <- test_func_closure(contrast)
        if(cluster_list$sign_mat["sign:-1",i] == 0){
          direction <- 1
        } else if(cluster_list$sign_mat["sign:+1",i] == 0){
          direction <- -1
        } else {
          direction <- NA
        }

        tmp <- selectiveModel::selected_model_inference(dat, fit_method = fit_method,
                                        test_func = test_func, num_samp = num_samp,
                                        direction = direction,
                                        ignore_jump = i, sigma = 1, cores = NA,
                                        verbose = F, param = list(burn_in = burn_in,
                                                                  lapse = 1))
        res[i+numSteps] <- direction
        res[i+2*numSteps] <- tmp$pval
      }
    }

    res
  }
}

fit_method_bs <- function(x){binseginf::bsfs(x, numSteps = numSteps)}
fit_method_fl <- function(x){binseginf::fLasso_fixedSteps(x, numSteps = numSteps)}

criterion_bs <- criterion_closure(fit_method_bs)
criterion_fl <- criterion_closure(fit_method_fl)

# criterion_bs(rule(paramMat[4,]), paramMat[4,], 1)

###########################

bs_res <- simulation::simulation_generator(rule = rule, criterion = criterion_bs,
                                           paramMat = paramMat_bs, trials = paramMat_bs[,"trials"],
                                           cores = 15, as_list = F,
                                           filepath = "main_powercurve_knownsigma_tmp.RData",
                                           verbose = T)
save.image("main_powercurve_knownsigma.RData")

fl_res <- simulation::simulation_generator(rule = rule, criterion = criterion_fl,
                                           paramMat = paramMat_fl, trials = paramMat_fl[,"trials"],
                                           cores = 15, as_list = F,
                                           filepath = "main_powercurve_knownsigma_tmp.RData",
                                           verbose = T)
save.image("main_powercurve_knownsigma.RData")
