rm(list=ls())
library(binseginf)
library(selectiveModel)

middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}
true_jumps <- c(100, 140)

.form_contrast <- function(jump_vec, n, location){
  jump_vec <- sort(unique(c(0, jump_vec, n)))
  contrast_vec <- rep(0, n)

  contrast_vec[(jump_vec[location]+1):jump_vec[location+1]] <- -1/(jump_vec[location+1] - (jump_vec[location]+1) + 1)
  contrast_vec[(jump_vec[location+1]+1):jump_vec[location+2]] <- 1/(jump_vec[location+2] - (jump_vec[location+1]+1) + 1)

  contrast_vec
}

n <- 200
numSteps <- 2
num_samp <- 4000
burn_in <- 1000
test_func <- selectiveModel::segment_difference
fit_method <- function(x){binseginf::fLasso_fixedSteps(x, numSteps = numSteps)}


set.seed(10)
counter <- 289
while(TRUE){
  print(counter)
  set.seed(10 * counter)
  y <- middle_mutation(lev = 4, n = n) + rnorm(n)
  fit <- fit_method(y)
  poly <- binseginf::polyhedra(fit)
  jump_vec <- sort(binseginf::jumps(fit))

  for(i in 1:numSteps){
    if(any(abs(jump_vec[i] - true_jumps) <= 2)){
      contrast <- .form_contrast(jump_vec, n, i)
      sat_pval <- binseginf::pvalue(y, poly, contrast, alternative = "one.sided")

      if(!is.nan(sat_pval) && sat_pval >= 0.1) {
        fit <- fit_method(y)
        sign_mat <-  binseginf::jump_sign(fit)
        direction <- sign_mat[which(sign_mat[,1] == sort(sign_mat[,1], decreasing = F)[i]),2]
        res <- selectiveModel::selected_model_inference(y, fit_method = fit_method,
                                                        test_func = test_func, num_samp = num_samp,
                                                        ignore_jump = i, sigma = 1, cores = NA,
                                                        direction = direction,
                                                        verbose = F, param = list(burn_in = burn_in,
                                                                                  lapse = 1),
                                                        return_samples = T)
        sel_pval <- res$pval

        stopifnot(abs(contrast - res$test_stat) <= 1e-6)

        if(sel_pval <= 0.05/2) {
          save.image("flasso_pvalue.RData")
          stop()
        }
      }
    }
  }

  counter <- counter + 1
}


