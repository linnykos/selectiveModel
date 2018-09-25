rm(list=ls())
library(binseginf)

middle_mutation <- function(lev, n){
  mn <- rep(0,n)
  mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
  mn
}
true_jumps <- c(100, 140)

numSteps <- 4

set.seed(10)
counter <- 3
while(TRUE){
  set.seed(10 * counter)
  y <- middle_mutation(lev = 4, n = 200) + rnorm(200)
  fit <- binseginf::fLasso_fixedSteps(y, numSteps)
  poly <- binseginf::polyhedra(fit)
  jump_vec <- sort(binseginf::jumps(fit))

  for(i in 1:numSteps){
    if(any(abs(jump_vec[i] - true_jumps) <= 2)){
      contrast <- as.matrix(binseginf::contrast_vector(fit, i))
      sat_pval <- binseginf::pvalue(y, poly, contrast, alternative = "one.sided")

      if(sat_pval >= 0.05/2) {
        print(paste0("Counter: ", counter))
        print(paste0("Jump: ", i))
        stop()
      }
    }
  }

  counter <- counter + 1
}

##################

# investigate this particular y
library(selectiveModel)
num_samp <- 4000
burn_in <- 1000
test_func <- selectiveModel::segment_difference
fit_method <- function(x){binseginf::fLasso_fixedSteps(x, numSteps = numSteps)}
fit <- fit_method(y)
sign_mat <-  binseginf::jump_sign(fit)
direction <- sign_mat[which(sign_mat[,1] == sort(sign_mat[,1], decreasing = F)[i]),2]
sel_pval <- selectiveModel::selected_model_inference(y, fit_method = fit_method,
                                     test_func = test_func, num_samp = num_samp,
                                     ignore_jump = i, sigma = 1, cores = NA,
                                     direction = direction,
                                     verbose = F, param = list(burn_in = burn_in,
                                                               lapse = 1))$pval
