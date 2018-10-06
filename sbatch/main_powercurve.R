rm(list=ls())
library(simulation, lib.loc = "/home/kevinl1/Rpackages/3.5")
library(binseginf, lib.loc = "/home/kevinl1/Rpackages/3.5")
library(selectiveModel, lib.loc = "/home/kevinl1/Rpackages/3.5")

args <- commandArgs(trailingOnly=TRUE)
#arguments:
## - 1st: type of signal (1 = 1-jump, 2 = 2-jump)
## - 2nd: signal to noise ratio, integers from 1 to 6 for (0, 0.25, 0.5, 1, 2, 4)
## - 3th: method (1 = bs, 2 = fl)
## - 4th: sigma (either 1 or 2 or 1 or NA)
## - 5th: ksteps (2 to 4)
## - 6th: decluttered (0 = no, 1 = yes)
## - 7th: trials

paramMat <- matrix(args, ncol = length(args))
colnames(paramMat) <- c("Type", "SnR", "method", "sigma", "ksteps", "decluttered", "trials")
paramMat[,"SnR"] <- c(0, 0.25, 0.5, 1, 2, 4)[paramMat[,"SnR"]]
paramMat[,"sigma"] <- c(1,NA)[paramMat[,"sigma"]]
args <- paramMat[1,]


######
num_samp <- 4000
burn_in <- 4000
n <- 200

edge_mutation <- function(lev, n=200){
  mn = rep(0,n)
  mn[seq(from=n-40+1, to=n)] = lev
  return(mn)
}

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

if(args["decluttered"] == 0){
  declutter_func <- function(x){selectiveModel::declutter(x, sign_vec = rep(1, length(x)),
                                                          how_close = 0)$jump_vec}
} else {
  declutter_func <- function(x){selectiveModel::declutter(x, sign_vec = rep(1, length(x)),
                                                          how_close = 2)$jump_vec}
}

#########

rule <- function(vec){
  if(vec["Type"] == 1){
    edge_mutation(lev = vec["SnR"], n = n) + stats::rnorm(n)
  } else {
    middle_mutation(lev = vec["SnR"], n = n) + stats::rnorm(n)
  }
}

criterion_closure <- function(fit_method){
  function(dat, vec, y){
    fit <- fit_method(dat)
    sign_mat <- binseginf::jump_sign(fit)
    if(vec["decluttered" == 0]){
      tmp <- unlist(lapply(true_jumps, function(x){x + c(-2:2)}))
      cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                                how_close = 0,
                                                desired_jumps = tmp)
    } else {
      cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                                how_close = 2,
                                                desired_jumps = true_jumps)
    }

    res <- rep(NA, 3*vec["ksteps"])
    len <- length(cluster_list$jump_vec)
    res[1:len] <- cluster_list$jump_vec
    names(res) <- c(paste0("Jump ", 1:vec["ksteps"]), paste0("Direction ", 1:vec["ksteps"]),
                    paste0("Pvalue ", 1:vec["ksteps"]))

    for(i in 1:len){
      if(cluster_list$target_bool[i]){
        set.seed(10*y)
        contrast <- contrast_from_cluster(cluster_list, n, i)
        test_func <- test_func_closure(contrast)
        if(cluster_list$sign_mat["sign:-1",i] == 0){
          direction <- 1
        } else if(cluster_list$sign_mat["sign:+1",i] == 0){
          direction <- -1
        } else {
          direction <- NA
        }

        tmp <- selectiveModel::selected_model_inference(dat, fit_method = fit_method,
                                        test_func = test_func,
                                        declutter_func = declutter_func,
                                        num_samp = num_samp,
                                        direction = direction,
                                        ignore_jump = i,
                                        sigma = vec["sigma"],
                                        verbose = F, param = list(burn_in = burn_in,
                                                                  lapse = 1))
        res[i+vec["ksteps"]] <- direction
        res[i+2*vec["ksteps"]] <- tmp$pval
      }
    }
    res
  }
}

if(args["method"] == 1){
  fit_method <- function(x){binseginf::bsfs(x, numSteps = args["ksteps"])}
} else {
  fit_method <- function(x){binseginf::fLasso_fixedSteps(x, numSteps = args["ksteps"])}
}

criterion <- criterion_closure(fit_method)

###########################

res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                           paramMat = paramMat, trials = paramMat[,"trials"],
                                           cores = 28, as_list = F, verbose = T)
save.image(paste0("/home/kevinl1/selectivemodel/sbatch/results/main_powercurve", paste0(args, collapse = "-"), ".RData"))
