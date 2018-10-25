rm(list=ls())
library(simulation, lib.loc = "/home/kevinl1/Rpackages/3.5")
library(binseginf, lib.loc = "/home/kevinl1/Rpackages/3.5")
library(selectiveModel, lib.loc = "/home/kevinl1/Rpackages/3.5")
library(data.tree, lib.loc = "/home/kevinl1/Rpackages/3.5")
library(hash, lib.loc = "/home/kevinl1/Rpackages/3.5")
library(lpSolve, lib.loc = "/home/kevinl1/Rpackages/3.5")

args <- commandArgs(trailingOnly=TRUE)
#arguments:
## - 1st: type of signal (1 = 1-jump, 2 = 2-jump)
## - 2nd: signal to noise ratio, integers from 1 to 6 for (0, 0.25, 0.5, 1, 2, 4)
## - 3th: method (1 = bs, 2 = fl)
## - 4th: sigma (either 1 or 2 or 1 or NA)
## - 5th: ksteps (2 to 4)
## - 6th: decluttered (0 = no, 1 = yes)
## - 7th: array index

paramMat <- matrix(as.numeric(args), ncol = length(args))
colnames(paramMat) <- c("Type", "SnR", "method", "sigma", "ksteps", "decluttered", "array")
paramMat[,"SnR"] <- c(0, 0.25, 0.5, 1, 2, 4)[paramMat[,"SnR"]]
paramMat[,"sigma"] <- c(1,NA)[paramMat[,"sigma"]]
args <- paramMat[1,]
print(paramMat)
paramMat_save <- paramMat

#######

# determine the seed
load("/home/kevinl1/selectivemodel/sbatch/preamble_jump.RData")
tab <- table(c(which(paramMat[,"Type"] == args["Type"]), which(paramMat[,"SnR"] == args["SnR"]),
               which(paramMat[,"method"] == args["method"]), which(paramMat[,"ksteps"] == args["ksteps"])))
idx <- as.numeric(names(tab)[which(tab == 4)])
stopifnot(length(idx) == 1)
jump_mat <- res[[idx]]

if(args["Type"] == 2) true_jumps <- c(100, 140) else true_jumps <- 160
indices <- which(apply(jump_mat, 2, function(x){any(sapply(true_jumps, function(y){abs(x-y) <= 2}))}))
stopifnot(args["array"]+1 <= length(indices))
seed_vec <- c(indices[args["array"]] : (indices[args["array"]+1]-1))
stopifnot(length(seed_vec) >= 1)

print(paste0("Start from ", min(seed_vec), " to ", max(seed_vec)))


## zz = sapply(res[1:12], function(k){ if(!is.matrix(k)) k <- matrix(k, nrow = 1); length(which(apply(k, 2, function(x){abs(x-160) <= 2})))})
## zz = sapply(res[13:48], function(k){ if(!is.matrix(k)) k <- matrix(k, nrow = 1); length(which(apply(k, 2, function(x){any(sapply(c(100,140), function(y){abs(x-y) <= 2}))})))})

######
paramMat <- paramMat_save
num_samp <- 4000
burn_in <- 1000
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
                                                          how_close = 5)$jump_vec}
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
    print(y)
    fit <- fit_method(dat)
    sign_mat <- binseginf::jump_sign(fit)
    if(vec["decluttered"] == 0){
      tmp <- unlist(lapply(true_jumps, function(x){x + c(-2:2)}))
      cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                                how_close = 0,
                                                desired_jumps = tmp)
    } else {
      cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                                how_close = 5,
                                                desired_jumps = true_jumps)
    }

    res <- rep(NA, 3*vec["ksteps"]+2)
    len <- length(cluster_list$jump_vec)
    res[1:len] <- cluster_list$jump_vec
    names(res) <- c(paste0("Jump ", 1:vec["ksteps"]), paste0("Direction ", 1:vec["ksteps"]),
                    paste0("Pvalue ", 1:vec["ksteps"]), "Fingerprint", "Seed")

    for(i in 1:len){
      if(cluster_list$target_bool[i]){
        set.seed(10*y)
        contrast <- selectiveModel:::contrast_from_cluster(cluster_list, n, i)
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

    res[3*vec["ksteps"]+1] <- sum(dat)
    res[3*vec["ksteps"]+2] <- y
    res
  }
}

if(args["method"] == 1){
  fit_method <- function(x){binseginf::bsfs(x, numSteps = args["ksteps"])}
} else {
  fit_method <- function(x){binseginf::fLasso_fixedSteps(x, numSteps = args["ksteps"])}
}

criterion <- criterion_closure(fit_method)
## set.seed(1); criterion(rule(paramMat[1,]), paramMat[1,], 1)

###########################

folder_name <- paste0("/home/kevinl1/selectivemodel/sbatch/results/", paste0(args[1:6], collapse = "-"))
dir.create(folder_name, showWarnings = FALSE)

res <- lapply(seed_vec, function(seed){
  set.seed(seed)
  criterion(rule(paramMat[1,]), paramMat[1,], seed)
})

save.image(paste0(folder_name, "/", args["array"], ".RData"))
