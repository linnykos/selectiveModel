rm(list=ls())

args <- commandArgs(trailingOnly=TRUE)

args <- as.numeric(args)

vec <- stats::rnorm(args[1], mean = args[2])

folder_name <- paste0("/home/kevinl1/selectivemodel/sbatch/results/", paste0(args[1], collapse = "-"))

dir.create(folder_name, showWarnings = FALSE)
save.image(paste0(folder_name, "/", args[2], ".RData"))
