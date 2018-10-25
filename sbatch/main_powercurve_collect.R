rm(list=ls())

intended_param <- as.matrix(expand.grid(2, c(0,.25,.5,1,2,4), 2, NA, 3, 1))
colnames(intended_param) <- c("Type", "SnR", "method", "sigma", "ksteps", "decluttered")

arg_mat <- apply(intended_param, 1, function(x){paste(x, collapse = "-")})

all_res <- vector("list", length(arg_mat))

for(i in 1:length(arg_mat)){
  tmp <- vector("list", 10000)

  files_vec <- list.files(path = paste0("/home/kevinl1/selectivemodel/sbatch/results/", arg_mat[i]),
                          full.names = T)

  #resort
  file_number <- as.numeric(sapply(files_vec, function(x){
    split <- strsplit(x, split = "/|\\.")[[1]]
    as.numeric(split[length(split)-1])
  }))
  files_vec <- files_vec[order(file_number, decreasing = F)]

  #start loading in
  counter <- 1
  for(j in 1:length(files_vec)){
    load(files_vec[j])
    tmp[counter:(counter+length(res)-1)] <- res

    counter <- counter+length(res)
  }

  len_vec <- sapply(tmp, length)
  tmp <- tmp[which(len_vec == max(len_vec))]

  tmp <- do.call(cbind, tmp)
  print(paste0("Done with iteration ", i, ": ", ncol(tmp), " columns"))
  all_res[[i]] <- tmp
}

#labeling
res <- all_res
names(res) <- arg_mat

#create the filename
string_vec <- rep(NA, ncol(intended_param))
for(i in 1:ncol(intended_param)){
  vec <- intended_param[,i]
  if(length(unique(vec)) == 1){
    string_vec[i] <- as.character(unique(vec))
  } else {
    string_vec[i] <- paste0(min(vec), "_", max(vec))
  }
}

filename <- paste0("/home/kevinl1/selectivemodel/sbatch/results/", paste0(string_vec, collapse = "-"), ".RData")
save(res, intended_param, file = filename)
