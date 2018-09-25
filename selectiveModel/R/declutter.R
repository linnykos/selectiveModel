.declutter <- function(jump_vec, sign_vec, how_close = 3,
                       desired_jumps = NA){
  stopifnot(length(jump_vec) == length(sign_vec), all(sign_vec %in% c(-1,1)))

  # first sort
  idx <- order(jump_vec, decreasing = F)
  jump_vec <- jump_vec[idx]; sign_vec <- sign_vec[idx]

  # cluster
  diff_vec <- diff(jump_vec)
  bool_vec <- (diff_vec <= how_close)
  segments <- .consecutive_boolean(bool_vec)

  # locate desired jumps if there are any
  jump_list <- lapply(1:nrow(segments), function(x){
    jump_vec[segments[x,1]:segments[x,2]]
  })
  if(!any(is.na(desired_jumps))){
    idx_vec <- sapply(jump_list, function(x){
      any(sapply(desired_jumps, function(k){any(abs(x - k) <= how_close)}))
    })
  } else {
    idx_vec <- rep(NA, length(jump_list))
  }

  # collect jumps and signs
  jump_median <- sapply(jump_list, function(x){
    floor(median(x))
  })

  sign_mat <- sapply(1:nrow(segments), function(x){
    tmp <- sign_vec[segments[x,1]:segments[x,2]]
    c(length(which(tmp == -1)), length(which(tmp == 1)))
  })

  colnames(sign_mat) <- as.character(jump_median)
  rownames(sign_mat) <- c("sign:-1", "sign:1")

  list(jump_vec = jump_median,
       sign_mat = sign_mat,
       target_bool = idx_vec)
}

.consecutive_boolean <- function(vec){
  stopifnot(length(vec) > 0, all(vec %in% c(TRUE, FALSE)))

  idx <- sort(which(vec))
  mat <- matrix(NA, ncol = 2, nrow = 0)

  if(length(idx) > 1){
    idx_plus <- idx+1
    idx_minus <- idx-1

    start_position <- setdiff(idx, idx_plus)
    end_position <- setdiff(idx, idx_minus)

    mat <- rbind(mat, cbind(start_position, end_position + 1))
  }

  # add the falses back in
  idx <- 1:(length(vec)+1)
  if(nrow(mat) >= 1){
    for(i in idx){
      if(length(intersect(which(i >= mat[,1]), which(i <= mat[,2]))) == 0){
        mat <- rbind(mat, rep(i, 2))
      }
    }
  } else {
    mat <- cbind(idx, idx)
  }

  # sort
  mat <- mat[order(mat[,1]),, drop = F]

  mat
}

####

.contrast_from_cluster <- function(cluster_list, n, location){
  stopifnot(location <= length(cluster_list))
  stopifnot(n > 0, n %% 1 == 0, location > 0, location %% 1 == 0)
  stopifnot(all(sort(cluster_list$jump_vec) == cluster_list$jump_vec))

  jump_vec <- c(0, cluster_list$jump_vec, n-1)
  contrast_vec <- rep(0, n)

  contrast_vec[(jump_vec[location]+1):jump_vec[location+1]] <- -1/(jump_vec[location+1] - (jump_vec[location]+1) + 1)
  contrast_vec[(jump_vec[location+1]+1):jump_vec[location+2]] <- 1/(jump_vec[location+2] - (jump_vec[location+1]+1) + 1)

  contrast_vec
}
