#' Projection of vector onto another vector
#'
#' Returns the component of \code{vec1} that is orthogonal to \code{vec2}
#'
#' @param vec1 vector
#' @param vec2 vector
#' @param tol small positive number
#'
#' @return vector
.projection <- function(vec1, vec2, tol = 1e-6){
  stopifnot(length(vec1) == length(vec2))

  d <- length(vec1)
  vec2 <- vec2/.l2norm(vec2)
  as.numeric((diag(d) - vec2%*%t(vec2))%*%vec1)
}

#' Projection of vector onto rows of a matrix
#'
#' Returns the component of \code{vec} that is orthogonal to all the rows
#' of \code{mat}
#'
#' @param vec vector
#' @param mat matrix
#'
#' @return vector
.projection_matrix <- function(vec, mat){
  stopifnot(length(vec) == ncol(mat), ncol(mat) >= nrow(mat))
  n <- length(vec)
  proj_mat <- diag(n) - t(mat) %*% solve(mat %*% t(mat)) %*% mat

  as.numeric(proj_mat %*% vec)
}

#' Sample unit vectors from the null space of a matrix
#'
#' This function assumes the row space of \code{mat} is the number of
#' rows in \code{mat}. If \code{num_vec} is \code{NA}, then the function
#' construct a basis for the entire null space of \code{mat}.
#'
#' @param mat matrix
#' @param num_vec positive integer or \code{NA}.
#'
#' @return a matrix of vectors, with number of rows equal to \code{ncol(mat)}
#' and number of columns equal to \code{num_vec}
.sample_nullspace <- function(mat, num_vec = NA){
  stopifnot(num_vec > 0, num_vec %% 1 == 0)
  stopifnot(nrow(mat) + num_vec <= ncol(mat))

  n <- ncol(mat)
  vec_mat <- sapply(1:num_vec, function(x){
    vec <- stats::rnorm(n)
    .projection_matrix(vec, mat)
  })

  if(num_vec > 1){
    for(i in 2:num_vec){
      vec_mat[,i] <- .projection_matrix(vec_mat[,i], t(vec_mat[,1:(i-1),drop = F]))
    }
  }

  sapply(1:ncol(vec_mat), function(x){vec_mat[,x]/.l2norm(vec_mat[,x])})
}

.determine_rank <- function(mat, tol = 1e-6){

}

.change_basis <- function(point, center, basis_mat){

}
