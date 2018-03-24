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
#' construct a basis for the entire null space of \code{mat}, which is
#' assumed to have dimension equal to the rank of \code{mat}.
#'
#' If \code{null} is \code{TRUE}, then this function samples from the nullspace.
#' Otherwise, sample from the rowspace.
#'
#' @param mat matrix
#' @param num_vec positive integer or \code{NA}
#' @param null boolean
#' @param tol small positive number
#'
#' @return a matrix of vectors, with number of rows equal to \code{ncol(mat)}
#' and number of columns equal to \code{num_vec}
.sample_matrix_space <- function(mat, num_vec = NA, null = T, tol = 1e-6){
  k <- Matrix::rankMatrix(mat)
  n <- ncol(mat)

  # form basis of either rowspace or nullspace
  if(null){
    res <- matrix(stats::rnorm(n*(n-k)), ncol = n-k)
    res <- apply(res, 2, .projection_matrix, mat = mat)
    res <- svd(t(res))$v
  } else {
    res <- svd(mat)$v
  }

  if(is.na(num_vec)) return(res)

  weights <- matrix(stats::rnorm(num_vec*ncol(res)), nrow = ncol(res))
  vec_mat <- res %*% weights

  #make the vectors orthogonal
  if(num_vec > 1){
    for(i in 2:num_vec){
      vec_mat[,i] <- .projection_matrix(vec_mat[,i], mat = t(vec_mat[,1:(i-1), drop = F]))
    }
  }

  sapply(1:ncol(vec_mat), function(x){vec_mat[,x]/.l2norm(vec_mat[,x])})
}

#' Change basis from lower dimension to higher dimension
#'
#' This function assumes \code{point} is in the coordinate basis in a lower
#' dimensional space, and translates \code{point} into a higher-dimensional
#' space parameterized by the basis matrix \code{basis_mat}, where each
#' column represents a different basis vector. Each basis vector must have
#' l2 norm of 1 and be orthogonal to one another. Additionally, we can shift this
#' higher dimensional basis by \code{center}.
#'
#' @param point vector of length k
#' @param center vector of length n
#' @param basis_mat matrix of size nxk
#'
#' @return vector of length n
.change_basis <- function(point, center, basis_mat, scaling){
  stopifnot(scaling > 0)
  stopifnot(sum(abs(t(basis_mat)%*%basis_mat) - diag(ncol(basis_mat))) < 1e-6)
  stopifnot(length(center) == nrow(basis_mat), ncol(basis_mat) == length(point))

  center + scaling * as.numeric(basis_mat %*% point)
}

#' Rotation matrix
#'
#' Creates a rotation matrix such that when multiplied by a vector, the
#' component of said vector in the plane spanned by \code{vec1} and \code{vec2}
#' are preserved, and the component in the direction of \code{vec1} is rotated
#' towards \code{vec2}.
#'
#' @param vec1 vector
#' @param vec2 vector
#'
#' @return square matrix with number columns and rows equal to the length of \code{vec1}
#' @source \url{https://math.stackexchange.com/questions/598750/finding-the-rotation-matrix-in-n-dimensions}
.rotation_matrix <- function(vec1, vec2){
  stopifnot(length(vec1) == length(vec2))
  d <- length(vec1)

  vec1 <- vec1/.l2norm(vec1)
  vec2 <- vec2/.l2norm(vec2)

  if(sum(abs(vec1 - vec2)) <= 1e-6) return(diag(d))

  vcomp <- .projection(vec2, vec1)
  vcomp <- vcomp/.l2norm(vcomp)

  val <- as.numeric(vec1%*%vec2)
  mat <- matrix(c(val, sqrt(1-val^2), -sqrt(1-val^2), val), ncol = 2)
  basis <- cbind(vec1, vcomp)

  diag(d) - vec1%*%t(vec1) - vcomp%*%t(vcomp) + basis %*% mat %*% t(basis)
}

.bijection_matrix <- function(segments){

}
