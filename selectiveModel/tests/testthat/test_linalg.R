context("Test linear algebra")

## .projection is correct

test_that(".projection works", {
  set.seed(10)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".projection preserves orthogonal vectors", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)
  res2 <- .projection(res, vec2)

  expect_true(sum(abs(res - res2)) < 1e-6)
})

test_that(".projection removes parallel vectors", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- 2*vec1
  res <- .projection(vec1, vec2)

  expect_true(sum(abs(res)) < 1e-6)
})

test_that(".projection reduces the norm", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)

  expect_true(.l2norm(vec1) >= .l2norm(res))
})

################################

## .projection_matrix is correct

test_that(".projection_matrix works", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".projection_matrix returns a vector that is orthogonal to rows of mat", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  for(i in 1:nrow(mat)){
    expect_true(abs(res %*% mat[i,]) < 1e-6)
  }
})

test_that(".projection_matrix reduces the norm", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  expect_true(.l2norm(vec) >= .l2norm(res))
})

test_that(".projection_matrix preserves orthogonal vectors", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)
  res2 <- .projection_matrix(res, mat)

  expect_true(sum(abs(res - res2)) < 1e-6)
})

test_that(".projection_matrix removes parallel vectors", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- 2*mat[1,]
  res <- .projection_matrix(vec, mat)

  expect_true(sum(abs(res)) < 1e-6)
})

#########################

## .sample_matrix_space is correct

test_that(".sample_matrix_space works", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  res <- .sample_matrix_space(mat, 2)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(10,2)))
})

test_that(".sample_matrix_space gives vectors that are orthogonal to mat", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  res <- .sample_matrix_space(mat, 3)

  for(i in 1:3){
    for(j in 1:3){
      expect_true(abs(res[,i]%*%mat[j,]) < 1e-6)
    }
  }
})

test_that(".sample_matrix_space gives vectors that are orthogonal", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec_mat <- .sample_matrix_space(mat, 3)
  res <- t(vec_mat)%*%vec_mat

  expect_true(sum(abs(res - diag(3))) < 1e-6)
})

test_that(".sample_matrix_space can determine num_vec automatically", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec_mat <- .sample_matrix_space(mat)

  expect_true(all(dim(vec_mat) == c(10,7)))

  res <- t(vec_mat)%*%vec_mat
  expect_true(sum(abs(res - diag(7))) < 1e-6)
})

test_that(".sample_matrix_space samples from rowspace properly", {
  set.seed(10)
  mat <- .segments(10, c(4,8))
  res <- .sample_matrix_space(mat, null = F)

  #check to see if nothing remains when projecting against res
  for(i in 1:3){
    tmp <- .projection_matrix(mat[i,], t(res))
    expect_true(sum(abs(tmp)) < 1e-6)
  }
})

test_that(".sample_matrix_space similary gives proper outputs for rowspace", {
  set.seed(10)
  mat <- .segments(10, c(1,3,9))
  res <- .sample_matrix_space(mat, null = F)

  expect_true(all(dim(res) == c(10,4)))

  tmp <- t(res) %*% res
  expect_true(sum(abs(tmp - diag(4))) < 1e-6)
})
