context("Test Gaussian")

## .transform_gaussian is correct

test_that(".transform_gaussian works", {
  set.seed(10)
  cov_mat <- diag(10)
  cov_mat[1:5,1:5] <- 0.9
  cov_mat[6:10,6:10] <- 0.9
  diag(cov_mat) <- 1
  gaussian <- .gaussian(rep(0, 10), cov_mat)

  vec1 <- rnorm(10); vec2 <- rnorm(10)
  rotation <- .rotation_matrix(vec1, vec2)
  shift <- rnorm(10)

  res <- .transform_gaussian(gaussian, shift = shift, rotation = rotation)

  expect_true(class(res) == "gaussian")
  expect_true(length(res$mean) == length(gaussian$mean))
  expect_true(all(dim(res$covariance) == dim(gaussian$covariance)))
})


############################

## .conditional_gaussian is correct

test_that(".conditional_gaussian works", {
  set.seed(10)
  cov_mat <- diag(10)
  cov_mat[1:5,1:5] <- 0.9
  cov_mat[6:10,6:10] <- 0.9
  diag(cov_mat) <- 1
  gaussian <- .gaussian(rep(0, 10), cov_mat)

  res <- .conditional_gaussian(gaussian, c(5,10))

  expect_true(class(res) == "gaussian")
  expect_true(length(res$mean) == 8)
  expect_true(all(dim(res$covariance) == 8))
})

###############################

## .sampler_truncated_gaussian is correct

test_that(".sampler_truncated_gaussian works", {
  set.seed(1)
  res <- .sampler_truncated_gaussian(.gaussian(0, 1), -1, 1)

  expect_true(length(res) == 1)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(res <= 1)
  expect_true(res >= -1)
})

test_that(".sampler_truncated_gaussian gives uniform p-values under ks test, easy case", {
  trials <- 100
  num <- 1000
  gaussian <- .gaussian(0, 1)

  pval1 <- sapply(1:trials, function(x){
    set.seed(x)
    vec1 <- sapply(1:trials, function(x){
      .sampler_truncated_gaussian(gaussian, -1, 1)
    })

    vec2 <- sapply(1:trials, function(x){
      while(TRUE){
        y <- rnorm(1)
        if(y >= -1 & y <= 1) break()
      }
      y
    })

    ks.test(vec1, vec2)$p.value
  })

  pval2 <- sapply(1:trials, function(x){
    set.seed(x)
    vec1 <- sapply(1:trials, function(x){
      .sampler_truncated_gaussian(gaussian, -1, 1)
    })

    vec2 <- sapply(1:trials, function(x){
      while(TRUE){
        y <- rnorm(1)
        if(y >= 1 & y <= 2) break()
      }
      y
    })

    ks.test(vec1, vec2)$p.value
  })

  expect_true(sum(abs(quantile(pval1) - seq(0, 1, length.out = 5))) <=
                sum(abs(quantile(pval2) - seq(0, 1, length.out = 5))))
})

test_that(".sampler_truncated_gaussian gives uniform p-values under ks test, scaled", {
  trials <- 100
  num <- 1000
  gaussian <- .gaussian(3, 4)

  pval1 <- sapply(1:trials, function(x){
    set.seed(x)
    vec1 <- sapply(1:trials, function(x){
      .sampler_truncated_gaussian(gaussian, 1, 5)
    })

    vec2 <- sapply(1:trials, function(x){
      while(TRUE){
        y <- rnorm(1, mean = 3, sd = sqrt(4))
        if(y >= 1 & y <= 5) break()
      }
      y
    })

    ks.test(vec1, vec2)$p.value
  })

  pval2 <- sapply(1:trials, function(x){
    set.seed(x)
    vec1 <- sapply(1:trials, function(x){
      .sampler_truncated_gaussian(gaussian, 1, 5)
    })

    vec2 <- sapply(1:trials, function(x){
      while(TRUE){
        y <- rnorm(1, mean = 3)
        if(y >= 1 & y <= 5) break()
      }
      y
    })

    ks.test(vec1, vec2)$p.value
  })

  expect_true(sum(abs(quantile(pval1) - seq(0, 1, length.out = 5))) <=
                sum(abs(quantile(pval2) - seq(0, 1, length.out = 5))))
})


test_that(".sampler_truncated_gaussian samples works when both bounds above mean", {
  trials <- 100
  num <- 1000
  gaussian <- .gaussian(3, 4)

  pval1 <- sapply(1:trials, function(x){
    set.seed(2*x)
    vec1 <- sapply(1:trials, function(x){
      .sampler_truncated_gaussian(gaussian, 4, 5)
    })

    vec2 <- sapply(1:trials, function(x){
      while(TRUE){
        y <- rnorm(1, mean = 3, sd = sqrt(4))
        if(y >= 4 & y <= 5) break()
      }
      y
    })

    ks.test(vec1, vec2)$p.value
  })

  pval2 <- sapply(1:trials, function(x){
    set.seed(2*x)
    vec1 <- sapply(1:trials, function(x){
      .sampler_truncated_gaussian(gaussian, 1, 5)
    })

    vec2 <- sapply(1:trials, function(x){
      while(TRUE){
        y <- rnorm(1, mean = 5, sd = sqrt(4))
        if(y >= 4 & y <= 5) break()
      }
      y
    })

    ks.test(vec1, vec2)$p.value
  })

  expect_true(sum(abs(quantile(pval1) - seq(0, 1, length.out = 5))) <=
                sum(abs(quantile(pval2) - seq(0, 1, length.out = 5))))
})


test_that(".sampler_truncated_gaussian samples works when both bounds below mean", {
  trials <- 100
  num <- 1000
  gaussian <- .gaussian(3, 4)

  pval1 <- sapply(1:trials, function(x){
    set.seed(3*x)
    vec1 <- sapply(1:trials, function(x){
      .sampler_truncated_gaussian(gaussian, 0, 2)
    })

    vec2 <- sapply(1:trials, function(x){
      while(TRUE){
        y <- rnorm(1, mean = 3, sd = sqrt(4))
        if(y >= 0 & y <= 2) break()
      }
      y
    })

    ks.test(vec1, vec2)$p.value
  })

  pval2 <- sapply(1:trials, function(x){
    set.seed(3*x)
    vec1 <- sapply(1:trials, function(x){
      .sampler_truncated_gaussian(gaussian, 0, 2)
    })

    vec2 <- sapply(1:trials, function(x){
      while(TRUE){
        y <- rnorm(1, mean = 3.5)
        if(y >= 0 & y <= 2) break()
      }
      y
    })

    ks.test(vec1, vec2)$p.value
  })

  expect_true(sum(abs(quantile(pval1) - seq(0, 1, length.out = 5))) <=
                sum(abs(quantile(pval2) - seq(0, 1, length.out = 5))))
})

