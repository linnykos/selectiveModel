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

# source: https://onlinecourses.science.psu.edu/stat505/node/43
test_that(".conditional_gaussian performs the correct calculation", {
  gaussian <- .gaussian(mean = c(175, 71), covariance = matrix(c(550, 40, 40, 8), nrow = 2))
  res <- .conditional_gaussian(gaussian, 70)

  expect_true(abs(res$mean - 170) < 1e-6)
  expect_true(abs(res$covariance - 350) < 1e-6)
})

test_that(".conditional_gaussian works in our pipeline", {
  n <- 10
  y <- rep(0, n)
  v <-  c(1, 1, rep(0,n-2)); v <- v/.l2norm(v)
  gaussian <- .gaussian(rep(0, n), diag(n))

  line <- .line(y, v)
  polyhedra <- structure(list(gamma = matrix(c(rep(-1, 10), rep(1, 10)), nrow = 2, byrow = T),
                         u = c(-1, -1)), class = "polyhedra")
  interval <- .intersect_polyhedron_line(polyhedra, line)
  rotation <- .rotation_matrix(v, c(1, rep(0, n-1)))
  c <- y + interval[1]*v
  d <- y + interval[2]*v

  gaussian <- .transform_gaussian(gaussian, c, rotation)

  univariate <- .conditional_gaussian(gaussian, rep(0, n-1))

  expect_true(abs(univariate$covariance - 1) < 1e-6)
  expect_true(abs(univariate$mean - abs(diff(interval))/2) < 1e-6)
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

