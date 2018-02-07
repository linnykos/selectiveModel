.rejection_sample <- function(lower){
  while(TRUE){
    z <- stats::rnorm(1)
    if(z >= lower) break()
  }
  z
}

.bad_sampler <- function(lower){
  upper <- lower + 10

  while(TRUE){
    z <- stats::runif(1, lower, upper)

    if(0 >= lower & 0 <= upper){
      thres <- exp(-z^2/2)
    } else if(upper < 0){
      thres <- exp((upper^2 - z^2)/2)
    } else {
      thres <- exp((lower^2 - z^2)/2)
    }

    u <- stats::runif(1)
    if(u <= thres) break()
  }

  z
}

trials <- 10000
lower <- 4
rejection_samp <- sapply(1:trials, function(x){
  if(x %% floor(trials/10) == 0) print('*')
  set.seed(x)
  .rejection_sample(lower)
})
print("done")
fancy_samp <- sapply(1:trials, function(x){
  if(x %% floor(trials/10) == 0) print('*')
  set.seed(x)
  .sampler_truncated_gaussian_onesided(lower)
})
print("done2")
bad_samp <- sapply(1:trials, function(x){
  if(x %% floor(trials/10) == 0) print('*')
  set.seed(x)
  .bad_sampler(lower)
})

par(mfrow = c(1,2))
plot(sort(rejection_samp), sort(fancy_samp), asp = T, xlab = "Rejection sampling",
     ylab = "My new method", pch = 16)
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 2)

plot(sort(rejection_samp), sort(bad_samp), asp = T, xlab = "Rejection sampling",
     ylab = "My old method", pch = 16)
lines(c(-1e6, 1e6), c(-1e6, 1e6), col = "red", lwd = 2)
