#!/usr/bin/env Rscript

cosine_simil <- function (m) {
  # Remove empty columns
  m <- m[, -which(apply(m, 2, sum) == 0)]

  # Compute cosine
  # n <- forceSymmetric(m %*% t(m)) # numerator
  n <- m %*% t(m) # numerator
  d <- sqrt(diag(n))
  # d <- forceSymmetric(d %*% t(d)) # denominator
  d <- d %*% t(d) # denominator
  return(n / d)
}
