#!/usr/bin/env Rscript
library(microbenchmark)

m10    <- matrix(rnorm(    100), nrow=10,    ncol= 10)
m100   <- matrix(rnorm(  10000), nrow=100,   ncol= 100)
m1000  <- matrix(rnorm(1000000), nrow=1000,  ncol= 1000)
m2000  <- matrix(rnorm(4000000), nrow=2000,  ncol= 2000)

# Load WGCNA::cor
if (!requireNamespace("WGCNA")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
  install.packages("WGCNA")
}

# Load bigcor
source('libs/bigcor.R')

# Run benchmarks
microbenchmark(times = 10,
  stats::cor(m10, m10),
  WGCNA::cor(m10, m10),
  bigcor(m10, size = 10, verbose = FALSE)
)
microbenchmark(times = 10,
  stats::cor(m100, m100),
  WGCNA::cor(m100, m100),
  bigcor(m100, size = 100, verbose = FALSE)
)
microbenchmark(times = 10,
  stats::cor(m1000, m1000),
  WGCNA::cor(m1000, m1000),
  bigcor(m1000, size = 1000, verbose = FALSE)
)
microbenchmark(times = 10,
  stats::cor(m2000, m2000),
  WGCNA::cor(m2000, m2000),
  bigcor(m2000, size = 1000, verbose = FALSE)
)
