#!/usr/bin/env Rscript
source('utils.R')
source('cipher.leave_one_out.R')

# Load previous results
if (file.exists('.RData')) {
  load('.RData')
} else {
  results <- list()
}

# Run!
results$direct <- cipher.leave_one_out(
  read.table('../temp/inner_ppi.txt'),
  read.table('../temp/inner_pheno_sim.txt'),
  read.pheno2gene('../temp/inner_pheno2gene_direct.txt'),
  cluster.number = 1
)
results$extend <- cipher.leave_one_out(
  read.table('../temp/inner_ppi.txt'),
  read.table('../temp/inner_pheno_sim.txt'),
  read.pheno2gene('../temp/inner_pheno2gene_extend.txt'),
  cluster.number = 1
)

# Summary & save current results
summary.genome_wide(results)
save(results, file = '.RData')
