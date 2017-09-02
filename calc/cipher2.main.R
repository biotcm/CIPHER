#!/usr/bin/env Rscript
source('utils.R')
source('cipher2.prepare_data.R')
source('cipher2.genome_wide.R')

# Prepare data
if (!file.exists('../temp/.RData')) {
  cipher2.prepare_data(
    '../temp/inter_gene.txt',
    '../temp/inter_mesh.txt',
    '../temp/inter_pheno.txt',
    '../temp/inner_pheno2gene_direct.txt',
    '../temp/inner_pheno2mesh_freq.txt'
  )
}
load('../temp/.RData')

# Load previous results
if (file.exists('.RData')) {
  load('.RData')
} else {
  results <- list()
}

# Run!
results$current <- cipher2.genome_wide(
  pheno_sim, pheno2gene, gene_distances,
  gene2pheno, mesh2pheno, gene_scores, mesh_scores,
  cluster.number = 1
)

# Summary & save current results
summary.genome_wide(results)
save(results, file = '.RData')
