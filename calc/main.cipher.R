#!/usr/bin/env Rscript
source('utils.R')
source('cipher.leave_one_out.R')

# Load input data
cat('Loading protein-protein interactions...\n')
ppi <- read.table('../data/old_hprd/inner_ppi.txt')
# ppi <- read.table('../data/old_extended/inner_ppi.txt')
cat('Loading phenotype similarities...\n')
pheno_sim <- read.table('../data/old_hprd/inner_phenotype_similarity.txt')
# pheno_sim <- read.table('../data/old_extended/inner_phenotype_similarity.txt')
cat('Loading phenotype-gene relationships...\n')
pheno2gene <- read.pheno2gene('../data/old_hprd/inner_phenotype_gene_relation.txt')
# pheno2gene <- read.pheno2gene('../data/old_extended/inner_phenotype_gene_relation.txt')

# Load the results
if (file.exists('.Rdata')) {
  load('.Rdata')
} else {
  results <- list()
}

# Run main function
results$old_hprd <- cipher.leave_one_out(
# results$old_extend <- cipher.leave_one_out(
  ppi, pheno_sim, pheno2gene,
  cluster.number = 2
)

# Save the results
save(results, file = '.Rdata')

# Calculate the area under "curve"
step <- 0.0001
areas <- list()
for (name in names(results)) {
  areas[[name]] <- sum(sapply(
    seq(0, 1, step)[-1],
    function (x, cdf) { (cdf(x-step) + cdf(x)) / 2 * step },
    ecdf(results[[name]])
  ))
}
print(areas)

# # Plot the curves
# ggplot(melt(results)) + theme_bw() +
#   stat_ecdf(aes(x = value, color = L1), geom = 'line') +
#   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)
