#!/usr/bin/env Rscript
#
# iCIPHER workflow - calculate leave-one-out results
#

# Install dependencies
source('https://aidistan.github.io/gist/R/use.packages.R')
use.packages('igraph', 'pbapply', 'ggplot2', 'reshape2')

# Load CIPHER family functions
source('cipher.leave_one_out.R')

# Load the results
if (file.exists('../temp/results.Rdata')) {
  load('../temp/results.Rdata')
} else {
  results <- list()
}

# Calculation
results$old.hprd <- cipher.leave_one_out(
  '../data/old_hprd/inner_phenotype_similarity.txt',
  '../data/old_hprd/inner_ppi.txt',
  '../data/old_hprd/inner_phenotype_gene_relation.txt',
  cluster.number = 2
)
results$old.extended <- cipher.leave_one_out(
  '../data/old_extended/inner_phenotype_similarity.txt',
  '../data/old_extended/inner_ppi.txt',
  '../data/old_extended/inner_phenotype_gene_relation.txt',
  cluster.number = 2
)

# Save the results
save(results, file = '../temp/results.Rdata')

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

# Plot the curves
ggplot(melt(results)) + theme_bw() +
  stat_ecdf(aes(x = value, color = L1), geom = 'line') +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)
