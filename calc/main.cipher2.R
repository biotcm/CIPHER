#!/usr/bin/env Rscript
#
# CIPHER2 workflow - calculate leave-one-out results
#
source('utils.R')
source('cipher.leave_one_out.R')
source('cipher2.leave_one_out.R')

# Prepare data
if (!file.exists('../temp/inner_pheno_sim.Rdata')) {
  calc.pheno_sim()
}

# Load input data
print.logging('info', 'Loading protein-protein interactions...')
ppi <- read.table('../temp/inner_ppi.txt')
print.logging('info', 'Loading phenotype similarities...')
load('../temp/inner_pheno_sim.Rdata')
print.logging('info', 'Loading phenotype-gene relationships...')
pheno2gene <- read.pheno2gene('../temp/inner_pheno2gene_direct.txt')
print.logging('info', 'Loading mesh-mesh interactions...')
mmi <- read.mmi('../temp/inter_mesh.txt')
print.logging('info', 'Loading phenotype-mesh relationships...')
pheno2mesh <- read.pheno2mesh('../temp/inner_pheno2mesh_freq.txt')

# !!!FOR TEST!!!
pheno_sim <- pheno_sim[1:100, 1:100]
pheno2gene <- reduce.pheno_list(pheno2gene, 100)
pheno2mesh <- reduce.pheno_list(pheno2mesh, 100)

# Load the results
if (file.exists('.Rdata')) {
  load('.Rdata')
} else {
  results <- list()
}

# Run main function
results$cipher <- cipher.leave_one_out(
  ppi, pheno_sim, pheno2gene,
  cluster.number = 1
)

results$cipher2 <- cipher2.leave_one_out(
  ppi, mmi, pheno_sim, pheno2gene, pheno2mesh,
  cluster.number = 1
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
  print.logging('info', name, areas[[name]])
}

# # Plot the curves
# ggplot(melt(results)) + theme_bw() +
#   stat_ecdf(aes(x = value, color = L1), geom = 'line') +
#   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)
