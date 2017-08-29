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
cat('Loading protein-protein interactions...\n')
ppi <- read.table('../temp/inner_ppi.txt')
cat('Loading phenotype similarities...\n')
load('../temp/inner_pheno_sim.Rdata')
cat('Loading phenotype-gene relationships...\n')
pheno2gene <- read.pheno2gene('../temp/inner_pheno2gene_direct.txt')
cat('Loading mesh-mesh interactions...\n')
mmi <- read.mmi('../temp/inter_mesh.txt')
cat('Loading phenotype-mesh relationships...\n')
pheno2mesh <- read.pheno2mesh('../temp/inner_pheno2mesh_freq.txt')

# !!!FOR TEST!!!
pheno_sim <- pheno_sim[1:1000, 1:1000]
pheno2gene <- reduce.pheno_list(pheno2gene, 1000)
pheno2mesh <- reduce.pheno_list(pheno2mesh, 1000)

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
}
print(areas)

# # Plot the curves
# ggplot(melt(results)) + theme_bw() +
#   stat_ecdf(aes(x = value, color = L1), geom = 'line') +
#   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)
