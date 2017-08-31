#!/usr/bin/env Rscript
#
# CIPHER2 workflow - calculate leave-one-out results
#
source('utils.R')

# Prepare data
if (!file.exists('../temp/inner_pheno_sim.Rdata')) {
  calc.pheno_sim()
}

# Load input data
ppi <- read.table('../temp/inner_ppi.txt')
load('../temp/inner_pheno_sim.Rdata')
pheno2gene <- read.pheno2gene('../temp/inner_pheno2gene_direct.txt')
mmi <- read.mmi('../temp/inter_mesh.txt')
pheno2mesh <- read.pheno2mesh('../temp/inner_pheno2mesh_freq.txt')

gene_num <- max(ppi)
mesh_num <- max(mmi)
pheno_num <- nrow(pheno_sim)

ppi_net <- graph_from_edgelist(as.matrix(ppi), directed = F)
gene_distances <- distances(ppi_net)
mmi_net <- graph_from_edgelist(as.matrix(mmi), directed = F)
mesh_distances <- distances(mmi_net)

gene2pheno <- matrix(data = 0, nrow = gene_num, ncol = pheno_num)
for (pheno_index in 1:pheno_num) {
  pheno_genes <- pheno2gene[[as.character(pheno_index)]]
  if (length(pheno_genes) == 0) {
    next()
  } else if (length(pheno_genes) == 1) {
    gene2pheno[,pheno_index] <- exp(-gene_distances[,pheno_genes]^2)
  } else {
    gene2pheno[,pheno_index] <- apply(exp(-gene_distances[,pheno_genes]^2), 1, sum)
  }
}

mesh2pheno <- matrix(data = 0, nrow = mesh_num, ncol = pheno_num)
for (pheno_index in 1:pheno_num) {
  pheno_meshs <- pheno2mesh[[as.character(pheno_index)]]
  if (length(pheno_meshs) == 0) {
    next()
  } else if (length(pheno_meshs) == 1) {
    mesh2pheno[,pheno_index] <- exp(-mesh_distances[,pheno_meshs]^2)
  } else {
    mesh2pheno[,pheno_index] <- apply(exp(-mesh_distances[,pheno_meshs]^2), 1, sum)
  }
}

timestamp() ##------ Thu Aug 31 19:02:31 2017 ------##
gene_scores <- cor(pheno_sim, t(gene2pheno), use = 'pairwise.complete.obs')
timestamp() ##------ Thu Aug 31 20:21:47 2017 ------##
mesh_scores <- cor(pheno_sim, t(mesh2pheno), use = 'pairwise.complete.obs')
timestamp() ##------ Thu Aug 31 22:20:14 2017 ------##
