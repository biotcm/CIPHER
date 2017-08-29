#!/usr/bin/env Rscript
library(igraph)
library(pbapply)
library(reshape2)

# ppi        - protein-protein interactions
# mmi        - mesh-mesh interactions
# pheno_sim  - phenotype similarities
# pheno2gene - phenotype-gene relationships
# pheno2mesh - phenotype-mesh relationships
cipher2.leave_one_out = function (
  ppi, mmi, pheno_sim, pheno2gene, pheno2mesh,
  cluster.number = 1,
  print.timestamp = TRUE
) {

  #
  # Prepare necessary data
  #

  if (print.timestamp)
    timestamp(prefix = '', suffix = '\tpreparing necessary data...')

  gene_num <- max(ppi)
  mesh_num <- max(mmi)
  pheno_num <- nrow(pheno_sim)

  ppi_net <- graph_from_edgelist(as.matrix(ppi), directed = F)
  gene_distances <- distances(ppi_net)
  mmi_net <- graph_from_edgelist(as.matrix(mmi), directed = F)
  mesh_distances <- distances(mmi_net)

  #
  # Calculate gene-phenotype closeness
  #

  if (print.timestamp)
    timestamp(prefix = '', suffix = '\tcalculating gene2pheno closeness...')

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

  #
  # Calculate mesh-phenotype closeness
  #

  if (print.timestamp)
    timestamp(prefix = '', suffix = '\tcalculating mesh2pheno closeness...')

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

  iterate <- function (sim) {
    old_mad <- 0 # Max Abs Difference

    for (i in 1:10) {
      old_sim <- sim

      pheno_gene_score <- cor(sim, t(gene2pheno), use = 'pairwise.complete.obs')
      pheno_gene_score[is.na(pheno_gene_score)] <- 0
      sim <- 0.95 * sim + 0.05 * calc.cosine_sim(pheno_gene_score)
      pheno_mesh_score <- cor(sim, t(mesh2pheno), use = 'pairwise.complete.obs')
      pheno_mesh_score[is.na(pheno_gene_score)] <- 0
      sim <- 0.95 * sim + 0.05 * calc.cosine_sim(pheno_mesh_score)

      mad <- max(abs(sim - old_sim))
      if (abs(mad - old_mad) < 0.1) return(sim)
    }

    return(sim)
  }

  #
  # Leave-one-out test
  #

  if (print.timestamp)
    timestamp(prefix = '', suffix = '\tstart leave-one-out test...')

  leave_one_out_results <- pbapply(melt(pheno2gene), 1, function (indexes) {
    gene_index <- as.integer(indexes[1])
    pheno_index <- as.integer(indexes[2])
    pheno_genes <- pheno2gene[[indexes[2]]]
    left_pheno_genes <- setdiff(pheno_genes, gene_index)

    if (length(left_pheno_genes) == 0) {
      gene2pheno[,pheno_index] <- 0
    } else if (length(left_pheno_genes) == 1) {
      gene2pheno[,pheno_index] <- exp(-gene_distances[,left_pheno_genes]^2)
    } else {
      gene2pheno[,pheno_index] <- apply(exp(-gene_distances[,left_pheno_genes]^2), 1, sum)
    }

    # Iterate
    if (!interactive() & print.timestamp)
      timestamp(prefix = '', suffix = paste('', pheno_index, gene_index, 'iterating', sep = '\t'))
    pheno_sim <- iterate(pheno_sim)

    gene_scores <- cor(pheno_sim[,pheno_index], t(gene2pheno), use = 'pairwise.complete.obs')
    gene_scores[is.na(gene_scores)] <- -1
    percentage <- sum(quantile(gene_scores, probs = seq(0, 1, 1e-04)) > gene_scores[gene_index]) * 1e-04

    if (!interactive() & print.timestamp)
      timestamp(prefix = '', suffix = paste('', pheno_index, gene_index, percentage, sep = '\t'))

    return(percentage)
  }, cl = cluster.number)

  if (print.timestamp)
    timestamp(prefix = '', suffix = '\tcipher.leave_one_out completed!')

  return(leave_one_out_results)
}
