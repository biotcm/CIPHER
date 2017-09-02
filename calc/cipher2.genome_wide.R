#!/usr/bin/env Rscript
library(igraph)
library(pbapply)
library(reshape2)

cipher2.genome_wide = function (
  pheno_sim, pheno2gene, gene_distances,
  gene2pheno, mesh2pheno, gene_scores, mesh_scores,
  cluster.number = 1, verbose = TRUE
) {

  #
  # Iterate function
  #

  sub_iterate <- function (pheno_sim, pheno_index) {
    sim <- pheno_sim[pheno_index,]
    old_mad <- 0
    old_sim <- sim

    for (i in 1:10) {
      gene_scores[pheno_index,] <- cor(sim, t(gene2pheno))
      gene_scores[is.na(gene_scores)] <- 0
      tmp <- apply(gene_scores, 1, function (v) sqrt(sum(v ** 2)))
      tmp <- gene_scores %*% gene_scores[pheno_index,] / tmp / tmp[pheno_index]
      sim <- 0.95 * sim + 0.05 * tmp

      mesh_scores[pheno_index,] <- cor(sim, t(mesh2pheno))
      mesh_scores[is.na(mesh_scores)] <- 0
      tmp <- apply(mesh_scores, 1, function (v) sqrt(sum(v ** 2)))
      tmp <- mesh_scores %*% mesh_scores[pheno_index,] / tmp / tmp[pheno_index]
      sim <- 0.95 * sim + 0.05 * tmp

      # Max Abs Difference
      mad <- max(abs(sim - old_sim))

      if (abs(mad - old_mad) < 0.01) {
        return(sim)
      } else {
        old_mad <- mad
        old_sim <- sim
      }
    }

    return(sim)
  }

  #
  # Leave-one-out test
  #

  if (verbose)
      print.logging('info', 'start leave-one-out test...')

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

    scores <- cor(sub_iterate(pheno_sim, pheno_index), t(gene2pheno), use = 'pairwise.complete.obs')
    scores[is.na(scores)] <- -1
    percentage <- sum(quantile(scores, probs = seq(0, 1, 1e-04)) > scores[gene_index]) * 1e-04

    if (!interactive() & verbose)
      print.logging('data', pheno_index, gene_index, percentage)

    return(percentage)
  }, cl = cluster.number)

  if (verbose)
    print.logging('info', 'cipher2.leave_one_out completed!')

  return(leave_one_out_results)
}
