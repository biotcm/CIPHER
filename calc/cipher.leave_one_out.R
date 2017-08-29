#!/usr/bin/env Rscript
library(igraph)
library(pbapply)
library(reshape2)

# ppi        - protein-protein interactions
# pheno_sim  - phenotype similarities
# pheno2gene - phenotype-gene relationships
cipher.leave_one_out = function (
  ppi, pheno_sim, pheno2gene,
  cluster.number = 1,
  print.timestamp = TRUE
) {

  #
  # Prepare necessary data
  #

  if (print.timestamp)
    timestamp(prefix = '', suffix = '\tpreparing necessary data...')

  gene_num <- max(ppi)
  pheno_num <- length(pheno_sim)

  ppi_net <- graph_from_edgelist(as.matrix(ppi), directed = F)
  gene_distances <- distances(ppi_net)

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
  # Leave-one-out test
  #

  if (print.timestamp)
    timestamp(prefix = '', suffix = '\tstart leave-one-out test...')

  leave_one_out_resolution <- 1e-04
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

    gene_scores <- cor(pheno_sim[,pheno_index], t(gene2pheno), use = 'pairwise.complete.obs')
    gene_scores[is.na(gene_scores)] <- -Inf
    percentage <- leave_one_out_resolution * sum(
      quantile(gene_scores, probs = seq(0, 1, leave_one_out_resolution))
      > gene_scores[gene_index]
    )

    if (!interactive() & print.timestamp)
      timestamp(prefix = '', suffix = paste('', pheno_index, gene_index, percentage, sep = '\t'))

    return(percentage)
  }, cl = cluster.number)

  if (print.timestamp)
    timestamp(prefix = '', suffix = '\tcipher.leave_one_out completed!')

  return(leave_one_out_results)
}
