#!/usr/bin/evn Rscript
library(igraph)
library(pbapply)
library(reshape2)

cipher.leave_one_out = function (
  phenotype_similarities_filepath,
  protein_protein_interactions_filepath,
  phenotype_gene_relations_filepath,
  cluster.number = 1,
  show.timestamp = TRUE
) {

  if (show.timestamp)
    timestamp()

  #
  # Load inputs
  #

  # Phenotype-phenotype similarities
  cat('Loading phenotype similarities...\n')
  phenotype_similarities <- read.table(phenotype_similarities_filepath)

  # Protein-protein interactions
  cat('Loading protein-protein interactions...\n')
  ppi <- read.table(protein_protein_interactions_filepath)
  ppi_net <- graph_from_edgelist(as.matrix(ppi), directed = F)
  gene_distances <- distances(ppi_net)

  # Phenotype-gene relationships (strings as keys)
  cat('Loading phenotype-gene relationships...\n')
  phenotype_gene_relationships <- list()
  fin <- file(phenotype_gene_relations_filepath)
  for (col in strsplit(readLines(fin), "\t")) {
    phenotype_gene_relationships[[col[1]]] <- unique(as.numeric(col[2:length(col)]))
  }
  close(fin)

  # Frequently-used numbers
  gene_num <- max(ppi)
  phenotype_num <- length(phenotype_similarities)

  #
  # Calculate gene-phenotype closeness
  #

  cat('Calculating gene-phenotype closeness...\n')
  gene2phenotype_closeness <- matrix(data = 0, nrow = gene_num, ncol = phenotype_num)

  for (phenotype_index in 1:phenotype_num) {
    phenotype_genes <- phenotype_gene_relationships[[as.character(phenotype_index)]]
    if (length(phenotype_genes) == 0) {
      next()
    } else if (length(phenotype_genes) == 1) {
      gene2phenotype_closeness[,phenotype_index] <-
        exp(-gene_distances[,phenotype_genes]^2)
    } else {
      gene2phenotype_closeness[,phenotype_index] <-
        apply(exp(-gene_distances[,phenotype_genes]^2), 1, sum)
    }
  }

  #
  # Leave-one-out test
  #

  if (show.timestamp)
    timestamp()

  leave_one_out_resolution <- 1e-04
  leave_one_out_results <- pbapply(melt(phenotype_gene_relationships), 1, function (indexes) {
    gene_index <- as.integer(indexes[1])
    phenotype_index <- as.integer(indexes[2])
    phenotype_genes <- phenotype_gene_relationships[[indexes[2]]]
    left_phenotype_genes <- setdiff(phenotype_genes, gene_index)

    if (length(left_phenotype_genes) == 0) {
      gene2phenotype_closeness[,phenotype_index] <- 0
    } else if (length(left_phenotype_genes) == 1) {
      gene2phenotype_closeness[,phenotype_index] <-
        exp(-gene_distances[,left_phenotype_genes]^2)
    } else {
      gene2phenotype_closeness[,phenotype_index] <-
        apply(exp(-gene_distances[,left_phenotype_genes]^2), 1, sum)
    }

    gene_scores <- cor(phenotype_similarities[,phenotype_index], t(gene2phenotype_closeness))
    gene_scores[is.na(gene_scores)] <- 0

    return(sum(quantile(gene_scores, probs = seq(0, 1, leave_one_out_resolution)) > gene_scores[gene_index]))
  }, cl = cluster.number)
  leave_one_out_results <- leave_one_out_results * leave_one_out_resolution

  if (show.timestamp)
    timestamp()

  return(leave_one_out_results)
}
