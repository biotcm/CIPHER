#!/usr/bin/evn Rscript
source('https://aidistan.github.io/gist/R/use.packages.R')
source('https://aidistan.github.io/gist/R/txtProgressBarETA.R')
use.packages('igraph', 'ggplot2')

#
# Load inputs
#

# Phenotype-phenotyp similarities
phenotype_similarities <- read.table('../tmp/inner_phenotype_similarity.txt')

# Protein-protein interactions
ppi <- read.table('../tmp/inner_ppi.txt')
ppi_net <- graph_from_edgelist(as.matrix(ppi), directed = F)
gene_distances <- distances(ppi_net)

# Phenotype-gene relationships (strings as keys)
phenotype_gene_relationships <- list()
fin <- file('../tmp/inner_phenotype_gene_relation.txt')
for (col in strsplit(readLines(fin), "\t")) {
  phenotype_gene_relationships[[col[1]]] <- as.numeric(c(col[2:length(col)]))
}
close(fin)

# Frequently-used numbers
gene_num <- max(ppi)
phenotype_num <- length(phenotype_similarities)

#
# Calculate gene-phenotype closeness
#

gene2phenotype_closeness <- matrix(data = 0, nrow = gene_num, ncol = phenotype_num)
for (phenotype_index in 1:phenotype_num) {
  phenotype_genes <- phenotype_gene_relationships[[as.character(phenotype_index)]]
  if (is.null(phenotype_genes)) {
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

leave_one_out_results <- c()
leave_one_out_resolution <- 1e-04

timestamp()
pb <- txtProgressBarETA(max = phenotype_num)

for (phenotype_index in 1:phenotype_num) {
  setTxtProgressBar(pb, phenotype_index)
  
  phenotype_genes <- phenotype_gene_relationships[[as.character(phenotype_index)]]
  backup_closeness <- gene2phenotype_closeness[,phenotype_index]
  
  if (is.null(phenotype_genes)) {
    next() # no known relationships to test
  } else if (length(phenotype_genes) == 1) {
    gene2phenotype_closeness[,phenotype_index] <- 0
    
    gene_score <- cor(phenotype_similarities[,phenotype_index], t(gene2phenotype_closeness))
    gene_score[is.na(gene_score)] <- 0
    
    leave_one_out_results <- c(leave_one_out_results,
                               sum(quantile(gene_score, probs = seq(0, 1, leave_one_out_resolution)) > gene_score[phenotype_genes]))
  } else {
    for (phenotype_gene_index in 1:length(phenotype_genes)) {
      new_phenotype_genes <- phenotype_genes[-phenotype_gene_index]
      if (length(new_phenotype_genes) == 1) {
        gene2phenotype_closeness[,phenotype_gene_index] <-
          exp(-gene_distances[,new_phenotype_genes]^2)
      } else {
        gene2phenotype_closeness[,phenotype_gene_index] <-
          apply(exp(-gene_distances[,new_phenotype_genes]^2), 1, sum)
      }
      
      gene_score <- cor(phenotype_similarities[,phenotype_index], t(gene2phenotype_closeness))
      gene_score[is.na(gene_score)] <- 0
      
      leave_one_out_results <- c(leave_one_out_results,
                                 sum(quantile(gene_score, probs = seq(0, 1, leave_one_out_resolution)) > gene_score[phenotype_genes]))
    }
  }
  
  gene2phenotype_closeness[,phenotype_index] <- backup_closeness
}

close(pb)
timestamp()

leave_one_out_results <- leave_one_out_results * leave_one_out_resolution

#
# Plot results
#

ggplot(data.frame(x = leave_one_out_results)) + theme_bw() +
  stat_ecdf(aes(x), geom = 'line') +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)
