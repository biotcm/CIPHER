#!/usr/bin/env Rscript
library(Matrix)
library(igraph)

cipher2.prepare_data <- function (
  gene_path, mesh_path, pheno_path,
  pheno2gene_path, pheno2mesh_path
) {
  calc.pheno_sim <- function () {
    # Results are always stored in `mat` to save the memory

    # Load pheno2mesh frequences
    mat <- read.table(pheno2mesh_path)
    mat <- as.matrix(sparseMatrix(
      i = mat[,1], j = mat[,2], x = mat[,3],
      dims = c(nrow(dict$pheno), nrow(dict$mesh))
    ))

    # Increase hypernyms
    for (i in rev(1:ncol(mat))) {
      if (sum(mat[,i]) == 0) next()
      p_m <- dict$mesh[i, 'parent']
      if (is.na(p_m)) next()
      p_i <- dict$mesh[p_m, 'index']

      mat[,p_i] <- mat[,p_i] + mat[,i] / sum(dict$mesh$parent == p_m, na.rm = T)
    }

    # Calculate global weights
    gw <- log2(ncol(mat) / apply(mat, 2, function (x) sum(x > 0)))
    gw[is.infinite(gw)] <- 0

    # Calculate local weights
    mat <- mat / apply(mat, 1, max)
    mat[is.nan(mat)] <- 0

    # Calculate phenotype-mesh profiles
    mat <- t(t(mat) * gw)

    # Calcualte phenotype similarities
    mat <- calc.cosine_sim(mat)

    # Save the result matrix
    write.table(
      mat, file = '../temp/inner_pheno_sim.txt',
      quote = F, sep = "\t", row.names = F, col.names = F
    )

    return(mat)
  }

  calc.gene2pheno <- function () {
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

    return(gene2pheno)
  }

  calc.mesh2pheno <- function () {
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

    return(mesh2pheno)
  }

  # Load dictionaries
  print.logging('info', 'loading dictionaries...')
  dict <- read.dict(gene_path, mesh_path, pheno_path)
  pheno2gene <- read.pheno2gene(pheno2gene_path)
  pheno2mesh <- read.pheno2mesh(pheno2mesh_path)

  # Load networks
  print.logging('info', 'loading networks...')
  ppi <- read.table('../temp/inner_ppi.txt')
  ppi_net <- graph_from_edgelist(as.matrix(ppi), directed = F)
  gene_distances <- distances(ppi_net)
  mmi <- read.mmi('../temp/inter_mesh.txt')
  mmi_net <- graph_from_edgelist(as.matrix(mmi), directed = F)
  mesh_distances <- distances(mmi_net)

  # Calculate phenotype similarities
  print.logging('info', 'calculating phenotype similarities...')
  pheno_sim <- calc.pheno_sim()

  # Define common-used constants
  gene_num <- max(ppi)
  mesh_num <- max(mmi)
  pheno_num <- nrow(pheno_sim)

  # Calculate closenesses
  print.logging('info', 'calculating gene2pheno closeness...')
  gene2pheno <- calc.gene2pheno()
  gene_scores <- cor(pheno_sim, t(gene2pheno), use = 'pairwise.complete.obs')
  print.logging('info', 'calculating mesh2pheno closeness...')
  mesh2pheno <- calc.mesh2pheno()
  mesh_scores <- cor(pheno_sim, t(mesh2pheno), use = 'pairwise.complete.obs')

  # Save all data
  print.logging('info', 'saving to temp/.RData ...')
  save(
    dict, pheno2gene, pheno2mesh,
    pheno_sim, gene_distances, mesh_distances,
    gene2pheno, gene_scores, mesh2pheno, mesh_scores,
    gene_num, mesh_num, pheno_num,
    file = '../temp/.RData'
  )
  print.logging('info', 'cipher2.prepare_data completed!')
}
