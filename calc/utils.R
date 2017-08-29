#!/usr/bin/env Rscript

# Load dependencies
invisible(sapply(c(
  'igraph',
  'pbapply',
  'ggplot2',
  'reshape2'
), function (pkg) {
  if (!requireNamespace(pkg))
    install.packages()

  unloadNamespace(pkg)
  attachNamespace(pkg)
}))

# For cosine similarity calculation
calc.cosine_sim <- function (m) {
  # Remove empty columns
  m <- m[, -which(apply(m, 2, sum) == 0)]

  # Compute cosine
  # n <- forceSymmetric(m %*% t(m)) # numerator
  n <- m %*% t(m) # numerator
  d <- sqrt(diag(n))
  # d <- forceSymmetric(d %*% t(d)) # denominator
  d <- d %*% t(d) # denominator
  return(n / d)
}

# For phenotype similarity calculation
calc.pheno_sim <- function () {
  library('Matrix')

  # Load dictionaries
  cat('Loading dictionaries...\n')
  dict = list(
    phenotypes = read.table('../temp/inter_pheno.txt', sep = "\t"),
    mesh_terms = read.table('../temp/inter_mesh.txt', sep = "\t", quote = "", na.strings = "", stringsAsFactors = F)
  )
  colnames(dict$phenotypes) <- c('index', 'mimNumber')
  rownames(dict$phenotypes) <- dict$phenotype$mimNumber
  colnames(dict$mesh_terms) <- c('index', 'mindex', 'term', 'parent')
  rownames(dict$mesh_terms) <- dict$mesh_terms$mindex

  # Load pheno2mesh
  cat('Loading pheno2mesh frequences...\n')
  pheno2mesh <- read.table('../temp/inner_pheno2mesh_freq.txt')
  pheno2mesh <- sparseMatrix(i = pheno2mesh[,1], j = pheno2mesh[,2], x = pheno2mesh[,3])
  pheno2mesh <- as.matrix(pheno2mesh)

  # Define some useful constants
  num_pheno <- min(nrow(dict$mesh_terms), ncol(pheno2mesh))

  # Increase hypernyms
  cat('Increasing hypernyms...\n')
  for (i in rev(1:num_pheno)) {
    if (sum(pheno2mesh[,i]) == 0) next()
    p_m <- dict$mesh_terms[i, 'parent']
    if (is.na(p_m)) next()
    p_i <- dict$mesh_terms[p_m, 'index']

    pheno2mesh[,p_i] <- pheno2mesh[,p_i] + pheno2mesh[,i] / sum(dict$mesh_terms$parent == p_m, na.rm = T)
  }

  # Calculate global weights
  cat('Calculating global weights...\n')
  gw <- log2(num_pheno / apply(pheno2mesh, 2, function (x) sum(x > 0)))
  gw[is.infinite(gw)] <- 0

  # Calculate local weights
  cat('Calculating local weights...\n')
  pheno2mesh <- pheno2mesh / apply(pheno2mesh, 1, max)
  pheno2mesh[is.nan(pheno2mesh)] <- 0

  # Calculate mesh profiles
  cat('Calculating pheno2mesh profiles...\n')
  pheno2mesh <- t(t(pheno2mesh) * gw)

  # Calcualte similarities
  cat('Calculating similarities...\n')
  pheno_sim <- calc.cosine_sim(pheno2mesh)

  # Save the result matrix
  # write.table(
  #   pheno_sim, file = '../temp/inner_pheno_sim.txt',
  #   quote = F, sep = "\t", row.names = F, col.names = F
  # )
  save(pheno_sim, file = '../temp/inner_pheno_sim.Rdata')
}

# For phenotype-gene relationship loading
read.pheno2gene <- function (path) {
  pheno2gene <- list()

  fin <- file(path)
  for (col in strsplit(readLines(fin), "\t")) {
    pheno2gene[[col[1]]] <- unique(as.numeric(col[2:length(col)]))
  }
  close(fin)

  return(pheno2gene)
}

# For phenotype-gene relationship reduction
reduce.pheno2gene <- function (pheno2gene, to = max(as.integer(names(pheno2gene)))) {
  for (n in names(pheno2gene)[as.integer(names(pheno2gene)) > to]) {
    pheno2gene[[n]] <- NULL
  }
  return(pheno2gene)
}
