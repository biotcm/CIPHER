#!/usr/bin/env Rscript
library('Matrix')

# Load dependencies
invisible(sapply(c(
  'ggplot2',
  'igraph',
  'pbapply',
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
  m <- m[, apply(m, 2, function (c) sum(c == 0) != length(c))]

  # Compute cosine
  n <- m %*% t(m) # numerator
  d <- sqrt(diag(n))
  d <- d %*% t(d) # denominator

  # Remove NaNs
  m <- n / d
  m[is.nan(m)] <- 0
  diag(m) <- 1

  return(m)
}

# For phenotype similarity calculation
calc.pheno_sim <- function () {
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

# For scoped logging printing
print.logging <- function(scope, ...) {
  cat(paste0(paste(date(), paste0('[', scope, ']'), ..., sep = "\t"), "\n"))
}

# For mesh-mesh interactions loading
read.mmi <- function (path) {
  mesh <- read.table(path, sep = "\t", quote = "", na.strings = "", stringsAsFactors = F)
  rownames(mesh) <- mesh[,2]

  mmi <- cbind(mesh[mesh[,2],1], mesh[mesh[,4],1])
  mmi <- mmi[!is.na(mmi[,1]) & !is.na(mmi[,2]),]
  return(mmi)
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

# For phenotype-mesh relationship loading
read.pheno2mesh <- function (path) {
  pheno2mesh <- list()

  tab <- read.table(path)
  for (ri in 1:nrow(tab)) {
    pheno_index <- as.character(tab[ri, 1])
    mesh_index <- as.numeric(tab[ri, 2])
    count <- tab[ri, 3]

    if (count < 5) next()
    if (is.null(pheno2mesh[[pheno_index]])) {
      pheno2mesh[[pheno_index]] <- c(mesh_index)
    } else {
      pheno2mesh[[pheno_index]] <- c(pheno2mesh[[pheno_index]], mesh_index)
    }
  }

  return(pheno2mesh)
}

# For phenotype relationship reduction
reduce.pheno_list <- function (pheno_list, to = max(as.integer(names(pheno_list)))) {
  for (n in names(pheno_list)[as.integer(names(pheno_list)) > to]) {
    pheno_list[[n]] <- NULL
  }
  return(pheno_list)
}
