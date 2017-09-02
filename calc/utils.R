#!/usr/bin/env Rscript

# Load dependencies
invisible(sapply(c(
  'ggplot2',
  'igraph',
  'pbapply',
  'reshape2'
), function (pkg) {
  if (!requireNamespace(pkg)) install.packages(pkg)
  try(attachNamespace('ggplot2'), silent = TRUE)
}))

# Load optional dependencies
if (requireNamespace('WGCNA')) attachNamespace('WGCNA')

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

# For scoped logging printing
print.logging <- function(scope, ...) {
  cat(paste0(paste(date(), paste0('[', scope, ']'), ..., sep = "\t"), "\n"))
}

# For dictionaries loading
read.dict <- function (gene_path, mesh_path, pheno_path) {
  dict = list(
    gene  = read.table(gene_path, sep = "\t"),
    mesh  = read.table(mesh_path, sep = "\t", quote = "", na.strings = "", stringsAsFactors = F),
    pheno = read.table(pheno_path, sep = "\t")
  )

  colnames(dict$gene) <- c('index', 'symbol')
  rownames(dict$gene) <- dict$gene$symbol
  colnames(dict$mesh) <- c('index', 'mindex', 'term', 'parent')
  rownames(dict$mesh) <- dict$mesh$mindex
  colnames(dict$pheno) <- c('index', 'mimNumber')
  rownames(dict$pheno) <- dict$pheno$mimNumber

  return(dict)
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

# For genome-wide leave-one-out result summary
summary.genome_wide <- function (results) {
  step <- 0.0001
  for (name in names(results)) {
    area <- sum(sapply(
      seq(0, 1, step)[-1],
      function (x, cdf) { (cdf(x-step) + cdf(x)) / 2 * step },
      ecdf(results[[name]])
    ))
    print.logging('info', name, area)
  }

  # Plot the curves
  data <- melt(results)
  colnames(data) <- c('Rank', 'Set')
  ggplot(data) + theme_bw() +
    stat_ecdf(aes(x = Rank, color = Set), geom = 'line') +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
    labs(x = NULL, y = NULL)
}
