#!/usr/bin/env Rscript
library('Matrix')
source('utils/cosine_simil.R')

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

# Load data files
cat('Loading data files...\n')
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
lw <- pheno2mesh / apply(pheno2mesh, 1, max)
lw[is.nan(lw)] <- 0

# Calcualte similarities
cat('Calculating similarities...\n')
phenotype_similarities <- cosine_simil(t(t(lw) * gw))

# Save the result matrix
write.table(
  phenotype_similarities, file = '../temp/inner_pheno_simil.txt',
  quote = F, sep = "\t", row.names = F, col.names = F
)
