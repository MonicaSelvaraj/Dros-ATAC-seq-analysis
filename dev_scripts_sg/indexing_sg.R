if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")
library(Rsubread)

genome <- "dros ref genome/dmel-all-aligned-r6.42.fasta"
indexForSubread <- gsub("\\.fa$", "", genome)

buildindex(indexForSubread, genome, indexSplit = FALSE)
