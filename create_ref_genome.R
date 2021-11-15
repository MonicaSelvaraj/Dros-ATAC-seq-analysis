#STEP 1 - Creating a reference genome
#This script indexes the Drosophila genome and maps the forward and reverse atac reads to the genome 

#!/usr/bin/env Rscript

setwd("/lab/solexa_lehmann/Dros-ATAC-seq-analysis")

#Import libraries to use in the script
library(Rsubread)

#read genome file path
genome <- "/nfs/genomes/d.melanogaster_aug_14/fasta_whole_genome/dm6.fa"

#set working file path
indexForSubread <- ("/lab/solexa_lehmann/Dros-ATAC-seq-analysis")

#build index
buildindex(indexForSubread, genome, indexSplit = FALSE)

