#This script counts the number of mapped reads
library(Rsubread)
sortedBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/noDups_noMito.bam"
pmapped <- propmapped(sortedBAM)
pmapped
