#This script counts the number of mapped reads
library(Rsubread)
sortedBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_ATAC.bam"
pmapped <- propmapped(sortedBAM)
pmapped
