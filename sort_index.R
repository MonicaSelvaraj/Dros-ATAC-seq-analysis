#This script sorts and indexes mapped ATAC reads

library(Rsamtools)

outBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/ATAC.bam"

#make an output folder for sorted bam. CHECK THIS - it ended up one index up
sortedBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_ATAC.bam"

#sort bam file
sortBam(outBAM, gsub("\\.bam", "", basename(sortedBAM)))

#index bam file
indexBam(sortedBAM)
