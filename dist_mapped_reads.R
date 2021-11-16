#!/usr/bin/env Rscript
sortedBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_ATAC.bam"
library(Rsamtools)
library(ggplot2)
library(magrittr)

pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/dist_mapped_reads_plot.pdf", width = 5, height = 100)

idxstatsBam(sortedBAM) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) + 
    geom_bar(stat = "identity") + coord_flip()

dev.off()