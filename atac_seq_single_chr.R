#!/usr/bin/env Rscript

library(Rsubread)
library(Rsamtools)
library(ggplot2)
library(magrittr)
library(GenomicAlignments)
library(dplyr)

#genome_2L <- "/nfs/genomes/d.melanogaster_aug_14/fasta/chrM.fa"
#indexForSubread <- ("/lab/solexa_lehmann/Dros-ATAC-seq-analysis")
#buildindex(indexForSubread, genome_2L, indexSplit = FALSE)

#read1 <- "/lab/solexa_public/LEHMANN/211103_WIGTC-MISEQ_K3LFG/FASTQ/L460_1_S1_L001_R1_001.fastq.gz"
#read2 <- "/lab/solexa_public/LEHMANN/211103_WIGTC-MISEQ_K3LFG/FASTQ/L460_1_S1_L001_R2_001.fastq.gz"

#outBAM_2L <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/ATAC_2L.bam"

#align(indexForSubread, readfile1 = read1, readfile2 = read2, output_file = outBAM_2L,
#      nthreads = 2, type = 1, unique = TRUE, maxFragLength = 2000)

sortedBAM_2L <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_ATAC_M.bam"
#sortBam(outBAM_2L, gsub("\\.bam", "", basename(sortedBAM_2L)))
#indexBam(sortedBAM_2L)
#pmapped <- propmapped(sortedBAM_2L)
#pmapped

pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/dist_mapped_reads_mito_plot_2.pdf", width = 40, height = 5)

#so here 17% of the reads mapped
idxstatsBam(sortedBAM_2L) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) +
  geom_bar(stat = "identity") + coord_flip()

dev.off()

