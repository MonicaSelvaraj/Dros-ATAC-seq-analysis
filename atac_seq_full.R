#!/usr/bin/env Rscript
setwd("/lab/solexa_lehmann/Dros-ATAC-seq-analysis")

library(Rsubread)
library(Rsamtools)
library(ggplot2)
library(magrittr)
library(GenomicAlignments)
library(dplyr)

#genome_hum <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Greenleaf_original_ATAC_data/hg19_Genome.fa"
#indexForSubread <- ("/lab/solexa_lehmann/Dros-ATAC-seq-analysis")
#buildindex(indexForSubread, genome_hum, indexSplit = FALSE)

#read1 <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Greenleaf_original_ATAC_data/SRR891269/SRR891269_1.fastq.gz"
#read2 <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Greenleaf_original_ATAC_data/SRR891269/SRR891269_2.fastq.gz"

#outBAM_hum <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/ATAC_hum.bam"

#align(indexForSubread, readfile1 = read1, readfile2 = read2, output_file = outBAM_hum,
#      nthreads = 2, type = 1, unique = TRUE, maxFragLength = 2000)

#sortedBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_noDups_noMito_ATAC.bam"
sortedBAM_hum <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_ATAC_hum.bam"
#sortBam(outBAM_hum, gsub("\\.bam", "", basename(sortedBAM_hum)))
#indexBam(sortedBAM_hum)

#pmapped <- propmapped(sortedBAM_hum)
#pmapped

#pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/dist_mapped_reads_hum.pdf", width = 10, height = 10)

#so here 17% of the reads mapped
#idxstatsBam(sortedBAM_hum) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) +
#  geom_bar(stat = "identity") + coord_flip()

#dev.off()

#retrieve insert sizes

atacReads <- readGAlignmentPairs(sortedBAM_hum, param = ScanBamParam(mapqFilter = 1, 
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
        "mapq", "isize"), which = GRanges("chr2", IRanges(1, 63025520))))

atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)

#plot fragment length with ggplot

fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes,
    Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)),
    Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) +
    geom_line()

pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/fragment_length_plot_hum_c2.pdf", width = 5, height = 5)

fragLenPlot + theme_bw()
#looks very monosome-ish

dev.off()

#pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/fragment_length_log_hum_c2.pdf", width = 5, height = 5)
#plot on a log scale
#fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
#still very monosome-ish
#dev.off()

