#!/usr/bin/env Rscript
setwd("/lab/solexa_lehmann/Dros-ATAC-seq-analysis")

#Import libraries used by the script
library(Rsubread)
library(Rsamtools)
library(ggplot2)
library(magrittr)
library(GenomicAlignments)
library(dplyr)

#---------------------------------------------------------------------------------------
#STEP 1: Create a reference genome. 
#Notes:
#a. This step indexes the Drosophila genome. 
#b. This only needs to be done once per reference genome. 
#c. There is no visible output for this step.

#TODO: Change reference genome file path here. File path for common ref's below:
#Human genome file path: "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Greenleaf_original_ATAC_data/hg19_Genome.fa"
#Drosophila genome file path: "/nfs/genomes/d.melanogaster_aug_14/fasta_canonical_whole_genome/dm6_can.fa"

#Read genome file path
genome <- "/nfs/genomes/d.melanogaster_aug_14/fasta_canonical_whole_genome/dm6_can.fa"
#set working file path
indexForSubread <- ("/lab/solexa_lehmann/Dros-ATAC-seq-analysis")
#build index
buildindex(indexForSubread, genome_hum, indexSplit = FALSE)

#---------------------------------------------------------------------------------------
#STEP 2: Align ATAC-seq paired end reads to indexed reference genome. 
#Notes:
#a. Output: ATAC.bam, ATAC.bam.summary, ATAC.bam.indel.vcf

#TODO: Change file path to forward and reverse reads here.
#The data is usually in the following home folder: "/lab/solexa_public/LEHMANN/"

read1 <- "/lab/solexa_public/LEHMANN/211123_WIGTC-MISEQ_K56TW/FASTQ/trimmed/L460_2_S1_L001_R1_001.fastq.gz"
read2 <- "/lab/solexa_public/LEHMANN/211123_WIGTC-MISEQ_K56TW/FASTQ/trimmed/L460_2_S1_L001_R2_001.fastq.gz"

#make an output path
outBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/ATAC.bam"

#align reads to generate bam file
align(indexForSubread, readfile1 = read1, readfile2 = read2, output_file = outBAM,
      nthreads = 2, type = 1, unique = TRUE, maxFragLength = 2000)

#----------------------------------------------------------------------------------------

#STEP 3: Sort and index mapped ATAC reads
#Notes:
#Input: ATAC.bam
#Output: Sorted_ATAC.bam

#make an output path for the sorted file
#sortedBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_ATAC.bam"

#sort bam file
#sortBam(outBAM, gsub("\\.bam", "", basename(sortedBAM)))

#index bam file
#indexBam(sortedBAM)

#----------------------------------------------------------------------------------------

#Note: Steps 4, 5, and 6 are commands from the BaRC ATAC-seq pipeline.
#Instead of running the commands directly in the terminal we are going to embed the commands in this R script using the system() function. 

#----------------------------------------------------------------------------------------

#STEP 4: Remove reads with low quality score
#Notes:
#a. Input file: Sorted_ATAC.bam
#b. Output file: MAPQ30.bam 
#system("alignmentSieve -b Sorted_ATAC.bam --minMappingQuality 30 --samFlagInclude 2 -o MAPQ30.bam")

#----------------------------------------------------------------------------------------

#STEP 5: Remove duplicates
#Notes:
#a. Input file: MAPQ30.bam
#b. Output file: noDups.bam, M=MAPQ30.marked_dup_metrics.txt
#system("java -jar /usr/local/share/picard-tools/picard.jar MarkDuplicates I=MAPQ30.bam O=noDups.bam M=MAPQ30.marked_dup_metrics.txt REMOVE_DUPLICATES=true")

#----------------------------------------------------------------------------------------

#STEP 6: Remove reads mapped to mitochondria
#Notes:
#a. Input file: noDups.bam
#b. Output file: noDups.noMito.bam
#system("samtools view -h noDups.bam | grep -v chrM | samtools view -b -h -f 0x2 - | samtools sort - > noDups.noMito.bam")

#----------------------------------------------------------------------------------------

#STEP 7: Sort and index file after quality control steps
#Notes:
#a. Input file: noDups.noMito.bam
#b. Output file: Sorted_postqc_ATAC.bam

#Set file path for the output file post quality control 
#outBAM_postqc <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/noDups.noMito.bam"

#Set file path for the post quality control sorted file
#sortedBAM_postqc <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_postqc_ATAC.bam"

#sort bam file
#sortBam(outBAM_postqc, gsub("\\.bam", "", basename(sortedBAM_postqc)))

#index bam file
#indexBam(sortedBAM_postqc)

#----------------------------------------------------------------------------------------

#Step 8: Count the number of reads that mapped post quality control 

#pmapped <- propmapped(sortedBAM_postqc)
#pmapped

#----------------------------------------------------------------------------------------

#Step 9: Plot the distribution of mapped reads. 

#pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/dist_mapped_reads.pdf", width = 10, height = 10)

#idxstatsBam(sortedBAM_postqc) %>% ggplot(aes(seqnames, mapped, fill = seqnames)) +
#  geom_bar(stat = "identity") + coord_flip()

#dev.off()

#----------------------------------------------------------------------------------------

#Step 10: Plot fragment lengths

#retrieve insert sizes
#atacReads <- readGAlignmentPairs(sortedBAM_postqc, param = ScanBamParam(mapqFilter = 1, 
#    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
#        "mapq", "isize")))

#atacReads_read1 <- GenomicAlignments::first(atacReads)
#insertSizes <- abs(elementMetadata(atacReads_read1)$isize)

#plot fragment length with ggplot
#fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes,
#    Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)),
#    Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) +
#    geom_line()

#pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/fragment_lengths.pdf", width = 5, height = 5)
#fragLenPlot + theme_bw()
#dev.off()

#Fragment length plot with y axis in log scale
#pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/fragment_lengths_ylog.pdf", width = 5, height = 5)
#plot on a log scale
#fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
#dev.off()

#----------------------------------------------------------------------------------------
