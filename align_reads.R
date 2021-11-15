#forward/reverse read file paths
read1 <- "/lab/solexa_public/LEHMANN/211103_WIGTC-MISEQ_K3LFG/FASTQ/L460_1_S1_L001_R1_001.fastq.gz"
read2 <- "/lab/solexa_public/LEHMANN/211103_WIGTC-MISEQ_K3LFG/FASTQ/L460_1_S1_L001_R2_001.fastq.gz"

#make an output path
outBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/ATAC.bam"

#align reads to generate bam file
align(indexForSubread, readfile1 = read1, readfile2 = read2, output_file = outBAM, 
      nthreads = 2, type = 1, unique = TRUE, maxFragLength = 2000)
