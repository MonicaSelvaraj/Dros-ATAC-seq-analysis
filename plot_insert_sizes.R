library(GenomicAlignments)
library(magrittr)
library(dplyr)
library(ggplot2)

sortedBAM <- "/lab/solexa_lehmann/Dros-ATAC-seq-analysis/Sorted_noDups_noMito_ATAC.bam"

#retrieve insert sizes
atacReads <- readGAlignmentPairs(sortedBAM, param = ScanBamParam(mapqFilter = 1, 
                                                                 flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
                                                                                                                                    "mapq", "isize")))
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
#print(insertSizes)

#loop to count number of fragments below and above 150
#count_less <- 0
#count_more <- 0
#total <- length(insertSizes)
#for (i in insertSizes){
#	if (i <= 150){ 
#		count_less = count_less + 1
#	}
#	else{ 
#		count_more = count_more + 1
#	}
#}
#print(count_less/total)
#print(count_more/total)
		

#this is somewhat opaque to me but ... these numbers are matching, so hopefully the paired reads are the same length?

#plot fragment length with ggplot
#fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
#                                                            Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
#                                                                                     Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
#  geom_line()

#pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/fragment_length_plot_postprocessing_noMito.pdf", width = 5, height = 5)

#fragLenPlot + theme_bw()
#looks very monosome-ish

#dev.off()

#pdf(file="/lab/solexa_lehmann/Dros-ATAC-seq-analysis/fragment_length_log_plot_postprocessing_noMito.pdf", width = 5, height = 5)
#plot on a log scale
#fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
#still very monosome-ish
#dev.off()
