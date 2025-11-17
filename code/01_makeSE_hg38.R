#------------------------------------------------------------------------------
# 01_makeSE_hg38.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
countToTpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#make summrized experiment
counts <- fread("data/counts_exon_uns.txt")
len <- counts$Length
gene <- counts$Geneid
counts <- counts[, -c(1:6)] %>% 
  `colnames<-`(., gsub("_hg38_Aligned.sortedByCoord.out.bam", "", colnames(counts[, -c(1:6)]))) %>%
  as.matrix %>% `rownames<-`(., gene) 
tpms <- countToTpm(counts, len)
se <- SummarizedExperiment(assays = list(counts = counts, log2tpm = log2(tpms+1)),
                           rowData = DataFrame(symbol = gene),
                           colData = DataFrame(sample =str_split(gsub("20230928R-", "", colnames(counts)), "_", simplify = T)[,1]))
se <- se[,c(1:6)]

colData(se)$sampleType <- c(rep("CNO",3),rep("PBS",3))

saveRDS(se, "rds/01_se.rds")
