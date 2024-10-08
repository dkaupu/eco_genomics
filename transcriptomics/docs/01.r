## Code for analyzing RNAsed data using DESeq2

library(DESeq2)
library(ggplot2)

setwd("~/projects/eco_genomics/transcriptomics/")

## importing counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                           header = T, row.names = 1)   ## 21 samples, 119428 trinity-predicted transcripts

##DESeq doesn't like decimals, gotta round
countsTableRound <- round(countsTable) 
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", ## conditions
                    header = T, stringsAsFactors = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds,
                              design = ~ DevTemp + FinalTemp)
