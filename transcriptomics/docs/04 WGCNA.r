## Script for analyzing and visualizing gene correlation networks

library(DESeq2)
library(ggplot2)
library(WGCNA); options(stringsAsFactors = F);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)
options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/transcriptomics/")

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", 
                    header = T, stringsAsFactors = T, row.names = 1)

traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, row.names = 1)

# Filter the matrix to just BASE data (only data we have traits for)
filtered_count_matrix_BASEonly <- countsTable[, conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE",]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)

###### Detecting Outliers ######

#detect outlier genes
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)

table(gsg$goodGenes) # false = 37235, true = 82203
table(gsg$goodSamples) # 7

# filter out bad genes
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes==T,]
dim(data_WGCNA)  

# use clustering with a tree dendrogram to identify outlier samples
htree <- hclust(dist(t(data_WGCNA)), method = "average")
plot(htree)  # even though sN2C5 looks like an outlier, goodGenes says its cool

# PCA -- outlier detection method
pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x
pca_data <- as.data.frame(pca_data) # transform to dataframe  

pca.var <- pca$sdev^2  
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0("PC1: ", pca.var.percent[1], "%"), 
       y = paste0("PC2: ", pca.var.percent[2],"%")) # still leaving everything in :P
  
##### Normalization #####

colData <- row.names(filtered_sample_metadata_BASEonly)

## running DESeq2 without knowledge of grouping

dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sample_metadata_BASEonly,
                                    design = ~1 ) # there are no specified groups
#removing low-counts transcripts
dds_WGCNA_75 <-  dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >= 6,] # the sum of the counts for each gene has to be >= 15 for at least 6 of the samples
nrow(dds_WGCNA_75) # filtered down to 29559 transcripts

dds_norm <- vst(dds_WGCNA_75) # performing variance stabilization

# get and save normalizes counts to use below
norm.counts <- assay(dds_norm) %>%
  t() # transpose it

##### Network Construction #####

## Choosing a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

## Call the network topology analysis function
sft <-  pickSoftThreshold(norm.counts,
                          powerVector = power,
                          networkType = "signed", ## focusing on genes positively correlated with e/o (regardless of direction of regulation)
                          verbose = 5)


sft.data <- sft$fitIndices

# plot to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") + # what 80% correlation line is
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") + # what 80% correlation line is
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

## WGCNA needs to pick degree of relationship; how tightly should these genes be to one another?
## scale free topology a descrition of a network; network of clusters tightly and loosely related
## if you change the strength of correlation, how does it affect clustering?
    ## as you increase power/strength of relationship with genes they're closer
## connectivity = avg number of connections a gene has, decreases as power increases

## choose a power that maximizes biological relevance while still maintaining connectivity between genes