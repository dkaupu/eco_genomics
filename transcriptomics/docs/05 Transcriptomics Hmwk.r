library(DESeq2)
library(ggplot2)
library(WGCNA); options(stringsAsFactors = F);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)
options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/transcriptomics/")

### Reading in all our background data
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", 
                    header = T, stringsAsFactors = T, row.names = 1)
traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, row.names = 1)
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE",]
filtered_count_matrix_BASEonly <- countsTable[, conds$FinalTemp == "BASE"]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)
colData <- row.names(filtered_sample_metadata_BASEonly)

gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes==T,]
dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sample_metadata_BASEonly,
                                    design = ~1 )
dds_WGCNA_75 <-  dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >= 6,] # the sum of the counts for each gene has to be >= 15 for at least 6 of the samples
dds_norm <- vst(dds_WGCNA_75) # performing variance stabilization
norm.counts <- assay(dds_norm) %>%
  t() # transpose it
norm.counts[] <- sapply(norm.counts, as.numeric) # making it numberical data table

## In my analysis I've chosen soft power thresholds of 26 and 30, both written to outputs. I'm just gonna load them in
bwnet26 <- readRDS("~/projects/eco_genomics/transcriptomics/outputs/bwnet26.rds")
bwnet30 <- readRDS("~/projects/eco_genomics/transcriptomics/outputs/bwnet30.rds")

cor <- temp_cor # resetting cor function to baseR's cor function instead of WGCNA's

############## soft_power=26 ################
##Exploring Module Eigengenes 

module_eigengenes26 <- bwnet26$MEs
head(module_eigengenes26)
dim(module_eigengenes26)

table(bwnet26$colors) # getting the number of genes for each module; names each module by color 

plotDendroAndColors(bwnet26$dendrograms[[1]], cbind(bwnet26$unmergedColors, bwnet26$colors), # plot the dendrogram and the module colors
                    c("unmerged", "merged"),
                    dendroLabels = F,
                    addGuide = T,
                    hang = 0.03,
                    guideHang = 0.05)

## Correlation of modules with traits 
## Defining the numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

## testing for a correlation between module eigengenes and trait data
module.trait.corr26 <- cor(module_eigengenes26, traitData, use = "p") ## use pearson's correlation
module.trait.corr.pvals26 <- corPvalueStudent(module.trait.corr26, nSamples) ## calculating pvals for each correlation
heatmap.data26 <- merge(module_eigengenes26, traitData, by = "row.names") ## VIsualizing module trait association as a heatmap
heatmap.data26 <- heatmap.data26 %>% ## addressing error of row.names not being numeric
  column_to_rownames(var = "Row.names") # transposes

############## soft_power=30 ################
##Exploring Module Eigengenes 

module_eigengenes30 <- bwnet30$MEs
head(module_eigengenes30)
dim(module_eigengenes30)

table(bwnet30$colors) # getting the number of genes for each module; names each module by color 


plotDendroAndColors(bwnet30$dendrograms[[1]], cbind(bwnet30$unmergedColors, bwnet30$colors), # plot the dendrogram and the module colors
                      c("unmerged", "merged"),
                      dendroLabels = F,
                      addGuide = T,
                      hang = 0.03,
                      guideHang = 0.05)

## testing for a correlation between module eigengenes and trait data
module.trait.corr30 <- cor(module_eigengenes30, traitData, use = "p") ## use pearson's correlation
module.trait.corr.pvals30 <- corPvalueStudent(module.trait.corr30, nSamples) ## calculating pvals for each correlation
heatmap.data30 <- merge(module_eigengenes30, traitData, by = "row.names") ## VIsualizing module trait association as a heatmap
heatmap.data30 <- heatmap.data30 %>% ## addressing error of row.names not being numeric
  column_to_rownames(var = "Row.names") # transposes

##### Making our two heatmaps #####

a3 <- CorLevelPlot(heatmap.data26,
             x = names(heatmap.data26)[52:54], ## changes based on number 
             y = names(heatmap.data26)[1:51],  ## of eigengenes
             col = c("royalblue4", "skyblue", "white", "pink", "red"),
             xlab("Soft Power Threshold = 26"))
a4 <- CorLevelPlot(heatmap.data30,
             x = names(heatmap.data30)[46:48], ## changes based on number 
             y = names(heatmap.data30)[1:45],  ## of eigengenes
             col = c("royalblue4", "skyblue", "white", "pink", "red"),
             xlab("Soft Power Threshold = 30"))

combined_plot <- grid.arrange(a3, a4, nrow = 1, ncol=2)

ggsave("~/projects/eco_genomics/transcriptomics/figures/combined_heatmaps.png", combined_plot, width = 20, height = 10)

## I did it again with pretty colors :P

a5 <- CorLevelPlot(heatmap.data26,
                   x = names(heatmap.data26)[52:54], ## changes based on number 
                   y = names(heatmap.data26)[1:51],  ## of eigengenes
                   col = c("olivedrab", "olivedrab3", "seashell", "violetred1", "violetred4"),
                   xlab("Soft Power Threshold = 26"))
a6 <- CorLevelPlot(heatmap.data30,
                   x = names(heatmap.data30)[46:48], ## changes based on number 
                   y = names(heatmap.data30)[1:45],  ## of eigengenes
                   col = c("olivedrab", "olivedrab3", "seashell", "violetred1", "violetred4"),
                   xlab("Soft Power Threshold = 30"))

combined_plot <- grid.arrange(a5, a6, nrow = 1, ncol=2)
