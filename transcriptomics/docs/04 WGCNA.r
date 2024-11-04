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
save(sft.data,file = "~/projects/eco_genomics/transcriptomics/outputs/sftdata.Rda")
sft.data <- as.data.frame("~/projects/eco_genomics/transcriptomics/outputs/sftdata.Rda")

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

combined_plot <- grid.arrange(a1, a2, nrow = 2)
ggsave("~/projects/eco_genomics/transcriptomics/figures/power.png", combined_plot, width = 10, height = 8)

## WGCNA needs to pick degree of relationship; how tightly should these genes be to one another?
## scale free topology; a description of a network; network of clusters tightly and loosely related
    ## how well do these powers match a scale-free topology vs scale topology
## if you change the strength of correlation, how does it affect clustering?
    ## as you increase power/strength of relationship with genes they're closer
## connectivity = avg number of connections a gene has, decreases as power increases

## choose a power that maximizes biological relevance while still maintaining connectivity between genes
## choosing 26!

soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor # sets the temp_cor functions to use WGCNA's correlation function

norm.counts[] <- sapply(norm.counts, as.numeric) # making it numberical data table

# create network and identify modules based on the parameters given
bwnet26 <- blockwiseModules(norm.counts, 
                            maxBlockSize = 30000,
                            TOMType = "signed", #allow correlations to be only positive or pos & neg (unsigned)
                            power = soft_power,
                            mergeCutHeight = 0.25,
                            numericLabels = F,
                            randomSeed = 1234,
                            verbose = 3)

saveRDS(bwnet26, file = "~/projects/eco_genomics/transcriptomics/outputs/bwnet26.rds")

# to load bwnet file in later, use:
# bwnet26 <- readRDS("~/projects/eco_genomics/transcriptomics/outputs/bwnet26.rds")

cor <- temp_cor # resetting cor function to baseR's cor function instead of WGCNA's

###### Exploring Module Eigengenes ######

module_eigengenes <- bwnet26$MEs
head(module_eigengenes)
dim(module_eigengenes)

# get the number of genes for each module
table(bwnet26$colors) # names each module by color; pulls out color section of table

# plot the dendrogram and the module colors
plotDendroAndColors(bwnet26$dendrograms[[1]], cbind(bwnet26$unmergedColors, bwnet26$colors),
                    c("unmerged", "merged"),
                    dendroLabels = F,
                    addGuide = T,
                    hang = 0.03,
                    guideHang = 0.05)

###### Correlation of modules with traits #######

## Defining the numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

## testing for a correlation between module eigengenes and trait data
module.trait.corr <- cor(module_eigengenes, traitData, use = "p") ## use pearson's correlation

module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples) ## calculating pvals for each correlation

heatmap.data <- merge(module_eigengenes, traitData, by = "row.names") ## VIsualizing module trait association as a heatmap
head(heatmap.data)

heatmap.data <- heatmap.data %>% ## addressing error of row.names not being numeric
  column_to_rownames(var = "Row.names") # transposes
names(heatmap.data)

## Making our heatmap of correlations
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[42:44], ## changes based on number 
             y = names(heatmap.data)[1:41],  ## of eigengenes
             col = c("blue4", "skyblue", "white", "pink", "red"))





