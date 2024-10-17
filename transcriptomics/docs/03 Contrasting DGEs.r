library(DESeq2)
library(ggplot2)
library(pheatmap)
library(eulerr)
options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/transcriptomics/")

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)
countsTableRound <- round(countsTable) 
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", 
                    header = T, stringsAsFactors = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds,
                              design = ~ DevTemp + FinalTemp)
dds <- dds[rowSums(counts(dds)>= 10) >=15,]
dds <- DESeq(dds)

##### Naming our groups here 

dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group ## make groups based on above factors, making every level of group factors variable
dds <- DESeq(dds)

resultsNames(dds) 
  # Intercept,group_D18A33_vs_D18A28,group_D18BASE_vs_D18A28,
  # group_D22A28_vs_D18A28, group_D22A33_vs_D18A28, group_D22BASE_vs_D18A28

###### 1. Comparing baseline GE between dev treatment groups ######
res_D18_BASE_D22_BASE <- results(dds, contrast=c("group","D18BASE","D22BASE"), alpha = 0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]## removing N/As
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE) ## 556 genes upregulated, 1379 downregulated

# making list of which genes in our comparisons of interest are DGE
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < 0.05,])

# looking at magnitude of up and downregulation
plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))

###### 2. Comparing GE between same dev temp, A28 final temp ######
res_D18_A28_D22_A28 <- results(dds, contrast=c("group","D18A28","D22A28"), alpha = 0.05)
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),]
head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28) ## 146 upregulated, 150 downregulated

degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < 0.05,])

plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))

########## 3. Comparing GE between same dev temp, A33 final temp ######
res_D18_A33_D22_A33 <- results(dds, contrast=c("group","D18A33","D22A33"), alpha = 0.05)
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),]
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),]
head(res_D18_A33_D22_A33)
summary(res_D18_A33_D22_A33) ## 47 upregulated, 31 downregulated

degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj < 0.05,])

plotMA(res_D18_A33_D22_A33, ylim=c(-4,4))

length(degs_D18_BASE_D22_BASE) # 1935 differentially expressed genes
length(degs_D18_A28_D22_A28) # 296
length(degs_D18_A33_D22_A33) # 78; 

########### Looking at overlap of our 3 constrasts #########
## which genes are differentially expressed in multiple contrasts?

length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A28_D22_A28)) # 107
length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A33_D22_A33)) # 44
length(intersect(degs_D18_A28_D22_A28,degs_D18_A33_D22_A33)) # 29
length(interesct(degs,))

nested_intersection <- intersect(degs_D18_BASE_D22_BASE,degs_D18_A28_D22_A28)
length(intersect(degs_D18_A33_D22_A33, nested_intersection)) ## 23 btw all groups
