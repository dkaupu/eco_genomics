## Code for analyzing RNAsed data using DESeq2

library(DESeq2)
library(ggplot2)
options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/transcriptomics/")

## importing counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                           header = T, row.names = 1)   ## 21 samples, 119428 trinity-predicted transcripts

##DESeq doesn't like decimals, gotta round
countsTableRound <- round(countsTable) 
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", ## conditions, like our metadata
                    header = T, stringsAsFactors = T, row.names = 1)

################ Exploring counts matrix ###############

colSums(countsTableRound) ## how many reads for each sample?
mean(colSums(countsTableRound)) ## gimme the mean,18454529, typically want this order # of reads

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound), #labeling it by column names of var
        cex.names = 0.5, las = 2, ylim=c(0,30000000)) # cex.names = size of labels, las=2 vertical labs
abline(h=mean(colSums(countsTableRound)), col="blue4", lwd=2) # lwd= line width

rowSums(countsTableRound) ## number of counts per gene
mean(rowSums(countsTableRound)) ## avg 3244.730 reads per sample
median(rowSums(countsTableRound)) ## median 64 reads per sample

apply(countsTableRound, 2, mean) #find mean gene expression level across columns (col=2)
                                 #avg reads/transcripts for each sample 

################## DESeq Time ##############


dds <- DESeqDataSetFromMatrix(countData = countsTableRound, 
                              colData = conds,
                              design = ~ DevTemp + FinalTemp) # what we're testing for DGE
dim(dds) # its all still there, woo!

dds <- dds[rowSums(counts(dds)>= 10) >=15,] ## 15 or more samples need 10 counts (of DGE trscrpts) per sample
nrow(dds) ## went from 119438 transcripts to 35527; # of tr with >= to 10 reads in 15 or more samples

# Running the DESeq model to test for global differential gene expression #

dds <- DESeq(dds) ## our data is normalized! accounts for diff in reads per sample

# what's our results?
resultsNames(dds) # Intercept, DevTemp_D22_vs_D18, FinalTemp_A33_vs_A28, FinalTemp_BASE_vs_A28

################# PCA Analysis ##############
## visualizing global gene expression patterns using PCA
## first need to transform data for plotting using variance stabilization, DESeq has a function!

vsd <- vst(dds, blind = F) ## normalized for PCA

pcaData <- plotPCA(vsd, intgroup = c("DevTemp","FinalTemp"),returnData = T) # group col created for dev and fin temp
percentVar <- round(100*attr(pcaData,"percentVar")) ## variance for PC1 49, PC2 15

final_temp_colors <- c("BASE"="grey", "A28" = "deeppink2", "A33" = "darkred")
shapes_choose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp)) +
            geom_point(size=5) +
            scale_shape_manual(values = shapes_choose) +
            scale_color_manual(values = final_temp_colors) +
            labs(x = paste0("PC1: ", percentVar[1], "%"),
                 y = paste0("PC2: ", percentVar[2], "%")) +
            theme_bw(base_size = 16)

p ## PC2 explains finaltemp well, PC1 loosely separating 18 and 22
## gene expression data has strong patterns of response to the treatments
