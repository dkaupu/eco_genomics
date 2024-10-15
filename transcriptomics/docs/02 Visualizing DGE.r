library(DESeq2)
library(ggplot2)
library(pheatmap)
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

###

resultsNames(dds) # Intercept, DevTemp_D22_vs_D18, FinalTemp_A33_vs_A28, FinalTemp_BASE_vs_A28

res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha = 0.05) ## pulling out DT22v18 results
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),] ## order by padj
head(res_D22vsD18) ## most significant (padj) difference of expression

############ Counts Plot ###############

## looking at counts of a specific top gene we're interested in to validate the model's working
d <- plotCounts(dds, gene = "TRINITY_DN140854_c0_g5_i2", int = (c("DevTemp", "FinalTemp")), returnData = T)
d # for the gene, counts for each sample

p <- ggplot(d, aes(x=DevTemp,y=count, 
                   color=DevTemp, shape=FinalTemp)) +
                   theme_minimal() +
                   theme(text = element_text(size=20), panel.grid.major=element_line(colour="gray")) 
p <- p + geom_point(position = position_jitter(w=0.2, h=0), size=3)
p

############# MA Plot #############
## amount of logfoldchange/DGE vs avg gene expression, so much overdispersion!
## high-density (highly expressed) genes prob for replication, metabolism, etc.
plotMA(res_D22vsD18, ylim=c(-4,4))

############# Volcano plot ##############

res_df <- as.data.frame(res_D22vsD18) ## converting DEseq results obj into data frame

##adding column to dataframe to label if significantly DGE or not
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                             "Significant", "Not Significant")

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=Significant)) +
               geom_point(alpha=0.8) + ## each point is a gene 
               scale_color_manual(values = c("slateblue","darkseagreen")) +
               labs(x="Log2 Fold Change", y="=log10 Adjusted P-value", title= "Volcano Plot") +
               theme_minimal() +
               theme(legend.position = "top") +
               geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") +
               geom_vline(xintercept=c(-1,1), linetype="dashed", color="gray") ## at least doubling of gene expression to show DGE

############## Heatmap #########

vsd <- vst(dds, blind = F)

topgenes  <- head(rownames(res_D22vsD18), 20) ## looking at 20 top genes (ordered by sf) in results file
mat <- assay(vsd)[topgenes,] ## counts matrix that's undergone variance stabilization
df <- as.data.frame(colData(dds)[,c("DevTemp","FinalTemp")])

pheatmap(mat, annotation_col=df, #color is magnitude of expression
      show_rownames= F, 
      cluster_cols = T, # col = samples
      cluster_rows = T) # rows = genes

## phenogram trying to group samples based on similar gene expression patterns 
## (row maybe same family of genes, col shows overall affect of treatment if clustered)


