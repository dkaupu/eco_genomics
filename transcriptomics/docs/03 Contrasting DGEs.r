library(DESeq2)
library(ggplot2)
library(pheatmap)
library(eulerr)
library(dplyr)
library(tidyr)
library(gridExtra)
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

########### Creating Euler Plot #########

length(degs_D18_BASE_D22_BASE) # 1935 differentially expressed genes
length(degs_D18_A28_D22_A28) # 296
length(degs_D18_A33_D22_A33) # 78; 

## which genes are differentially expressed in multiple contrasts?
length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A28_D22_A28)) # 107
length(intersect(degs_D18_BASE_D22_BASE,degs_D18_A33_D22_A33)) # 44
length(intersect(degs_D18_A28_D22_A28,degs_D18_A33_D22_A33)) # 29

nested_intersection <- intersect(degs_D18_BASE_D22_BASE,degs_D18_A28_D22_A28)
length(intersect(degs_D18_A33_D22_A33, nested_intersection)) ## 23 btw all groups

## calculating number of unique genes in each contrast
1935-107-44+23 # 1807 unique genes DGE at BASE b/w dev18 and dev28
296-107-29+23 # 183 unique A28 " "
78-44-29+23 # 28 unique A33 " "

107-23 # 84 uniquely shared btw BASE and A28
44-23 # 21 btw BASE & A33
29-23 # 6 btw A28 & A33

#### making the plot!
myEuler <- euler(c("BASE"=1807,"A28"=183,"A33"=28,"BASE&A28"=84,
                   "BASE&A33"=21,"A28&A33"=6,"BASE&A28&A33"=23))

plot(myEuler, lty=1:3, quantities=T)

#### Creating scatter plots ####

#### 1. Scatter plot of responses to final temp A28 when copepods develop @ 18 v 22 ####
# contrast D18_A28vBASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast = c("group","D18BASE","D18A28"), alpha=0.05))

# contrast D22_A28vBASE
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast = c("group","D22BASE","D22A28"), alpha=0.05))

# merge dataframes
res_df28 <- merge (res_D18_BASEvsA28, res_D22_BASEvsA28, by = "row.names", suffixes = c(".18",".22"))
rownames(res_df28) <- res_df28$Row.names
res_df28 <- res_df28[, -1] #removing 1st column; Row.names

# Define color mapping logic with the mutate function
res_df28 <- res_df28 %>%
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise3", ## genes significant in 18 degree contrast & upregulated in 28
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta4", # stat = test statistic, correlated to up/downregulation (if LFC was also pos, stat is pos)
    padj.22 < 0.05 & stat.22 < 0 ~ "chartreuse3",
    padj.22 < 0.05 & stat.22 > 0 ~ "coral3"
  ))

## count number of points per fill color
color_counts28 <- res_df28 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions28 <- data.frame(
  fill = c("chartreuse3","magenta4","coral3","turquoise3"),
  x_pos = c(1, 5, 0, -7.5),
  y_pos = c(-5, 0, 9, 3))

label_data28 <- merge(color_counts28, label_positions28, by = "fill")

plot28 <- ggplot(res_df28, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  geom_text(data = label_data28, aes(x = x_pos, y = y_pos, label = count, color = fill), size =5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "antiquewhite4") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "antiquewhite4") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x = "Log2FoldChange 28 vs. BASE at 18",
       y = "Log2FoldChange 28 vs. BASE at 22",
       title = "How does response to 28C vary by DevTemp?") +
  theme_minimal()
## in comparison to base, developing at 22 is upregulated

#### 2. Scatter plot of responses to final temp A33 when copepods develop @ 18 v 22 ####

res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast = c("group","D18BASE","D18A33"), alpha=0.05))

res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast = c("group","D22BASE","D22A33"), alpha=0.05))

res_df33 <- merge (res_D18_BASEvsA33, res_D22_BASEvsA33, by = "row.names", suffixes = c(".18",".22"))
rownames(res_df33) <- res_df33$Row.names
res_df33 <- res_df33[, -1] 

res_df33 <- res_df33 %>%
  mutate(fill = case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "turquoise3", 
    padj.18 < 0.05 & stat.18 > 0 ~ "magenta4",
    padj.22 < 0.05 & stat.22 < 0 ~ "chartreuse3",
    padj.22 < 0.05 & stat.22 > 0 ~ "coral3"
  ))

color_counts33 <- res_df33 %>%
  group_by(fill) %>%
  summarise(count = n())

label_positions33 <- data.frame(
  fill = c("chartreuse3","magenta4","coral3","turquoise3"),
  x_pos = c(1, 5, 0, -4),
  y_pos = c(-5, 0, 7, 2))

label_data33 <- merge(color_counts33, label_positions33, by = "fill")

plot33 <- ggplot(res_df33, aes(x = log2FoldChange.18, y = log2FoldChange.22, color = fill)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  geom_text(data = label_data33, aes(x = x_pos, y = y_pos, label = count, color = fill), size =5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "antiquewhite4") +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "antiquewhite4") +
  xlim(-10,10) + ylim(-10,10) +
  labs(x = "Log2FoldChange 33 vs. BASE at 18",
       y = "Log2FoldChange 33 vs. BASE at 22",
       title = "How does response to 33C vary by DevTemp?") +
  theme_minimal()

##### Putting the two scatter plots together in a two-panel plot! #####

combined_plot <- grid.arrange(plot28, plot33, ncol=2) #two columns, one row

#lets save it
ggsave("~/projects/eco_genomics/transcriptomics/figures/combined_scatter_plot.png", combined_plot, width = 12, height = 6)








