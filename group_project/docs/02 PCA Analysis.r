library(vcfR)
library(tidyverse)
library(qqman)
library(LEA)
library(ggplot2)
library(gridExtra)
library(ggfortify)
X11.options(type="cairo")
options(bitmapType = "cairo")

capitula <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/PNW_EU_NE_capitulummeasurements.csv")
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal.csv", row.names = "X")
meta2 <- meta[,9:11] ## only climatic variables

climate_pca <- prcomp(x = meta2,K=5, scale. = TRUE)


######## PC1 v PC2
PC1vPC2 <- ggplot(as.data.frame(climate_pca$x), 
                  aes(x = PC1, y = PC2, color=meta$region, shape=meta$continent)) +
                  geom_point(alpha=1) +
                  labs(title = "PC1 vs PC2", x = "PC1", y = "PC2", color="Region", shape="Continent") +
                  xlim(-4, 2) + ylim(-3, 2) 
                  ##geom_smooth(aes(fill = meta$region, color = meta$region), show.legend = F)
                  ##geom_polygon(aes(fill=meta$region, color=meta$region, alpha=0.3), show.legend = F)
  
autoplot(climate_pca, data=meta, colour = "region", loadings = T, loadings.label = T)

ggsave("~/projects/eco_genomics/group_project/figures/climate_pca.png", combined_plot, width=12, height=6, units="in")


