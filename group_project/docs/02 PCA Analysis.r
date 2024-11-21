library(tidyverse)
library(qqman)
library(ggplot2)
library(gridExtra)
library(ggfortify)
X11.options(type="cairo")
options(bitmapType = "cairo")

capitula <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/PNW_EU_NE_capitulummeasurements.csv")
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal19.csv", row.names = "X")
metaclim <- meta[,9:27] ## only climatic variables

climate_pca <- prcomp(x = metaclim, scale. = TRUE)


######## PC1 v PC2
### this was originally for our first 3 variables, autoplot was used next to look at loadings
PC1vPC2 <- ggplot(as.data.frame(climate_pca$x), 
                  aes(x = PC1, y = PC2, color=meta$region, shape=meta$continent)) +
                  geom_point(alpha=1) +
                  labs(title = "PC1 vs PC2", x = "PC1", y = "PC2", color="Region", shape="Continent") +
                  xlim(-4, 2) + ylim(-3, 2) 
                  ##geom_smooth(aes(fill = meta$region, color = meta$region), show.legend = F)
                  ##geom_polygon(aes(fill=meta$region, color=meta$region, alpha=0.3), show.legend = F)
  
autoplot(climate_pca, 
         data=meta, 
         colour = "region", 
         loadings = T, 
         loadings.label = T,
         loadings.label.vjust = c(),
         loadings.label.hjust = c(),
         loadings.colour = c("tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","darkslateblue",
                             "darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue"),
         loadings.label.colour = c("tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","darkslateblue",
                                    "darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue"))

## will later edit position of each label for easy reading; Temp vars colored red and Prec vars blue

## Variables in order:
## "TempM","MDiurnalR","Isotherm","TempSeason","MxTempWM","MnTempCM","TempR","MTempWtQ","MTempDQ",            
## "MTempWmQ","MTempCQ","Prec","PrecWM","PrecDM","PrecSeason","PrecWtQ","PrecDQ","PrecWmQ","PrecCQ"

ggsave("~/projects/eco_genomics/group_project/figures/climate_loadings1.png", width=6, height=4, units="in")

############################## doin' some grouping stuffs ####

meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal19.csv", row.names = "X")
meta <- plyr::arrange(meta, region)

notA = c("NEU","CEU","WEU","NE")
metaA <- replace_with_na_if(data = meta, .predicate = is.character, condition = ~.x %in% notA)
metaA <- na.omit(metaA)
write.csv(metaA, "/gpfs1/cl/pbio3990/GroupProjects/Cent_climadapt/metaA.csv")