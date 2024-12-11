library(tidyverse)
library(qqman)
library(ggplot2)
library(gridExtra)
library(ggfortify)
X11.options(type="cairo")
options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal19.csv", row.names = "X")
metaclim <- meta[,9:27] ## only climatic variables

climate_pca <- prcomp(x = metaclim, scale. = TRUE)

######## PC1 v PC2
### !!!!!! this was originally for our first 3 variables, autoplot was used next to look at loadings
##PC1vPC2 <- ggplot(as.data.frame(climate_pca$x), 
##                  aes(x = PC1, y = PC2, color=meta$region, shape=meta$continent)) +
##                  geom_point(alpha=1) +
##                  labs(title = "PC1 vs PC2", x = "PC1", y = "PC2", color="Region", shape="Continent") +
##                  xlim(-4, 2) + ylim(-3, 2) 
##                  ##geom_smooth(aes(fill = meta$region, color = meta$region), show.legend = F)
##                  ##geom_polygon(aes(fill=meta$region, color=meta$region, alpha=0.3), show.legend = F)

a
## we've coloured by region and coloured the loadings and loadings label by whether 
## they are a temp variable (red) or if they are a precipitation variable (blue)
autoplot(climate_pca, 
         data=meta, 
         colour = "region", 
         loadings = T, 
         loadings.label = T,
         loadings.label.vjust = c(-1,-1,-1,-1,-1,-1,0.5,-1,1,-1,-0.5,2,-0.2,2,-1,0.75,1,2,-0.5),
         loadings.label.hjust = c(0,0.5,0,0.5,0,0.75,0,0,-0.5,0,0.5,0,1.1,1,0,1.1,-0.5,0,0),
         loadings.colour = c("tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","darkslateblue",
                             "darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue"),
         loadings.label.colour = c("tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","darkslateblue",
                                    "darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue"))

ggsave("~/projects/eco_genomics/group_project/figures/climate_loadings1.png", width=12, height=6, units="in")

## I'm also going to be making one WITHOUT labels for ease of reading and looking at points (I couldn't figure out loadings alpha/opacity)
autoplot(climate_pca, 
         data=meta, 
         colour = "region")
         ##loadings = T,
         ##loadings.colour = c("tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","darkslateblue",
         ##                                  "darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue"),
         ##loadings.label.colour = c("tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","tomato","darkslateblue",
         ##                          "darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue","darkslateblue"))

ggsave("~/projects/eco_genomics/group_project/figures/climate_loadings3.png", width=12, height=6, units="in")   


## Variables in order:
## "TempM","MDiurnalR","Isotherm","TempSeason","MxTempWM","MnTempCM","TempR","MTempWtQ","MTempDQ", 9           
## "MTempWmQ","MTempCQ","Prec","PrecWM","PrecDM","PrecSeason","PrecWtQ","PrecDQ","PrecWmQ","PrecCQ"
         