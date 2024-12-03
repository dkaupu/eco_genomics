library(tidyverse)
library(qqman)
library(ggplot2)
library(gridExtra)
library(ggfortify)
library(naniar)
X11.options(type="cairo")
options(bitmapType = "cairo")

meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal19.csv", row.names = "X")
meta <- plyr::arrange(meta, region)

##### We're gonna be dviding meta into three clusters of 
##### similar regions of climate based on our clim pca

## metaA = PNW and SEU
notA = c("NEU","CEU","WEU","NE")
metaA <- replace_with_na_if(data = meta, .predicate = is.character, condition = ~.x %in% notA)
metaA <- na.omit(metaA)
write.csv(metaA, "/gpfs1/cl/pbio3990/GroupProjects/Cent_climadapt/metaA.csv")
dim(metaA) ## 123 specimen


## metaB = NEU, CEU, and WEU
notB = c("PNW","SEU","NE")
metaB <- replace_with_na_if(data = meta, .predicate = is.character, condition = ~.x %in% notB)
metaB <- na.omit(metaB)
write.csv(metaB, "/gpfs1/cl/pbio3990/GroupProjects/Cent_climadapt/metaB.csv")
dim(metaB) ## 165 specimen

## metaC = NE
notC = c("PNW","SEU","NEU","CEU","WEU")
metaC <- replace_with_na_if(data = meta, .predicate = is.character, condition = ~.x %in% notC)
metaC <- na.omit(metaC)
write.csv(metaC, "/gpfs1/cl/pbio3990/GroupProjects/Cent_climadapt/metaC.csv")
dim(metaC) ## 341 specimen

### All observations add up to our total observations !!!
