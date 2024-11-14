library(raster)
library(geodata)
library(tidyverse)
library(dplyr)
## pull out BIO1, BIO7, BIO12 for our latitudes and longitudes
setwd("C:/Users/dakot/")
meta <- read.csv("C:/Users/dakot/Downloads/meta.csv")

clim <- geodata::worldclim_global(var = 'bio', res = 10, download = T, path='data')

TempM <- clim$wc2.1_10m_bio_1
TempR <- clim$wc2.1_10m_bio_7
Prec <- clim$wc2.1_10m_bio_12

id <- meta$id
lats <- meta$latitude
lons <- meta$longitude

samples <- data.frame(lons, lats, row.names=id)
temp.data <- samples 
temp.data$TempM <- terra::extract(x = TempM, y = data.frame(samples, row.names = id))
temp.data$TempR <- terra::extract(TempR, samples)
temp.data$Prec <- terra::extract(Prec, samples)

## write.csv(temp.data, "C:/Users/dakot/Downloads/metfinal.csv")
## temp.data <- read.csv("C:/Users/dakot/Downloads/metfinal.csv", row.names = "X")

temp.data_ID <- cbind.data.frame(samples, metfinal)
myVars = c("lons","lats","TempM","TempR","Prec")
temp.data_ID <- temp.data_ID[,c(myVars)]
write.csv(temp.data_ID, "C:/Users/dakot/Downloads/tempdata_ID.csv")
temp.data_ID <- read.csv("C:/Users/dakot/Downloads/tempdata_ID.csv")

#### joining metadata and climate data ####
metafinal <- meta %>%
  inner_join(temp.data_ID, by=c(id))
metafinal <- left_join(meta, temp.data_ID, join_by(id == id))

write.csv(metafinal,"C:/Users/dakot/Downloads/metafinal.csv")

###### I FINALLY DID IT!!!!!!!!!! ######

meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal.csv", row.names = "X")


