### NOTE: All of this was done locally as the libraries listed weren't available here
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

####### doing all 19 bioclimatic variables ######

TempM <- clim$wc2.1_10m_bio_1
MDiurnalR <- clim$wc2.1_10m_bio_2
Isotherm <- clim$wc2.1_10m_bio_3
TempSeason <- clim$wc2.1_10m_bio_4
MxTempWM <- clim$wc2.1_10m_bio_5
MnTempCM <- clim$wc2.1_10m_bio_6
TempR <- clim$wc2.1_10m_bio_7
MTempWtQ <- clim$wc2.1_10m_bio_8
MTempDQ <- clim$wc2.1_10m_bio_9
MTempWmQ <- clim$wc2.1_10m_bio_10
MTempCQ <- clim$wc2.1_10m_bio_11
Prec <- clim$wc2.1_10m_bio_12
PrecWM <- clim$wc2.1_10m_bio_13
PrecDM <- clim$wc2.1_10m_bio_14
PrecSeason <- clim$wc2.1_10m_bio_15
PrecWtQ <- clim$wc2.1_10m_bio_16
PrecDQ <- clim$wc2.1_10m_bio_17
PrecWmQ <- clim$wc2.1_10m_bio_18
PrecCQ <- clim$wc2.1_10m_bio_19

id <- meta$id
lats <- meta$latitude
lons <- meta$longitude

samples <- data.frame(lons, lats, row.names=id)
temp.data <- samples 
temp.data$TempM <- terra::extract(TempM, samples)
temp.data$MDiurnalR <- terra::extract(MDiurnalR, samples)
temp.data$Isotherm <- terra::extract(Isotherm, samples)
temp.data$TempSeason <- terra::extract(TempSeason, samples)
temp.data$MxTempWM <- terra::extract(MxTempWM, samples)
temp.data$MnTempCM <- terra::extract(MnTempCM, samples)
temp.data$TempR <- terra::extract(TempR, samples)
temp.data$MTempWtQ <- terra::extract(MTempWtQ, samples)
temp.data$MTempDQ <- terra::extract(MTempDQ, samples)
temp.data$MTempWmQ <- terra::extract(MTempWmQ, samples)
temp.data$MTempCQ <- terra::extract(MTempCQ, samples)
temp.data$Prec <- terra::extract(Prec, samples)
temp.data$PrecWM <- terra::extract(PrecWM, samples)
temp.data$PrecDM <- terra::extract(PrecDM, samples)
temp.data$PrecSeason <- terra::extract(PrecSeason, samples)
temp.data$PrecWtQ <- terra::extract(PrecWtQ, samples)
temp.data$PrecDQ <- terra::extract(PrecDQ, samples)
temp.data$PrecWmQ <- terra::extract(PrecWmQ, samples)
temp.data$PrecCQ <- terra::extract(PrecCQ, samples)

write.csv(temp.data, "C:/Users/dakot/Downloads/metfinalall.csv") ## here I needed to clean col names and remove each variable's unique identifier
metfinalall <- read.csv("C:/Users/dakot/Downloads/metfinalall.csv")
colnames(metfinalall)

temp.data_ID19 <- cbind.data.frame(samples, metfinalall)
myVars = c("TempM","MDiurnalR","Isotherm","TempSeason","MxTempWM","MnTempCM","TempR","MTempWtQ","MTempDQ",            
           "MTempWmQ","MTempCQ","Prec","PrecWM","PrecDM","PrecSeason","PrecWtQ","PrecDQ","PrecWmQ","PrecCQ")
temp.data_ID19 <- temp.data_ID19[,c(myVars)]
write.csv(temp.data_ID19, "C:/Users/dakot/Downloads/tempdata_ID19.csv") ## here in excel I manually changed id col name to "id"
temp.data_ID19 <- read.csv("C:/Users/dakot/Downloads/tempdata_ID19.csv")

#### joining metadata and climate data ####
metafinal19 <- left_join(meta, temp.data_ID19, join_by(id == id)) ## here in excel I double checked lat, long and remove lat, long, and x (: 

write.csv(metafinal19,"C:/Users/dakot/Downloads/metafinal19.csv")
