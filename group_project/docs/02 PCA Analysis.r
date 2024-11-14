library(vcfR)
library(tidyverse)
library(qqman)
X11.options(type="cairo")
options(bitmapType = "cairo")

capitula <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/PNW_EU_NE_capitulummeasurements.csv")
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")
vcf.thin <- distance_thin(vcf, min.distance=500) 
meta <- read.csv("~/projects/eco_genomics/group_project/outputs/metafinal.csv", row.names = "X")

geno <- vcf2geno(input.file="/gpfs1/home/d/k/dkaupu/vcf_final.filtered.thinned.vcf",
                 output.file="outputs/vcf_final.filtered.thinned.geno")


CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE) 
CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

######## PCA graphin' #######
ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +  ## close AFTER aes, layer + after
  geom_point(alpha=1) + ## alpha= transparency, could try change to see clusters better
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent") # +
#      xlim(-10,10) + ylim(-10,10)

ggsave("figures/CentPCA_PC1vPC2.pdf", width=6, height=6, units="in")


