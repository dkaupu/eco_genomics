##################### PCA ###################

library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

X11.options(type="cairo")
options(bitmapType = "cairo")

setwd("~/projects/eco_genomics/population_genomics/")

vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

## Thinning out the SNPs for Linkage Disequilibrium (LD) before running PCA 
## & Admixture to satisfy assumption of independence among loci

vcf.thin <- distance_thin(vcf, min.distance=500) ## eleminate any SNP if its 500bp closer to next
      ## 15454 SNPS => 3646 SNPS, a lot of them were close together due to methodology of GBS

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]), ]
dim(meta2)

write.vcf(vcf.thin,"outputs/vcf_final.filtered.thinned.vcf.gz")

## hide the uncompressed vcf file too big for github outside of repo

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")  

## system similar to running it in shell
## uncompress x found in "" and move to > titled y

geno <- vcf2geno(input.file="/gpfs1/home/d/k/dkaupu/vcf_final.filtered.thinned.vcf",
                 output.file="outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE) ## explicitly LEA library pca function

plot(CentPCA$projections,
     col=as.factor(meta2$region), ## color by region
     legend("bottomright", legend=as.factor(unique(meta2$region)), 
                           fill=as.factor(unique(meta2$region))))
