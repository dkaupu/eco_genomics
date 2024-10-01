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

CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject") ## reload outputs

show(CentPCA) ## shows basic info of object
plot(CentPCA) ## y = eigen value/ effectivity of PCA line

# plot(CentPCA$projections,
#     col=as.factor(meta2$region), ## color by region
#     legend("bottomright", legend=as.factor(unique(meta2$region)), 
#                           fill=as.factor(unique(meta2$region))))

######################### 9/26/2024 #######################

## PC1 v PC2
ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +  ## close AFTER aes, layer + after
       geom_point(alpha=1) + ## alpha= transparency, could try change to see clusters better
       labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent") # +
#      xlim(-10,10) + ylim(-10,10)

ggsave("figures/CentPCA_PC1vPC2.pdf", width=6, height=6, units="in")

## PC2 v PC3
ggplot(as.data.frame(CentPCA$projections),
       aes(x=V2, y=V3, color=meta2$region, shape=meta2$continent)) +  ## close AFTER aes, layer + after
  geom_point(alpha=1) + ## alpha= transparency, could try change to see clusters better
  labs(title="Centaurea genetic PCA", x="PC2", y="PC3", color="Region", shape="Continent")

############################## 10/1/2024 #####################

pdf("figures/Admixture_K5.pdf", width=10, height=5)
barplot(as.matrix(t(myKQmeta[ , 1:myK])), ## t = transpose, graph [ , ] all indv for values 1-myK/5
          border = NA,
          space = 0,
          col=my.colors[1:myK],
          xlab = "Geographic Regions", 
          ylab = "Ancestry Proportions",
          main = paste0("Ancestry matrix K=", myK))
axis (1,
        at=1:length(myKQmeta$region),
        labels=myKQmeta$region,
        tick = F,
        cex.axis = 0.5,
        las = 3) ## horizontal labels
> dev.off()