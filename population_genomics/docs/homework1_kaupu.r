library(vcfR)
library(SNPfiltR)
library(tidyverse)
library(LEA)

options(bitmapType = "cairo")

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz") 
head(vcf)

############################ filtering by depth ########################

vcf.filt <- hard_filter(vcf, depth=3) ## if you dont have more than 3 reads, set to N/A

max_depth(vcf.filt)   ## huge tail at end may suggest paralogy, still should set max (2*avg)

vcf.filt <- max_depth(vcf.filt, maxdepth = 60) #filter out genotypes with >60 reads/SNP

## This has all been kept the same from previous filterings

############################ filtering by missingness ########################

meta <- read.csv("metadata/meta4vcf.csv", header=T)

meta2 <- meta[,c(1,4)]  # give me all rows of meta, and columns 1 & 4 (ID and region)

names(meta2) <- c("id","pop") ## changing region to "pop" for compatibility
meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

## Originally we filtered for cutting off missingness at 0.75 missing data or greater, I'll also be doing 
## I'll also be doing 1 and 0.50; variables just renamed accordingly

vcf.filt.indMiss100 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=1)

vcf.filt.indMiss100 <- filter_biallelic(vcf.filt.indMiss100) ## should assume only 2 alleles for a SNP, if more somethings wrong
vcf.filt.indMiss100 <- min_mac(vcf.filt.indMiss100, min.mac = 1) ## number of times you see alternate allele to keep it

vcf.filt.indSNPMiss100 <- missing_by_snp(vcf.filt.indMiss100,
                                      cutoff=0.5)

write.vcf(vcf.filt.indSNPMiss100, 
         "~/projects/eco_genomics/population_genomics/outputs/vcf_final100.filtered.vcf.gz")

dim(vcf.filt.indSNPMiss100) # variants = SNP, cols = s# of indv

############################ DIVERSITY ANALYSIS ########################

## calling in indv vcf files (50, 75, 100) & meta data
vcf50 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final50.filtered.vcf.gz")
vcf75 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final75.filtered.vcf.gz")
vcf100 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final100.filtered.vcf.gz")

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

## find all meta ids also present in colnames of vcf files and make a new variable
## luckily we only have three files; I'll just do it manually
meta50 <- meta[meta$id %in% colnames(vcf50@gt[,-1]),] 
meta75 <- meta[meta$id %in% colnames(vcf75@gt[,-1]),] 
meta100 <- meta[meta$id %in% colnames(vcf100@gt[,-1]),] 

## calculate diversity stats using the genetic_diff fxn in vcfR, again x3
vcf.div50 <- genetic_diff(vcf50,
                        pops = as.factor(meta50$region), 
                        method = "nei")
vcf.div75 <- genetic_diff(vcf75,
                        pops = as.factor(meta75$region), 
                        method = "nei")
vcf.div100 <- genetic_diff(vcf100,
                        pops = as.factor(meta100$region), 
                        method = "nei")

str(vcf.div50)

chr.main50 <- unique(vcf.div50$CHROM)[1:8] ## what's unique in CHROM? CM=chromosome, take first 8 entries
chr.main75 <- unique(vcf.div75$CHROM)[1:8]
chr.main100 <- unique(vcf.div100$CHROM)[1:8]

chrnum50 <- as.data.frame(cbind(chr.main50, seq(1,8, 1))) ## number vars in chr.main 1-8, intervals of 1
chrnum75 <- as.data.frame(cbind(chr.main75, seq(1,8, 1)))
chrnum100 <- as.data.frame(cbind(chr.main100, seq(1,8, 1)))

vcf.div.MHplot50 <- left_join(chrnum50, vcf.div50, join_by(chr.main50==CHROM)) ## merging vcf.div w/ chrnum
vcf.div.MHplot75 <- left_join(chrnum75, vcf.div75, join_by(chr.main75==CHROM))
vcf.div.MHplot100 <- left_join(chrnum100, vcf.div100, join_by(chr.main100==CHROM))
head(vcf.div.MHplot) 

vcf.div.MHplot50 <- vcf.div.MHplot50 %>% ## create new variable; chromosome_position number
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main50,"_",POS))
vcf.div.MHplot75 <- vcf.div.MHplot75 %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main75,"_",POS))
vcf.div.MHplot100 <- vcf.div.MHplot100 %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main100,"_",POS))

vcf.div.MHplot50$V2 = as.numeric(vcf.div.MHplot50$V2)   ## originally character data, now numeric
vcf.div.MHplot50$POS = as.numeric(vcf.div.MHplot50$POS)

vcf.div.MHplot75$V2 = as.numeric(vcf.div.MHplot75$V2)  
vcf.div.MHplot75$POS = as.numeric(vcf.div.MHplot75$POS)

vcf.div.MHplot100$V2 = as.numeric(vcf.div.MHplot100$V2)  
vcf.div.MHplot100$POS = as.numeric(vcf.div.MHplot100$POS)

########## calculating mean Hs and SD for every region ##########
########## calculating number of loci with Hs=0 and Hs>0 for every region #########

vcf.div.MHplot50 %>%  
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>%
  filter(value=0) %>% ## filter and keep values that ARE NOT != zero, let's just keep values that ARE polymorphic
  summarise(avg_Hs=mean(value), ## calculate y for values I'm giving you and call it "x"
            StdDev_Hs=sd(value), 
            N_HS=n())           ## n=SNPs 
vcf.div.MHplot50 %>%  
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>%
  filter(value!=0) %>% ## filter and keep values that ARE NOT != zero, let's just keep values that ARE polymorphic
  summarise(avg_Hs=mean(value), ## calculate y for values I'm giving you and call it "x"
            StdDev_Hs=sd(value), 
            N_HS=n())           ## n=SNPs 

vcf.div.MHplot75 %>%  
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>%
  filter(value=0) %>% 
  summarise(avg_Hs=mean(value), 
            StdDev_Hs=sd(value), 
            N_HS=n())          
vcf.div.MHplot75 %>%  
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>%
  filter(value!=0) %>% 
  summarise(avg_Hs=mean(value), 
            StdDev_Hs=sd(value), 
            N_HS=n())           

vcf.div.MHplot100 %>%  
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>%
  filter(value=0) %>% 
  summarise(avg_Hs=mean(value),
            StdDev_Hs=sd(value), 
            N_HS=n())         
vcf.div.MHplot100 %>%  
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>%
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), 
            StdDev_Hs=sd(value), 
            N_HS=n())         

################# PCA graphin' ####################
setwd("~/projects/eco_genomics/population_genomics/")
getwd()
vcf50 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final50.filtered.vcf.gz") 
vcf75 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final75.filtered.vcf.gz") 
vcf100 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final100.filtered.vcf.gz") 

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

## Thinning out the SNPs for Linkage Disequilibrium (LD) before running PCA
vcf.thin50 <- distance_thin(vcf50, min.distance=500) ## eleminate any SNP if its 500bp closer to next
vcf.thin75 <- distance_thin(vcf75, min.distance=500)
vcf.thin100 <- distance_thin(vcf100, min.distance=500)

meta.thin50 <- meta[meta$id %in% colnames(vcf.thin50@gt[,-1]), ]
meta.thin75 <- meta[meta$id %in% colnames(vcf.thin75@gt[,-1]), ]
meta.thin100 <- meta[meta$id %in% colnames(vcf.thin100@gt[,-1]), ]

write.vcf(vcf.thin50,"~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final50.filtered.thinned.vcf.gz")
write.vcf(vcf.thin75,"~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final75.filtered.thinned.vcf.gz")
write.vcf(vcf.thin100,"~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final100.filtered.thinned.vcf.gz")

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final50.filtered.thinned.vcf.gz > ~/vcf_final50.filtered.thinned.vcf")
system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final75.filtered.thinned.vcf.gz > ~/vcf_final75.filtered.thinned.vcf")
system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/homework/vcf_final100.filtered.thinned.vcf.gz > ~/vcf_final100.filtered.thinned.vcf")

geno50 <- vcf2geno(input.file="/gpfs1/home/d/k/dkaupu/vcf_final50.filtered.thinned.vcf", 
                   output.file="outputs/vcf_final50.filtered.thinned.geno") 
geno75 <- vcf2geno(input.file="/gpfs1/home/d/k/dkaupu/vcf_final75.filtered.thinned.vcf", 
                   output.file="outputs/vcf_final75.filtered.thinned.geno")
geno100 <- vcf2geno(input.file="/gpfs1/home/d/k/dkaupu/vcf_final100.filtered.thinned.vcf", 
                    output.file="outputs/vcf_final100.filtered.thinned.geno")

CentPCA50 <- LEA::pca("outputs/vcf_final50.filtered.thinned.geno", scale=TRUE) ## explicitly LEA library pca function
CentPCA75 <- LEA::pca("outputs/vcf_final75.filtered.thinned.geno", scale=TRUE)
CentPCA100 <- LEA::pca("outputs/vcf_final100.filtered.thinned.geno", scale=TRUE)

## CentPCA50 <- load.pcaProject("vcf_final50.filtered.thinned.pcaProject")    ## only if you want to reload after doing everything before
## CentPCA75 <- load.pcaProject("vcf_final75.filtered.thinned.pcaProject")
## CentPCA100 <- load.pcaProject("vcf_final100.filtered.thinned.pcaProject")

## PC1 v PC2 for all three!!!!

ggplot(as.data.frame(CentPCA50$projections),
       aes(x=V1, y=V2, color=meta.thin50$region, shape=meta.thin50$continent)) + 
       geom_point(alpha=0.75) + ## alpha= transparency
       labs(title="50%", x="PC1", y="PC2", color="Region", shape="Continent") # +

ggsave("figures/CentPCA50_PC1vPC2.pdf", width=6, height=6, units="in")

ggplot(as.data.frame(CentPCA75$projections),
       aes(x=V1, y=V2, color=meta.thin75$region, shape=meta.thin75$continent)) + 
       geom_point(alpha=0.75) +
       labs(title="75%", x="PC1", y="PC2", color="Region", shape="Continent") # +

ggsave("figures/CentPCA75_PC1vPC2.pdf", width=6, height=6, units="in")

ggplot(as.data.frame(CentPCA100$projections),
       aes(x=V1, y=V2, color=meta.thin100$region, shape=meta.thin100$continent)) + 
       geom_point(alpha=0.5) +
       labs(title="100%", x="PC1", y="PC2", color="Region", shape="Continent") # +

ggsave("figures/CentPCA100_PC1vPC2.pdf", width=6, height=6, units="in")

