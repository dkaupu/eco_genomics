library(vcfR)
library(SNPfiltR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz") 
head(vcf)

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")

chr1 <- create.chromR(name = "Chromosome 1", vcf=vcf, seq=dna, ann=gff) ## associating vcf w/ ref genome

plot(chr1) ## mapping quality good!

############################ filtering by depth ########################

DP <- extract.gt(vcf, element="DP", as.numeric =T)
dim(DP)

quantile(DP)

DP[DP==0] <- NA     ## Set values of 0 as true NA values
quantile(DP, na.rm=T) ## 50% of our samples had 11 reads, etc.

heatmap.bp(DP[1:1000,], rlabels=F, clabels=F)

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

## Originally we filtered for cutting off missingness at 0.75 missing data or greater

missingness=0.75 ## lets just make a variable to make it easier

vcf.filt.indMiss75 <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=missingness)

vcf.filt.indMiss75 <- filter_biallelic(vcf.filt.indMiss) ## should assume only 2 alleles for a SNP, if more somethings wrong
vcf.filt.indMiss75 <- min_mac(vcf.filt.indMiss, min.mac = 1) ## number of times you see alternate allele to keep it

vcf.filt.indSNPMiss75 <- missing_by_snp(vcf.filt.indMiss75,
                                      cutoff=0.5)

DP2 <- extract.gt(vcf.filt.indSNPMiss75,
                  element = "DP",
                  as.numeric=T)

heatmap.bp(DP2[1:1000,], rlabels=F, clabels=F) #loci 5001-10K

# write.vcf(vcf.filt.indSNPMiss, 
#          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

