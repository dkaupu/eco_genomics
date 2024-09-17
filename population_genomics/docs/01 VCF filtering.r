library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/") ## getwd in console to double-check
list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz") 
  ## vcf in console for overview
  ## CHROMs = chromosomes; 8 with 13 scaffolding
  ## variants = # of SNPS
head(vcf)
  ## fixed section lists each SNP
  ## genotype section matrix of SNP rows and genotype columns
      ## GT= genotype, PL= , DP=depth, AD=allele depth (ref # of reads,alt #)
      ## 0/1 = heterozygote, ./.=NA

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")
  ## gff file = where gene information (function, etc.) is stored

## create object that associates vcf file with reference genome
chr1 <- create.chromR(name = "Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chr1)
  ## mapping quality good!
  ## depth evenly distributed! (multiple humps suggest read mapping to wrong spot)

pdf(file = "~/projects/eco_genomics/population_genomics/figures/ChromoPlot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
  ## xlim "from what bp to what bp to plot"
  ## each dot a SNP, red bar for every gene located in genome
dev.off()
  ## closes pdf writer


############################ Starting 9/17/2024 ########################

DP <- extract.gt(vcf, element="DP", as.numeric =T)
dim(DP)
  ## 18233 rows/SNPS, 629 columns/individuals
DP[1:5,1:10]

quantile(DP)

DP[DP==0] <- NA     ## Set values of 0 as true NA values
quantile(DP, na.rm=T) ## 50% of our samples had 11 reads, etc.

## Visualizing matrix of DP and missingness

heatmap.bp(DP[1:1000,], rlabels=F, clabels=F)
  ##individuals=col, loci=rows, white=missing data
  ## bars=mean DP value

## FILTERING BY DEPTH

library(SNPfiltR)

vcf.filt <- hard_filter(vcf, depth=3)
    ## if you dont have more than 3 reads, set to N/A
    ## now 35% missing data

    max_depth(vcf.filt)
        ## huge tail at end may suggest paralogy, still should set max (2*avg)

vcf.filt <- max_depth(vcf.filt, maxdepth = 60) #filter out genotypes with >60 reads/SNP

## FILTERING BY MISSINGNESS

meta <- read.csv("metadata/meta4vcf.csv", header=T) #View(meta) in console to see table

meta2 <- meta[,c(1,4)]  # give me all rows of meta, and columns 1 & 4

names(meta2) <- c("id","pop") ## changing region to "pop" for compatibility
meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=0.75) # remove individuals with "x" or more missing data

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
  ## should assume only 2 alleles for a SNP, if more somethings wrong
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)
  ## number of times you see alternate allele to keep it

vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss,
                                      cutoff=0.5)

DP2 <- extract.gt(vcf.filt.indSNPMiss,
                  element = "DP",
                  as.numeric=T)

heatmap.bp(DP2[5001:10000,], rlabels=F, clabels=F) #plotting loci 5001-10K

write.vcf(vcf.filt.indSNPMiss, 
          "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

           