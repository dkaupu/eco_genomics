library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/") ## getwd in console to double-check

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz") 
  ## vcf in console for overview
  ## CHROMs = chromosomes; 8 with 13 scaffolding
  ## variants = # of SNPS
head(vcf)
  ## fixed section lists each SNP
  ## genotype section matrix of SNP rows and genotype columns
      ## last section shows how much reads were indicative of heterozygote (6,8)

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