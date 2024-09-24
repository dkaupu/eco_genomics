## Estimating diversity and genetic differentiation in the filtered Centaurea data

library(vcfR)
library(tidyverse)
library(qqman)

## !! plotting issue fix:
X11.options(type="cairo")
options(bitmapType = "cairo")

## read in VCF file 
vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

## read in metadata (pop of origin, region, continent, etc.)
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta)
dim(meta)

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),] 
  ## find all meta ids also present in colnames of vcf files and make a new variable
  ## -1; skip over right column; which is not individual. info

## calculate diversity stats using the genetic_diff fxn in vcfR
vcf.div <- genetic_diff(vcf,
                        pops = as.factor(meta2$region), ##choosing region as grouping factor
                        method = "nei")

str(vcf.div) ## shows structure of variable just made
chr.main <- unique(vcf.div$CHROM)[1:8] 
  ## what are the unique values in CHROM? CM=chromosome, JARY=scaffolds
  ## take first 8 entries

chrnum <- as.data.frame(cbind(chr.main, seq(1,8, 1))) 
  ## number vars in chr.main 1-8, intervals of 8

### merge vcf.div with chr.main; need to assign chrom #'s with integers

vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))
  ## join left frame chrnum with right frame vcf.div; chr.main corresponds to CHROM
head(vcf.div.MHplot) #our V2 is there!

vcf.div.MHplot <- vcf.div.MHplot %>%
                          filter(Gst>0) %>%
                          mutate(SNP=paste0(chr.main,"_",POS))
  ## create new variable chromosome_position number

vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2) ## originally character data, now numeric
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)

manhattan(vcf.div.MHplot,
          chr = "V2",
          bp = "POS",
          p = "Gst",
          col = c("blue4","orange3"),
          logp = F,
          ylab = "Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.5)) ## each dot a SNP, y=FST=allele frequency divergence btw pop.
    ## everything above line in top 1%; most differentiated part of these regions
    ## quantile value 0.5=median

####################### STARTING 9/24/2024 ###################

write.csv(vcf.div.MHplot,"~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_byRegion.csv",
          quote=F,
          row.names=F) ## export vcf.div.MHplot to spreadsheet, don't put name in quotes and give row names

names(vcf.div.MHplot) ## columns 4-9 have Hs values

vcf.div.MHplot %>%  
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%   ## stack columns 4 through 9 into one
  ggplot(aes(x=value, fill=name)) + ## ggplot w/ all Hs values, color them by diff. names
  geom_histogram(position="identity",alpha=0.5, bins=50) + ## alpha = opacity for every histogram
  labs (title= "Genome-wide expected heterozygosity (Hs)", fill="Regions",
        x="Gene diversity (Hs) within Regions",y="Counts of SNPs")
## fill = key, x=x-axis, y=y-axis
## ggplot = "layers" on features by +

ggsave("Histogram_GenomicDiversity_byRegion.pdf", 
       path="~/projects/eco_genomics/population_genomics/figures/")

vcf.div.MHplot %>%  
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>% 
  group_by(name) %>%
  filter(value!=0) %>% ## filter and keep values that ARE NOT != zero, let's just keep values that ARE polymorphic
  summarise(avg_Hs=mean(value), ## calculate y for values I'm giving you and call it "x"
            StdDev_Hs=sd(value), 
            N_HS=n())           ## n=SNPs 
    
    



