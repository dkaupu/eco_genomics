# Coding & Data Notes for Population Genetics Module

## Author: Dakota Kaupu

------------------------------------------------------------------------

**10 Sep 2024 - Intro to Centaurea GBS Data and working with VCF files**

-   Analyzing GBS data from 3 regions (EU,NE, PNW) w/ variant call format files (VCF's)

Using cluster VACC shell access:

```         
/gpfs1/cl/pbio3990/example_data/
ll example_data/
zcat file.gz|head -n X   ##(unzip file with X amount of lines)
    head = first 10 lines
    | = send output of one function to another
```

Reading fast q file:

-   Shows sequence and corresponding Q score
-   Q scores represented by:

```         
Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                  |         |         |         |         |
   Quality score: 0........10........20........30........40 
```

------------------------------------------------------------------------

**12 Sep 2024 - Viewing VCF files and filtering**

Bioinformatics Pipeline:

-   Trimmed adapters and low-quality sequence reads with fastp
-   Aligned reads with reference genome [*C. solstitialis*](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030169165.1/)
-   Alignments converted to binary, sorted, indexed, and genotyped
-   VCF file output filtered using vcftools


Using cluster VACC shell access to look at alignments:

```         
/gpfs1/cl/pbio3990/example_data/
?
samtools tview WV_9.sorted.bam ##? for key for samtools
cd ..
cd /variants
zcat Centaurea_filtered.vcf.gz | head -n 20
zcat Centaurea_filtered.vcf.gz | head -n 2000 | vim - ##cmd C to quit
```

**Created "01 VCF filtering.r"**
-   Created plot with vcfR
-   Created ChromoPlot.pdf 