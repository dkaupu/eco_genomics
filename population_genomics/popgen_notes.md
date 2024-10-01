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
spack load samtools
samtools tview WV_9.sorted.bam ##? for key for samtools
cd ..
cd /variants
zcat Centaurea_filtered.vcf.gz | head -n 20
zcat Centaurea_filtered.vcf.gz | head -n 2000 | vim - ##cmd C to quit
```

**Created "01 VCF filtering.r"** - Created plot with vcfR - Created ChromoPlot.pdf

------------------------------------------------------------------------

**17 Sep 2024 - Filtering VCF files**

-   Filter by depth (3-60 reads)
-   Filter by sample (0.75) and SNP missingness (0.5)
-   Filter out low-frequency alleles (only biallelic, min.mac=1)

See code for more notes

------------------------------------------------------------------------

**19 Sep 2024 - Estimating Genetic Diversity**

Created "02 Diversity Differentiation.r" in popgen docs

-   Loaded in vcfR, tidyverse, & qqman
-   Pulled final filtered vcf file - Pulled meta info from class
-   Created "meta2" with matching ids of vcf (to create same \# of variables)
-   Created "vcf.div" grouping by region
-   Created "chr.main" to pull top 8 unique chromosomes and "chrnum" to number top 8 1-8
-   Created "vcf.div.MHplot" to include chromosome #, added variable stating chromosome & position number
-   Created manhattan plot of 8 chromosomes against Fst among regions

See code for more notes

```         
## plotting issue fix:
X11.options(type="cairo")
options(bitmapType = "cairo")
```

See code for more notes

------------------------------------------------------------------------

**24 Sep 2024 - Hexp and PCA graphs**

Manhattan plot: F~ST~ shows percent average differentiation between population; plot may also show which SNPs have higher differentiation within the chromosome than others. Most SNPS in the graph have low F~ST~ showing neutral selection; outliers may show genes under selection

-   Created independent vcf spreadsheet w/ filtered and concatenated data
-   Created variable to make all H~S~ data into one column and directly graphed it
-   Exported histogram graph showing "Genome-wide expected Heterozygosity"
-   Analyzed non-zero values of H~S~ to reflect all polymorphic *individual* data (removed any SNPs that did NOT experience any heterozygosity)

PCA graph:

-   Filtered data to remove SNPs 500 bp or closer (var "vcf.thin")
-   Combined metadata "vcf_final.filtered.thinned.vcf.gz"
-   Unzipped gz file in system, created simple PCA graph

See code for more notes

------------------------------------------------------------------------

**26 Sep 2024 - PCAs!**

-   Created scree plot to analyze eigen values of PCAs
-   Created PCA plot using PC1&2, grouped by region and continent (PNW and CEU clustering!, overlap = possible colonization populations; PNW suggest unsampled ancestral *Centaurea* sp. or diverged from rest of sp.)
-   Created PCA plot using PC2&3, not a lot of clustering by PC3

```         
## Find PC % or percent variation explained by PCx

CentPCA%eigenvalues [1:5] # show me eigen values 1-5
sum(CentPCA$eigenvalues)
CentPCA$eigenvalues[1]/sum(CentPCA$eigenvalues)
```

See code for more notes

------------------------------------------------------------------------

**01 Oct 2024 - Admixture Analysis & Selection Outliers**

Evaluating cross-entropy and eigen values to consider best K value for admixture analysis

-   Created admixture analysis for K=5 grouped by regions, see some color clustering for PNW and CEU (Can be found in figures folder)
-   Created easily producable admixture analysis for any K groups; only need to change var myK
*i renamed the file before saving; it rolled back to 9/26 file :(, will use code uploaded to brightspace as placeholder*
 
Selection Outlier Analyses
-   Created Manhattan Plot of loci after filtering out rare alleles
-   Created PCAdapt Graph to scan for selection and filtered out highly selected BP's

See code for more notes
