# Coding & Data Notes for Transcriptomics Module

## Author: Dakota Kaupu

------------------------------------------------------------------------

**8 Oct 2024 -- Beginning of transcriptomics! (DESeq2)**

Downloading DESeq2 with NO project fix:

```         
library(BiocManager)
BiocManager::install("Biobase", dependencies=T, force=T)
library(DESeq2)
```

Quick recap on simple R functions:

```         
str() = shows x obsv. for y variables
head() = show first 10 samples
tail () = show last 10 samples
dim() = dimensions of file
```

See code for more notes :)

------------------------------------------------------------------------

**10 Oct 2024 -- Exploring in DESeq2**

-   reading in conditions and counts matrix
-   explored counts matrix; looked at reads per sample and reads per gene
-   filtered our data for transcripts with 10 or more reads in 15 or more samples (15 samples ensures we aren't filtering out just samples from one of our dev temps and keeping the other)
-   PCA analysis, accounting for dev and final temp

*Commenting out large chunks*: highlight code, cmd+shift+c

See code for more notes [:

------------------------------------------------------------------------

**15 Oct 2024 -- Visualizing DGE**

-   Created counts plot for dev22v18
-   Created MA plot for dev22v18
-   Created volcano plot, heat map ""

See code for more notes \<:

------------------------------------------------------------------------

**17 Oct 2024 -- Visualizing DGE cont.**

-   Compared devtemps of 18 and 22 between final temps of baseline, A28, and A33
-   Number of DGE genes decreased as final temp increased (at higher acute temp exposures, lower magnitude of DGE, suggesting converged lethal response instead of adapted tolerance)
-   Found shared DGEs between each groups (usually represented by venn diagram)

See code for more notes c:

------------------------------------------------------------------------

**22 Oct 2024 -- Euler plot and scatter plots**

-   Calculated all values for euler plot
-   Created Euler plot (venn diagram thing, except scaled by values)
-   Created scatter plots of copepod unique developmental temp's response to finaltemp (base v 28) (base v 33)

See code for more notes :\>

------------------------------------------------------------------------

**24 Oct 2024 -- Mean Connectivity & Scale Topology**

Review of plots we've made:

-   Bar plot for \# of reads per sample; shows success and variation in sequencing effort
-   PCA; visualizes variation among groups
-   MA plot shows average counts vs logFoldChange; difference in expression between 2 groups relating to relative expression
-   Counts/Point plot of DGE; normalize counts of specific genes \@ different vars (devTemp, regardless of FinalTemp)
-   Volcano plot shows LFC vs -logpval; shows significant up and downregulation
-   Euler plot; visualizes which DGE genes are shared/unique
-   Heatmap; see by color & saturation differences among genes or other matrices of data (correlation, etc.)
-   Scatter plots;  colored by contrast and significance

What we did today:

-   Reviewed scatter plots & combined into one picture
-   Filtered data to BASE data & good genes and samples
-   Made dendrogram AND PCA to look for sample outliers
-   Filtered further (check code) and stabilized
-   Graphed mean connectivity and scale topology (explanation in code)


See code for more notes C:

------------------------------------------------------------------------

**29 Oct 2024 -- Blockwise analyses & Heatmaps**

- By looking at mean connectivity and scale topology, we picked a power of 26
- ran blockwise analysis, but ran out of time

* dont be greedy, but ask for what you need*
See code for more notes ^-^

------------------------------------------------------------------------
