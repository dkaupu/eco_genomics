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

See code for more notes <:

------------------------------------------------------------------------

**17 Oct 2024 -- Visualizing DGE cont.**

-   Compared devtemps of 18 and 22 between final temps of baseline, A28, and A33
-   Number of DGE genes decreased as final temp increased (at higher acute temp exposures, lower magnitude of DGE, suggesting converged lethal response instead of adapted tolerance)
-   Found shared DGEs between each groups (usually represented by venn diagram)

See code for more notes c:

------------------------------------------------------------------------


