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
