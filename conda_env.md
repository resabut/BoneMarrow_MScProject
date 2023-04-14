# Environment set up
Creates a conda environment and installs R 4.2.0, and all the required packages.

## Conda environment set up
```bash
conda create -n MScProject39
conda activate MScProject39
conda install -c r r=4.2 r-essentials
R # starts the R console
```

## Install packages
Copy the code into the R console and run it.
```r
install.packages("knitr")
install.packages("rmarkdown")
install.packages("remotes")
install.packages("BiocManager")
remotes::install_version("Seurat", version = "4.3.0")
install.packages("HGNChelper")
install.packages("openxlsx")
install.packages("clustree")
remotes::install_github("stemangiola/CuratedAtlasQueryR")
BiocManager::install("GEOquery")
BiocManager::install("dittoSeq")

```


## Load packages
Test if the packages are installed correctly.
```r
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(dittoSeq)
library(dplyr)
library(clustree)
library(GEOquery)
library(CuratedAtlasQueryR)
```

## Session
Verify the R version and installed packages.
```r
sessionInfo()
```



