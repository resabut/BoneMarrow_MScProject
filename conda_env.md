# Environment set up
Creates a conda environment and installs R 4.2.0, and all the required packages.

## Conda environment set up
```bash
conda create -n MScProject39
conda activate MScProject39
conda install -c r r=4.2 r-essentials
conda install -c conda-forge r-monocle3
conda install r-duckdb
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
BiocManager::install("GEOquery")
BiocManager::install("dittoSeq")
remotes::install_github('satijalab/seurat-wrappers')
BiocManager::install("patchwork")
BiocManager::install(c("Biobase", "SingleCellExperiment", "batchelor", "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "S4Vectors", "SummarizedExperiment", "pcaMethods"))

remotes::install_github('cole-trapnell-lab/leidenbase')
remotes::install_github('cole-trapnell-lab/monocle3')
# this package is gone from CRAN in 4.2.3, it is needed to install monocle3
remotes::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')

# For curatedatlastqueryr

remotes::install_github("stemangiola/CuratedAtlasQueryR")

# for cell xy
 mamba install -c bioconda bioconductor-annotationdbi
remotes::install_github("phipsonlab/cellxy")

# for speckle
mamba install r-locfit bioconductor-edger
remotes::install_github("phipsonlab/speckle")

# FOR FAST INTEGRATION
mamba install r-rcppgsl r-rcppziggurat r-rfast
#### or maybe just 
mamba install r-rcpp
remotes::install_github("git@github.com:JinmiaoChenLab/FastIntegrate.git")
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
library(patchwork)
library(monocle3)

```

## Session
Verify the R version and installed packages.
```r
sessionInfo()
```



# for monocle3
This doesnt work for me. 
```bash
conda create -n r-monocle
conda activate r-monocle
conda install -c r r=4.2 r-essentials
```
## Install packages
Copy the code into the R console and run it.
```r
if (!any(rownames(installed.packages()) == "remotes")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("remotes")
}
library(remotes)

if (!any(rownames(installed.packages()) == "Seurat")){
  BiocManager::install("Seurat")
}
library(Seurat)

if (!any(rownames(installed.packages()) == "R.utils")){
  BiocManager::install("R.utils")
}
library(R.utils)

remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

if (!any(rownames(installed.packages()) == "patchwork")){
  BiocManager::install("patchwork")
}
library(patchwork)

## Monocle3 dependancies
BiocManager::install(c("Biobase", "SingleCellExperiment", "batchelor", "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "S4Vectors", "SummarizedExperiment", "pcaMethods"))

remotes::install_github('cole-trapnell-lab/leidenbase')
remotes::install_github('cole-trapnell-lab/monocle3')

## Test out the installation
library(monocle3)

## Mac users may also experience installation problems due to Xcode or gfortran.
## Instructions: https://cole-trapnell-lab.github.io/monocle3/docs/installation/

# If you had issues with leidenbase and gfortran, trying downloading and installing a newer gfortran
# https://gcc.gnu.org/wiki/GFortranBinaries

## Possible you may have to install Rtools if the above fails
## http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/#1

sessionInfo()
```



# for h5ad files
## Conda
```bash
 mamba install r-reticulate 
```
## R

```r
library(reticulate)
```

# for loom files
## conda
```bash
mamba install r-hdf5r
```
## R
```r
library(hdf5r)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
```


# r-integration
## conda
```bash
mamba create -n r-integration
mamba activate r-integration
mamba install r-base r-essentials
mamba install r-seurat r-seuratobject r-data.table r-matrix r-tictoc r-dplyr r-pbmcapply r-stringr
mamba install r-remotes
# install compilers
mamba install gcc
# biocmanager
mamba install r-biocmanager
mamba install bioconductor-genomeinfodbdata
mamba install bioconductor-scdblfinder
mamba install bioconductor-dittoseq
mamba install bioconductor-biocstyle bioconductor-cellbench bioconductor-scater
mamba install bioconductor-annotationdbi
mamba install bioconductor-slingshot
mamba install r-clustree
mamba install bioconductor-tradeseq
```

## R
```r
remotes::install_github("JinmiaoChenLab/FastIntegrate")

BiocManager::install("BiocParallel")

remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE, 
dependencies = "Suggest")
remotes::install_github("phipsonlab/cellxy")


```