# Environment set up
Creates a conda environment and installs R 4.2.0, and all the required packages.



## conda
Run this in the bash terminal 
```bash
mamba create -n r-integration
mamba activate r-integration
#install r
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
mamba install r-paletteer
mamba install r-vegan
mamba install bioconductor-topgo
mamba install r-ggridges
mamba install r-gtextras r-svglite
mamba install bioconductor-muscat
mamba install r-ggpubr bioconductor-sva
mamba install bioconductor-enrichplot bioconductor-clusterprofiler bioconductor-pathview
mamba install r-ggtext r-magick
```

## R
Next run this in R to install the appropriate 
```r
remotes::install_github("JinmiaoChenLab/FastIntegrate")
BiocManager::install("BiocParallel")
remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE, dependencies = "Suggest")
remotes::install_github("phipsonlab/cellxy")
remotes::install_github("saeyslab/muscatWrapper")
```

