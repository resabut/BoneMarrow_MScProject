# Sex differences in hematopoiesis, a systems immunology approach
This repository contains the necessary files to 

# BoneMarrow_MScProject

This is a private repository containing the Seurat single-cell RNAseq analysis of bone marrow data.

![graphical abstract](img/grapihcal_abstract.png "Graphical abstract")

## Project structure

The project has the following files and directories:

- **Main_dir/**: This directory contains the main scripts and notebooks for the analysis.
- **README.md**: This file provides an overview of the project and its contents.
- **conda_env.md**: This file lists the conda environment and packages used for the analysis for manual installation.
- **r-integration.txt**: This file contains the instruction for a automatic installation of all the packages in the conda environment.

## How to use

To use this project, follow these steps:

1. Clone or download this repository to your local machine.
2. Create and activate the conda environment (See installation section).
3. Go to the **Main_dir/** directory and run the scripts or notebooks as needed (See run section).


## Installation

Run the following commands to install the conda environment:

```bash
conda env create --file r-integration.txt
conda activate r-integration
```
Next install the R packages using the following command:

```r
remotes::install_github("JinmiaoChenLab/FastIntegrate")
BiocManager::install("BiocParallel")
remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE, dependencies = "Suggest")
remotes::install_github("phipsonlab/cellxy")
remotes::install_github("saeyslab/muscatWrapper")
```

## Run
Go to the **Main_dir/** directory and into every subdirectory and run the scripts or notebooks as needed.
Here's a brief overview of what they do:

