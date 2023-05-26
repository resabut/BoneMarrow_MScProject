

# load libraries
library(cellXY)
library(SingleCellExperiment)
library(org.Hs.eg.db)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dittoSeq)

SexAssign <- function(data, sex_col = "sex", genome = "Hs", sample_col = "sample"){
  print("Checking file format...")
  # Check file format
  if (class(data) == "Seurat"){
      orig.format <- "Seurat"
      data.seu <- data
      data.sce <- as.SingleCellExperiment(data)
      rm(data)
    }
    else if (class(data) == "SingleCellExperiment"){
      orig.format <- "SingleCellExperiment"
      data.sce <- data
      data.seu <- as.Seurat(data)
      rm(data)
    }
    else{
      stop("Input data must be Seurat or SingleCellExperiment object")
    }
  print(paste("Succesfully read input data. Input:", orig.format))
  # Get counts
  print("Computing counts...")
  data.counts <- counts(data.sce)
  # correct gene names if necessary
  # ann <- select(org.Hs.eg.db, keys=rownames(sc_data),
  #               columns=c("GENENAME", "SYMBOL"), keytype="SYMBOL")
  # m <- match(rownames(data.counts), ann$GENENAME)
  # rownames(data.counts) <- ann$SYMBOL[m]
  # classify sex of cells using cellXY
  print("Identifying sex of cells with cellXY...")
  sex_pred <- classifySex(data.counts, genome = genome)
  print("Adding metadata...")
  # add metadata to Seurat object
  data.seu <- AddMetaData(data.seu, metadata = sex_pred, col.name = "Pred_cell_sex")
  # change column names
  data.seu@meta.data -> metadata
  metadata[["sample"]] <- metadata[[sample_col]]
  metadata[["sex"]] <- metadata[[sex_col]]
  # compute sample sex
  metadata %>%
    group_by(sample, Pred_cell_sex, sex) %>%
    count() %>%
    group_by(sample) %>%
    mutate(prop = n/sum(n)) %>%
    pivot_wider(names_from = Pred_cell_sex, values_from = c(n, prop)) %>%
    mutate(sample_sex = ifelse(prop_Female >= 0.70 | prop_Female/prop_Male >= 2, "F",
                            ifelse(prop_Male >= 0.70 | prop_Male/prop_Female >= 2, "M",
                                  "inconclusive"))) %>%
    mutate(match = ifelse(sex == sample_sex,
                          "match", "mismatch")) -> sex_assigned_df
  # add metadata to Seurat object
  metadata %>%
    left_join(sex_assigned_df, by = c("sample" = "sample")) -> sex_assigned_metadata
  # add rownames lost during join
  rownames(sex_assigned_metadata) <- rownames(metadata)
  data.seu <- AddMetaData(data.seu, metadata = sex_assigned_metadata$sample_sex, col.name = "Predicted_sex")
  print(paste("Done! Returning object in", orig.format, "format..."))
# return object in original format
  if (orig.format == "Seurat"){
    return(data.seu)
  }
  else if (orig.format == "SingleCellExperiment"){
    # convert to SingleCellExperiment
    return(as.SingleCellExperiment(data.seu))
  }
}

