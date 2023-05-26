

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
  # Check if it is Seurat object
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
  # Get counts
  data.counts <- counts(data.sce)
  # correct gene names if necessary
  # ann <- select(org.Hs.eg.db, keys=rownames(sc_data),
  #               columns=c("GENENAME", "SYMBOL"), keytype="SYMBOL")
  # m <- match(rownames(data.counts), ann$GENENAME)
  # rownames(data.counts) <- ann$SYMBOL[m]
  # classify sex of cells using cellXY
  sex_pred <- classifySex(data.counts, genome = genome)
  # add metadata to Seurat object
  data.seu <- AddMetaData(data.seu, metadata = sex_pred, col.name = "Pred_cell_sex")
  # compute sample sex
  data.seu@meta.data %>%
    group_by({{sample_col}}, Pred_cell_sex, sex) %>%
    count() %>%
    group_by({{sample_col}}) %>%
    mutate(prop = n/sum(n)) %>%
    pivot_wider(names_from = Pred_cell_sex, values_from = c(n, freq)) %>%
    mutate(sample_sex = ifelse(freq_Female >= 0.70 | freq_Female/freq_Male >= 2, "F",
                           ifelse(freq_Male >= 0.70 | freq_Male/freq_Female >= 2, "M",
                                  "inconclusive"))) %>%
    mutate(match = ifelse({{sex_col}} == sample_sex,
                          "match", "mismatch")) -> sex_assigned_df
  # add metadata to Seurat object
  data.seu[[sample_col]] %>%
    left_join(sex_assigned_df, by = structure(names=sample_col, .Data=sample_col)) -> sex_assigned_metadata
  data.seu <- AddMetaData(data.seu, metadata = sex_assi_metadata$sample_sex, col.name = "Predicted_sex")
# return object in original format
  if (orig.format == "Seurat"){
    return(data.seu)
  }
  else if (orig.format == "SingleCellExperiment"){
    # convert to SingleCellExperiment
    return(as.SingleCellExperiment(data.seu))
  }
}

