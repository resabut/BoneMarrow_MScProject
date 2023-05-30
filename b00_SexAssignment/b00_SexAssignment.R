

# load libraries
library(cellXY)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dittoSeq)

SexAssign <- function(data, genome = "Hs", sex_col = "sex", sample_col = "sample", report = FALSE, label_plot = TRUE,
                      min.percent = 0.7, min.ratio = 2, assay = "RNA"){
  sex_sample_data <- NULL
  freq_plot <- NULL
  print("Checking file format...")
  # Check file format
  if (class(data) == "Seurat"){
      orig.format <- "Seurat"
      data.seu <- data
      data.sce <- as.SingleCellExperiment(data, assay = assay)
      rm(data)
    }
    else if (class(data) == "SingleCellExperiment"){
      orig.format <- "SingleCellExperiment"
      data.sce <- data
      data.seu <- as.Seurat(data)
      DefaultAssay(data.seu) <- assay
      rm(data)
    }
    else{
      stop("Input data must be Seurat or SingleCellExperiment object")
    }
  # Check if required columns are present
  if (!(sex_col %in% colnames(data.seu@meta.data) && sample_col %in% colnames(data.seu@meta.data))) {
    stop("Input data must contain columns for sex and sample IDs")
  }
  # Check if sample column contains valid values
  if (!all(!is.na(data.seu[[sample_col]]))) {
    stop("Sample column must not contain missing values")
  }
  print(paste("Succesfully read input data. Input:", orig.format))
  # Get counts
  print("Computing counts...")
  data.counts <- counts(data.sce)
  rm(data.sce)
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
    mutate(sample_sex = ifelse(prop_Female >= min.percent | prop_Female/prop_Male >= min.ratio, "F",
                            ifelse(prop_Male >= min.percent | prop_Male/prop_Female >= min.ratio, "M",
                                  "inconclusive"))) %>%
    mutate(match = ifelse(sex == sample_sex,
                          "match", "mismatch")) -> sex_assigned_df
  # add metadata to Seurat object
  metadata %>%
    left_join(sex_assigned_df, by = c("sample" = "sample")) -> sex_assigned_metadata
  # add rownames lost during join
  rownames(sex_assigned_metadata) <- rownames(metadata)
  data.seu <- AddMetaData(data.seu, metadata = sex_assigned_metadata$sample_sex, col.name = "Predicted_sex")
  # generate report
  if (report){
    print("Generating report...")
    sex_sample_data <- sex_assigned_df
    # make frequency plot
    bar_data <- as.data.frame(table(unlist(data.seu[["Pred_cell_sex"]]),
                                    unlist(data.seu[[sample_col]])))
    colnames(bar_data) <- c("Cell_sex", "Sample", "Count")
    freq_plot <- ggplot(bar_data, aes(x = Sample, y = Count, fill = Cell_sex)) +
      geom_col(position = "fill") +
      # set colors
      scale_fill_manual(values = c("purple", "yellow", "gray")) +
      # add table with sex assignment
      labs(x = "Sample", y = "Count", fill = "Cell sex") +
      theme_bw()
    if (label_plot){
        freq_plot <- freq_plot + geom_text(aes(label = round(Count, 2)), position = position_fill(vjust = 0.5))
        }

  }
    print(paste("Done! Returning object in", orig.format, "format..."))

  # return object in original format
  if (orig.format == "Seurat"){
    return(list(data.seu, sex_sample_data, freq_plot))
  }
  else if (orig.format == "SingleCellExperiment"){
    # convert to SingleCellExperiment
    return(list(as.SingleCellExperiment(data.seu), sex_sample_data, freq_plot))
  }
}

