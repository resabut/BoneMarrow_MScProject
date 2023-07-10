

# load libraries
library(cellXY)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dittoSeq)

SexAssign <- function(data, genome = "Hs", sex_col = "sex", sample_col = "sample",
                      do.report = FALSE, label_plot = TRUE, do.heatmap = TRUE,
                      min.percent = 0.7, min.ratio = 2, assay = "RNA"){
  sex_sample_data <- NULL
  freq_plot <- NULL
  x_y_heatmap <- NULL
  print("Checking file format...")
  # Check file format
  if (class(data) == "Seurat"){
      orig.format <- "Seurat"
      data.seu <- data
      data.sce <- as.SingleCellExperiment(data, assay = assay)
      rm(data)
    } else if (class(data) == "SingleCellExperiment"){
      orig.format <- "SingleCellExperiment"
      data.sce <- data
      data.seu <- as.Seurat(data)
      DefaultAssay(data.seu) <- assay
      rm(data)
    } else{
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
  print(paste("Successfully read input data. Input:", orig.format))
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
    summarise(n=n()) %>%
    group_by(sample) %>%
    mutate(prop = n/sum(n)) %>%
    pivot_wider(names_from = Pred_cell_sex, values_from = c(n, prop)) %>%
    # fill NA values in proportions and counts
    mutate(n_Female = ifelse(is.na(n_Female), 0, n_Female),
           n_Male = ifelse(is.na(n_Male), 0, n_Male),
           prop_Female = ifelse(is.na(prop_Female), 0, prop_Female),
           prop_Male = ifelse(is.na(prop_Male), 0, prop_Male)) %>%
    mutate(sample_sex = ifelse(prop_Female >= min.percent | prop_Female/prop_Male >= min.ratio , "F",
                            ifelse(prop_Male >= min.percent | prop_Male/prop_Female >= min.ratio, "M",
                                  "inconclusive"))) %>%
    mutate(match = ifelse(is.na(sample_sex), "Inconclusive", ifelse(is.na(sex), "New_annotation",
                          ifelse(sex == sample_sex,"Match", "Mismatch")))) -> sex_assigned_df
  # add metadata to Seurat object
  metadata %>%
    left_join(sex_assigned_df, by = c("sample" = "sample")) -> sex_assigned_metadata
  # add rownames lost during join
  rownames(sex_assigned_metadata) <- rownames(metadata)
  data.seu <- AddMetaData(data.seu, metadata = sex_assigned_metadata$sample_sex, col.name = "Predicted_sex")
  # generate report
  if (do.report){
    print("Generating report...")
    # table
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
      labs(x = "Sample", y = "Proportion of cells", fill = "Cell sex") +
      theme_bw()
    if (label_plot){
        freq_plot <- freq_plot + geom_text(aes(label = round(Count, 2)), position = position_fill(vjust = 0.5))
        }
    # print some summary stats
    match_sum <- sex_sample_data %>% dplyr::ungroup() %>% dplyr::count(match)
    print(match_sum)
    if(match_sum[which(match_sum == "Mismatch"), 2] > 0){ # if there are mismatches
      print(paste("Out of", nrow(sex_sample_data), "samples in total,", match_sum[which(match_sum == "Mismatch"), 2],
                "sample(s) sex prediction did not match the provided annotation!!!!"))
    }
    # make a heatmap
    if(do.heatmap){
    # genes located in the X chromosome that have been reported to escape
    # X-inactivation
    # http://bioinf.wehi.edu.au/software/GenderGenes/index.html
    Xgenes <- c("ARHGAP4", "STS", "ARSD", "ARSL", "AVPR2", "BRS3", "S100G",
                "CHM", "CLCN4", "DDX3X", "EIF1AX", "EIF2S3", "GPM6B",
                "GRPR", "HCFC1", "L1CAM", "MAOA", "MYCLP1", "NAP1L3",
                "GPR143", "CDK16", "PLXNB3", "PRKX", "RBBP7", "RENBP",
                "RPS4X", "TRAPPC2", "SH3BGRL", "TBL1X", "UBA1", "KDM6A",
                "XG", "XIST", "ZFX", "PUDP", "PNPLA4", "USP9X", "KDM5C",
                "SMC1A", "NAA10", "OFD1", "IKBKG", "PIR", "INE2", "INE1",
                "AP1S2", "GYG2", "MED14", "RAB9A", "ITM2A", "MORF4L2",
                "CA5B", "SRPX2", "GEMIN8", "CTPS2", "CLTRN", "NLGN4X",
                "DUSP21", "ALG13", "SYAP1", "SYTL4", "FUNDC1", "GAB3",
                "RIBC1", "FAM9C", "CA5BP1")

    # genes belonging to the male-specific region of chromosome Y (unique genes)
    # http://bioinf.wehi.edu.au/software/GenderGenes/index.html
    Ygenes <- c("AMELY", "DAZ1", "PRKY", "RBMY1A1", "RBMY1HP", "RPS4Y1", "SRY",
                "TSPY1", "UTY", "ZFY", "KDM5D", "USP9Y", "DDX3Y", "PRY", "XKRY",
                "BPY2", "VCY", "CDY1", "EIF1AY", "TMSB4Y", "CDY2A", "NLGN4Y",
                "PCDH11Y", "HSFY1", "TGIF2LY", "TBL1Y", "RPS4Y2", "HSFY2",
                "CDY2B", "TXLNGY", "CDY1B", "DAZ3", "DAZ2", "DAZ4")

    # plot average expression of y and x genes per sample
    # taken from https://github.com/ConsiglioLab/Sex_differences_in_PBMCs/blob/main/Obtain_data
    sample.averages <- AverageExpression(data.seu, features = c(Xgenes, Ygenes), return.seurat = TRUE,
                                         group.by = sample_col)
    Idents(sample.averages) <- names(Idents(sample.averages))
    # scale data for raw counts (contains all the genes)
    # sample.averages <- ScaleData(sample.averages, assay = assay, verbose = FALSE)
    # remove features with 0 expression in all samples
    mean_gene_expr <- rowSums(sample.averages[[assay]]@scale.data)
    var_genes <- names(mean_gene_expr[mean_gene_expr != 0])
    # get M - F sample order for heatmap
    sex_sample_data %>%
      arrange(sample_sex) %>%
      pull(sample) -> sample_order
    # order
    Idents(sample.averages) <- factor(sample.averages@active.ident, levels = sample_order)
    # plot heatmap
    DoHeatmap(sample.averages, features = var_genes, size = 2, draw.lines = FALSE,
              assay = assay) + NoLegend() -> x_y_heatmap
    }


  }
  print(paste("Done! Returning object in", orig.format, "format..."))
  print("Predicted sex is under 'Predicted_sex' metadata column")

  # return object in original format
  if (orig.format == "Seurat"){
    return(list(data.seu, sex_sample_data, freq_plot, x_y_heatmap))
  }
  else if (orig.format == "SingleCellExperiment"){
    # convert to SingleCellExperiment
    return(list(as.SingleCellExperiment(data.seu), sex_sample_data, freq_plot, x_y_heatmap))
  }
}