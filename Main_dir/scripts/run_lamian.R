run_lamian <- function(sce, data, permuiter = 5){
  cnt <- counts(sce)
  # select one lineage
  pt <- sce$slingPseudotime_1
  names(pt) <- colnames(cnt) # add name of cells
  pt <- pt[!is.na(pt)]
  selectedCell <- names(pt)
  sample <- data@meta.data %>%
    filter(rownames(.) %in% selectedCell) %>%
    pull(sample)
  sample <- data@meta.data %>%
          filter(rownames(.) %in% selectedCell) %>%
            dplyr::select(sample)
  selectedSample <- data@meta.data %>%
          filter(rownames(.) %in% selectedCell) %>%
          count(sample) %>%
            filter(n > 50) %>%
              pull(sample)
  selectedCell2 <-data@meta.data %>%
          filter(rownames(.) %in% selectedCell) %>%
          filter(sample %in% selectedSample) %>%
              rownames()
  # subset sce
  sce <- sce[,selectedCell2]
  data_sub <- subset(data, cells = selectedCell2)
  # create cell - sample table
  cellanno <- data_sub@meta.data %>%
      dplyr::select(ct, sample) %>%
      mutate(cell = as.character(rownames(.)),
             sample = as.character(sample)) %>%
      dplyr::select(cell, sample)
  # create design matrix
  sex_design <- data_sub@meta.data %>%
          distinct(Predicted_sex, sample) %>%
          mutate(intercept = 1) %>%
          mutate(sex_bin = ifelse(Predicted_sex == "M", 1,
                                  ifelse(Predicted_sex == "F", 0, NA)))
  rownames(sex_design) <- sex_design$sample
  sex_design %>% select(intercept, sex_bin) -> sex_design
  # generate expression matrix
  # only variable features
  data_sub <- FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = 2000)
  var_genes <- VariableFeatures(data_sub)
  GetAssayData(data_sub)[var_genes,]-> expr
  filtr_expr <- expr
  # filtr_expr <- expr[rowSums(expr > 10) > ncol(expr)/100, ]
  filtr_expr <- LogNormalize(filtr_expr)

  # Run lamian
  Res <- lamian_test(
    expr = as.matrix(filtr_expr),
    cellanno = cellanno,
    pseudotime = pt[selectedCell2],
    design = sex_design,
    verbose.output = TRUE,
    test.type = 'variable',
    testvar = 2,
    permuiter = permuiter,
    ncores = 10
    )
  return(Res)
}