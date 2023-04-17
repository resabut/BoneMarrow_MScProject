rm(list=ls())
# Load libraries
library(dplyr)
library(tidyr)
library(openxlsx)
# Load data
data <- read.table("data/marker_gene_lists/normal_BMA_markerInfor.txt", header = TRUE, sep = "\t")

# Parse data
# group by 3rd column and list all 4th column values
sc_type_format <- data %>% group_by(Celltype) %>% summarise("MarkerGene" = paste(MarkerGene, collapse = ","))
# add dataset column
sc_type_format <- cbind("Dataset" = data$Dataset, sc_type_format)
# rename columns according to sc-type format
colnames(sc_type_format) <- c("tissueType", "cellName", "geneSymbolmore1")
# write to xlsx
write.xlsx(sc_type_format,
           "data/marker_gene_lists/normal_BMA_markerInfor_sc-type.xlsx",
           sheetName = "Sheet1", row.names = FALSE)