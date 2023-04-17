rm(list=ls())
# works with R version 4.2.3

# supress warnings
options(warn=-1)
# Load libraries
suppressMessages(library(dplyr))
library(tidyr)
library(openxlsx)
# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# first argument is the input file pathway
input_file <- args[1]
# Load data
data <- read.table(input_file, header = TRUE, sep = "\t")

# Parse data
# group by 3rd column and list all 4th column values
sc_type_format <- data %>% group_by(Celltype) %>% summarise("MarkerGene" = paste(MarkerGene, collapse = ","))
# add dataset column and empty negative markers column
print("before")
sc_type_format <- cbind("Dataset" = data$Dataset, sc_type_format)
print("here")
sc_type_format <- cbind(sc_type_format, data.frame("NegativeMarkerGene" = "" ))
# rename columns according to sc-type format
colnames(sc_type_format) <- c("tissueType", "cellName", "geneSymbolmore1", "geneSymbolmore2")

# save output
# make new name by removing extension and adding -sc-type to the end of the input file name
output_file <- paste0(sub('\\...*$','', input_file), "-sc-type.xlsx")
# write to xlsx
write.xlsx(sc_type_format,
           output_file,
           sheetName = "Sheet1", rowNames = FALSE)

print("Done!")
print(paste0("Output file: ", output_file))