# Porthmeus
# 04.10.24

library(data.table)

# get the file to be read
arguments <-commandArgs(trailingOnly=TRUE)
print(arguments)
fl <- arguments[1]
sce <- readRDS(fl)
meta <- data.table(sce@meta.data)
stratCells <- meta[,.(PopSize =.N), by = .(get(arguments[2]), get(arguments[3]))]
colnames(stratCells)[1:2] <- c(arguments[2], arguments[3])
fwrite(stratCells, file = file.path("resources",
                                    paste0(tools::file_path_sans_ext(basename(arguments[1])),".csv")))

