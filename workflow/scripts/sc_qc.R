# Porthmeus
# 30.09.24

# quality control for single cell data
# logs
log.out <- file(snakemake@log[["out"]], open = "wt")
log.err <- file(snakemake@log[["err"]], open = "wt")
sink(log.out)
sink(log.err, type = "message")

# load packages
library(Seurat)
library(data.table)


isOutlier <- function(x,nmad =3){
    mad.x <- mad(x)
    out <- abs(x-median(x)) > nmad*mad.x 
    return(out)
}



# params
nmad <- snakemake@params[["nmad"]]
strat <- snakemake@params[["stratification"]]
celltype <- snakemake@params[["celltype"]]

# load data
sce <- readRDS(snakemake@input[[1]])

# do the filtering
metadata <- data.table(sce@meta.data, keep.rownames = TRUE)
metadata[,outlier.sampleCell := (isOutlier(log(nCount_RNA),nmad=nmad)|
                  isOutlier(log(nFeature_RNA),nmad = nmad) |
                  isOutlier(percent.mt,nmad=nmad)), by = .(get(strat), get(celltype))]
sce@meta.data[["outlier.sampleCell"]] <- metadata[,outlier.sampleCell]
sce <- sce[,!sce@meta.data[["outlier.sampleCell"]]]
saveRDS(sce, snakemake@output[["dataset"]])

# calculate stats
stat <- metadata[,.(
                    min_log_nCount_RNA = median(log(nCount_RNA)) - mad(log(nCount_RNA)),
                    max_log_nCount_RNA = median(log(nCount_RNA)) + mad(log(nCount_RNA)),
                    min_log_nFeature_RNA = median(log(nFeature_RNA)) - mad(log(nFeature_RNA)),
                    max_log_nFeature_RNA = median(log(nFeature_RNA)) + mad(log(nFeature_RNA)),
                    min_percent.mt = median(percent.mt) - mad(percent.mt),
                    max_percent.mt = median(percent.mt) + mad(percent.mt),
                    cells_removed = sum(outlier.sampleCell),
                    cells_kept = sum(!outlier.sampleCell)), by = .(stratification = get(strat), Cell = get(celltype))]
fwrite(stat, snakemake@output[["thresholds"]])

