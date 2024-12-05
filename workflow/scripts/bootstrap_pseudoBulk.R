# Porthmeus
# 04.10.24

# logs
log.out <- file(snakemake@log[["out"]], open = "wt")
log.err <- file(snakemake@log[["err"]], open = "wt")
sink(log.out)
sink(log.err, type = "message")

# libraries
require(Matrix)
require(Seurat)
require(data.table)
require(parallel)

sum_log <- function(x){
    # function to sum a matrix x by rows, if this matrix is log2 normalized
    log2(rowSums(2^x))
}

mean_log <- function(x){
    # function to take the mean of a matrix x by rows, if this matrix is log2 normalized
    log2(rowMeans(2^x))
}


# functions for bootstrapping
bootstrap_pseudoBulk<-function(sce, strat, strat_val, celltype, celltype_val,boots, k = 5, FUN = sum_log, matrix_slot= "data"){
    # sce - seurat object of the single cell, will assess
    # strat - column in the metat data for the stratification of the dataset, usually the samples within the set
    # strat_val - which sample should be used in the stratification column
    # celltype - the clustering of the data, usually the cell calling
    # celltype_val - which celltype specifically should be selected
    # boots - number of bootstraps
    # k - number of final pseudobulks created by bootstrapping
    # FUN - function for aggregation of values from individual cells of the bootstrapping
    # matrix_slot - which slot of the sce object should be accessed - default: "data" (corresponds to log normalized reads)
    # value - returns a sparse matrix with k columns and genes as rows


    all_cells <- data.table(sce@meta.data, keep.rownames=TRUE)[get(strat) == strat_val & get(celltype) == celltype_val,rn]
    boot_cells <- lapply(1:k, sample,x = all_cells, size = boots)
    
    mat <- slot(sce@assays$RNA, matrix_slot)
    if(deparse(FUN) == deparse(sum)){
        boot_mat <-  Matrix(sapply(boot_cells, function(x) rowSums(mat[,x])), sparse = TRUE)
    } else if(deparse(FUN) == depars(mean)){
        boot_mat <-  Matrix(sapply(boot_cells, function(x) rowMeans(mat[,x])), sparse =  TRUE)
    } else if(deparse(FUN) == depars(sum_log)){
        boot_mat <-  Matrix(sapply(boot_cells, function(x) sum_log(mat[,x])), sparse =  TRUE)
    } else if(deparse(FUN) == depars(mean_log)){
        boot_mat <-  Matrix(sapply(boot_cells, function(x) mean_log(mat[,x])), sparse =  TRUE)
    } else {
        boot_mat <- Matrix(sapply(boot_cells, function(x) apply(mat[,x],1, FUN)), sparse = TRUE)
    }
    colnames(boot_mat) <- paste(strat_val, celltype_val, 1:k, sep = "_")
    return(boot_mat)
}

# variables
stratification <-snakemake@config[["stratification"]]
celltype <-snakemake@config[["celltype"]]

# load seurat object and the table for the bootstrapping
sce <- readRDS(snakemake@input[["sce"]])
bs_table <- fread(snakemake@input[["bootstraps"]])

# do the actual bootstrapping and save the results
cl <- makeForkCluster(snakemake@threads)
pseudoBulk <- parSapply(cl,
                        #1:6,
                        1:nrow(bs_table),
                        function(x) {
                            bootstrap_pseudoBulk(sce = sce, 
                                                 strat = stratification,
                                                 strat_val = bs_table[x,get(stratification)],
                                                 celltype = celltype,
                                                 celltype_val = bs_table[x,get(celltype)],
                                                 boots = bs_table[x, bootstraps],
                                                 k = bs_table[x,suggested_k],
                                                 FUN = sum)})
stopCluster(cl)

pseudoBulk <- do.call(cbind, pseudoBulk)
writeMM(pseudoBulk, file = snakemake@output[["matrix"]])

# create a meta data file for the data
keep_cols <- c(stratification, celltype, "suggested_k")
meta_new <- bs_table[,..keep_cols]
by_cols <- c(stratification, celltype)
meta_new <- meta_new[,.(pseudoBulk_sample = 1:suggested_k,rn = paste(get(stratification), get(celltype), 1:suggested_k, sep = "_")), by =mget(by_cols)]
setkey(meta_new, "rn")
meta_new <- meta_new[colnames(pseudoBulk),]
fwrite(meta_new, file = snakemake@output[["meta"]])
