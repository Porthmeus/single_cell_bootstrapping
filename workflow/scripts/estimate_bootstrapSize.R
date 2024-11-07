# Porthmeus
# 02.10.24

# logs
log.out <- file(snakemake@log[["out"]], open = "wt")
log.err <- file(snakemake@log[["err"]], open = "wt")
sink(log.out)
sink(log.err, type = "message")

# snakemake info
#print(snakemake)

# libraries
library(Seurat)
library(Matrix)
library(data.table)
library(parallel)

min_max_norm <- function(x){
    return((x-min(x))/(max(x)-min(x)))
}

expectedJaccard <- function(N,n){
    # N total population size of the cell
    # n subsample drawn from population
    n/(2*N)
}

cellPopVar <- function(x1, x2, pseudo_count = 0.05){
    # x1,x2 - matrix containing gene expression valies with genes in rows and cells in columns for both cell populations
    # pseudo_count - this number is added to the matrix in order to avoid division by 0 and log2 of 0
    if(!all(dim(x1) == dim(x2))){
        stop("Matrices have different dimensions")
    }
    if(pseudo_count <= 0){
        stop("pseudo_count needs to >0")
    }
    means1 <- rowMeans(x1) + pseudo_count
    means2 <- rowMeans(x2) + pseudo_count
    R = sum(abs(log2(means1/means2)))/nrow(x1)
    return(R)
}

meanCellPopVar <- function(x, draws, n = 30, pseudo_count = 0.05, replace = FALSE){
    # x - the matrix of the cell population which should be sampled. Rows are genes, columns are cells
    # draws - how many cells should be drawn from the population. This integer must be smaller or equal ncol(x)/2, unless replace = TRUE.
    # n - how often should the sampling be done, to calculate the average mean cell population variance
    # pseudo_count - this number is added to the matrix in order to avoid division by 0 and log2 of 0
    # replace - should cells be unique to the two populations which are drawn, if false, draws might be larger than N/2

    if(draws > ncol(x)/2 & replace == FALSE){
        stop("Draws need to be smaller than half of the number of cells")
    }

    Rs <- rep(0.0, n)
    for(i in 1:n){
        sel <- split(sample(1:ncol(x),2*draws, replace = replace), c(1,2))
        Rs[[i]] <- cellPopVar(x[,sort(sel[[1]])], x[,sort(sel[[2]])], pseudo_count = pseudo_count)
    }
    R <- mean(Rs)
    return(R)
}
        

meanCellPopVar_curves <- function(x, steps = 100, ...){
    # x - the matrix of the cell population which should be sampled. Rows are genes, columns are cells
    # steps - the number of steps the cell population variance should be calculated with different randomly drawn population sizes. If an integer, the maximum population size is split into equally distributed numbers. Use a vector to define specific numbers of population size which should be sampled
    # ... - pass on variables to meanCellPopVar
    
    if(length(steps) == 1){
        steps_length <- steps
        steps <- unique(round(seq(0,ncol(x)/2,length.out = steps)[-1]))
        steps <- steps[steps>2 & steps <= ncol(x)/2]
    } else {
        steps_length <- length(steps)
        steps<- unique(steps[steps>2&steps <= ncol(x)/2])
    }

    if(steps_length > length(steps)){
        warning(paste0("Reduced sampling to ", length(steps)," steps, because there were to few cells in the sampling pool, or given steps were smaller than 2 or larger than ncol(x)/2"))
    }


    stepped_R <-rep(0.0, length(steps))
    for(i in 1:length(steps)){
        st <- steps[i]
        stepped_R[i] <- meanCellPopVar(x = x, draws = st, ...)
    }

    dat <- data.frame(draws = steps,
                      R = stepped_R,
                      relative_R=min_max_norm(stepped_R),
                      total_pop = ncol(x))
    return(dat)
}

# set variables
stratification = snakemake@params[["stratification"]]
celltype = snakemake@params[["celltype"]]
alpha = snakemake@params[["alpha"]]
min_cell_pop <- snakemake@params[["min_cell_pop"]]
threads = snakemake@threads

# get the data
sce <- readRDS(snakemake@input[["sce"]])
mat <- sce@assays$RNA@counts
# split the matrix by celltype x sample
sets <- split(1:ncol(mat), f = paste0(sce@meta.data[[stratification]],"___",sce@meta.data[[celltype]]))
sets_sel <- sets[sapply(sets, length) > min_cell_pop]

# estimate the bootstrapping value
cl <- makeForkCluster(threads)
cellPopVars <-parLapply(cl = cl, sets_sel, function(x) meanCellPopVar_curves(mat[,x]))
stopCluster(cl)

# sort the data of the estimation
cellPopVars <- lapply(names(cellPopVars), function(x) {cbind(
                                             stratification = gsub("^(.*)___(.*)$","\\1", x),
                                             celltype = gsub("^(.*)___(.*)$","\\2", x),
                                             cellPopVars[[x]])}
)
cellPopVars <- data.table(do.call(rbind, cellPopVars))
colnames(cellPopVars)[1:2] <-c(stratification, celltype) 

# save the variance estimation
fwrite(cellPopVars, file = snakemake@output[["stats"]])

# get the bootstrap values per sample and celltype
drawsPerPop <- cellPopVars[relative_R <=alpha,.(bootstraps = min(draws), total_pop = unique(total_pop)), by = .(get(stratification), get(celltype))]
drawsPerPop[,expected_Jaccard := expectedJaccard(N = total_pop, n = bootstraps)]
drawsPerPop[,suggested_k := round(1/expected_Jaccard)]
colnames(drawsPerPop)[c(1,2)] <- c(stratification, celltype)
fwrite(drawsPerPop, file = snakemake@output[["bootstraps"]])

#### test ###
#stratification <- "Sample"
#celltype <- "CellType"
#min_cell_pop <- 30
#alpha <- 0.05
#threads <- 2
#
#sce <- readRDS("results/QC/SampleData_qced.rds")
#a <- lapply(sets_sel[1:5], function(x) meanCellPopVar_curves(mat[,x]))
#cellPopVars <- a
