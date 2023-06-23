# CDSC: complete deconvolution from bulk gene expression by leveraging single-cell RNA-seq data

## 1. Introduction

CDSC is a complete deconvolution algorithm, which estimates cell type-specific gene expression profiles (GEPs) and cell type density simultaneously from bulk samples by leveraging single-cell gene expression data. 

Given bulk expression matrix and referenced single-cell gene expression data (or referenced GEPs), CDSC can compute cell-cell similarity matrix and gene-gene similarity matrix from the bulk matrix, and calculate the averaged gene expression of each cell type for reference from the single-cell data (or use the referenced GEPs directly). Then, CDSC performs deconvolution to infer cell type-specific GEPs and  cell type proportions of heterogeneous samples.

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.8020767

## 2. Quick start    
Depends:
    
    R (>= 4.1.0) 

Source:
    
    rm (list=ls ())
    source("CDSC.R")
    source("function_help.R")

### 2.1 Prepare data
    # if one generates simulated bulk data from single-cell data, please use:
    scData <- list(data = readRDS("XXX.rds"), full_phenoData = readRDS("XXX_phenoData.rds"))
    bulkData <- simulation(scData)
    
    # if one analyzes true data, please use:
    # bulkData <- trueData 

### 2.2 Deconvolution
Input data:
    data_bulk, the input bulk data;
    data_ref, the reference GEP matrix.

The parameters used in SCDC:    
    lambda1 is used to constrain the sample-sample similarity matrix;
    lambda2 is used to constrain the gene-gene similarity matrix;
    lambdaC is used to constrain the GEP matrix;
    k: number of cell types used for matrix factorization initialization.
    
    retult <- CDSC(data_bulk = bulkData$Indata$T,
                   data_ref = data_bulk$Indata$C_ref,  
                   k = dim(data_bulk$Indata$C_ref)[2], 
                   lambda1 = 1e-03, 
                   lambda2 = 0e+00,
                   lambdaC = 1000,
                   Ss = SM(t(scData$Indata$T)),
                   Sg = SM(scData$Indata$T))
### 2.3 Evalution

    ctlabels <- Row_label(data_bulk$Indata$C_ref,retult$c,leastnum = 3)
    rownames(retult$p) <- ctlabels
    colnames(retult$c) <- ctlabels
    getPearsonRMSE(retult$p, data_bulk$Indata$P)
