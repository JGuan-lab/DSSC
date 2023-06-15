# CDSC: complete deconvolution from bulk gene expression by leveraging single-cell RNA-Seq data

## 1. Introduction

CDSC is a complete deconvolution algorithm based on NMF. It can be used to infer cell-type-specific gene expression profiles (GEP) and cell-type proportions.

Given the bulk matrix and the reference GEP matrix, CDSC can compute cell similarity matrix and gene similarity matrix from the bulk matrix. Then, using the reference GEP matrix, CDSC performs deconvolution to calculate the cell type proportions and cell type-specific gene expression profiles of heterogeneous samples.

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.8020767

## Running the tests

### Depends:
    R (>= 4.1.0) 
### Input data:
    data_bulk: the input dropout bulk data.

    data_ref: the reference GEP matrix
### The parameter used in SCDC:

    #In this article, the parameters can be selected through grid search.
    
    parameter the vector of parameters. 
    lambda1 is the value of lambda1 in the mathematical model to limit the sample-sample similarity matrix;
    lambda2 is the value of lambda2 in the mathematical model to limit the gene-gene similarity matrix;
    lambdaC is the value of lambda3 in the mathematical model to limit the GEP matrix.
    
    k: Number of cell types used for matrix initialization，it can be obtained from the reference GEP matrix or manually inputted.
    
### Example:
#### Source
    rm (list=ls ())
    source("CDSC.R")
    source("function_help.R")

#### Read data
    scData <- list(data = readRDS("XXX.rds"), full_phenoData = readRDS("XXX_phenoData.rds"))

#### Simulation
    bulkData <- simulation(scData)
    # if you need true data, please insert here directly.
    bulkData <- trueData 

#### Deconvolution

    retult <- CDSC(data_bulk = bulkData$Indata$T,
                   data_ref = data_bulk$Indata$C_ref,  
                   k = dim(data_bulk$Indata$C_ref)[2], 
                   lambda1 = 1e-03, 
                   lambda2 = 0e+00,
                   lambdaC = 1000,
                   Ss = SM(t(scData$Indata$T)),
                   Sg = SM(scData$Indata$T))
#### Evalution

    ctlabels <- Row_label(data_bulk$Indata$C_ref,retult$c,leastnum = 3)
    rownames(retult$p) <- ctlabels
    colnames(retult$c) <- ctlabels
    getPearsonRMSE(retult$p, data_bulk$Indata$P)
