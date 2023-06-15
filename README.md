# CDSC: complete deconvolution from bulk gene expression by leveraging single-cell RNA-Seq data

## 1. Introduction

CDSC is a complete deconvolution algorithm based on NMF. It can be used to infer cell-type-specific gene expression profiles (GEP) and cell-type proportions.

Given the bulk matrix and the reference GEP matrix, CDSC can compute cell similarity matrix and gene similarity matrix from the bulk matrix. Then, using the reference GEP matrix, CDSC performs deconvolution to calculate the cell type proportions and cell type-specific gene expression profiles of heterogeneous samples.

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.8000867

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
    #example.R
    setwd(path)
    source('CDSC.R')
    source('CDSC_help.R')
    example <- list(data = readRDS("example.rds"), full_phenoData = readRDS("example_phenoData.rds"))
    example$simulate1 <- scSimulateSplit(Segerstolpe,
                                         leastNum=50, plotmarker = F,
                                         norm1 = "CPM",log2.threshold = 1)
    example <- scSimulateShift(Segerstolpe,"all",standardization=TRUE)
    
    result = CDSC_3(data_bulk = scData$Indata$T, 
                data_ref = scData$Indata$C_ref, 
                k = dim(scData$Indata$C_ref)[2], 
                lambda1 = 1e-03, 
                lambda2 = 0e+00,
                lambdaC = 1000,
                Ss = SM(t(scData$Indata$T)),
                Sg = SM(scData$Indata$T))
                
    saveRDS(result,file="example.rds")                 
