# DSSC: deconvolution from bulk gene expression by leveraging similarity and single-cell RNA-seq data

## 1. Introduction

DSSC is a deconvolution algorithm, which estimates cell type-specific gene expression profiles (GEPs) and cell type density simultaneously from bulk samples by leveraging single-cell gene expression data. 

Given bulk expression matrix and referenced single-cell gene expression data (or referenced GEPs), DSSC computes sample-sample similarity matrix and gene-gene similarity matrix from the bulk data, and calculates the averaged gene expression of each cell type for reference from the single-cell data (or uses the referenced GEPs directly). Then, DSSC performs deconvolution to infer cell type-specific GEPs and cell type proportions of heterogeneous samples.

The datasets analyzed in the paper are available at: https://doi.org/10.5281/zenodo.8020767

## 2. Quick start    
Depends:
    
    R (>= 4.1.0) 

Source:
    
    source("DSSC.R")
    source("function_help.R")

### 2.1 Prepare data
    # if you want to generate simulated bulk data from single-cell data, please use:
    scData <- list(data = readRDS("Segerstolpe.rds"), full_phenoData = readRDS("Segerstolpe_phenoData.rds"))
    bulkData <- simulation(scData)

    # if you want to analyze true data we provide, please use:
    bulkData <- readRDS("CellLines.rds")

    # if you want to use your own single-cell data to generate referenced GEP matrix, please use:
    C_ref <- scSimulateC(scData,leastNum = 0,plotmarker = F,norm1 = "none",log2.threshold=log2(2))$C
    # or you can use a given referenced GEP:
    C_ref <- 'your signature'
    bulkData$Indata <- list(T = 'your bulk data',
                            C_ref = C_ref,
                            P = 'your groundtruth')
    # To use DSSC correctly, please ensure that the genes of T and C_ref are the same.

### 2.2 Paramater
    #Here, we recommend using a 5-fold cross validation parameter tuning method to obtain appropriate parameters. Of course, you can also change the settings of the tuning process according to the actual situation
    #bulk: the input bulk data
    #dref: the reference GEP matrix
    #k_folds: k_folds cross validatioin, set by default to 5
    #seedd: random seed, set by default to 1234
    #TerCondition: iteration termination condition, set by default to 10^-8,
    #lambda1: the range of parameter lambda1, set by default to (0,10^-3,10^-2,10^-1,1,10)
    #lambda2: the range of parameter lambda2, set by default to (0,10^-3,10^-2,10^-1,1,10)
    #lambdaC: the range of parameter lambdaC, set by default to (0,10^-1,10^0,10^1,10^2,10^3)
    #max_num: maximum number of iterations, default set to 1500
    para_table <- cross_validation(bulk = bulkData$Indata$T,
                                ref = bulkData$Indata$C_ref,
                                k_folds = 5,
                                seedd = 1234)
    #you can chose other strategy, like "which.max(para_table$PCC.T)"
    lambda1 <- para_table$lambda1[which.min(para_table$RMSE.C)]
    lambda2 <- para_table$lambda2[which.min(para_table$RMSE.C)]
    lambdaC <- para_table$lambdaC[which.min(para_table$RMSE.C)]

### 2.3 Deconvolution
    
    #data_bulk: the input bulk data
    #data_ref: the reference GEP matrix
    #lambda1 is used to constrain the sample-sample similarity matrix
    #lambda2 is used to constrain the gene-gene similarity matrix
    #lambdaC is used to constrain the GEP matrix
    #k: number of cell types used for matrix factorization initialization
    
    retult <- DSSC(data_bulk = bulkData$Indata$T,
                   data_ref = bulkData$Indata$C_ref,  
                   lambda1 = lambda1, 
                   lambda2 = lambda2,
                   lambdaC = lambdaC,
                   Ss = SM(t(scData$Indata$T)),
                   Sg = SM(scData$Indata$T))
    ctlabels <- Row_label(bulkData$Indata$C_ref,retult$c,leastnum = 3)
    rownames(retult$p) <- ctlabels
    colnames(retult$c) <- ctlabels
                   
### 2.4 Evaluation

    getPearsonRMSE(retult$p, data_bulk$Indata$P)
