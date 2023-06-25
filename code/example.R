## Source
source("CDSC.R")
source("function_help.R")

## Prepare data
# if one generates simulated bulk data from single-cell data, please use:
scData <- list(data = readRDS("Segerstolpe.rds"), full_phenoData = readRDS("Segerstolpe_phenoData.rds"))
bulkData <- simulation(scData)

# if one analyzes true data, please use:
bulkData <- trueData 
# bulkData <- readRDS("CellLines.rds")

# The following code is available if you want to build the bulkData yourself:
# You can calculate the reference GEP matrix from scData
C_ref <- scSimulateC(scData,leastNum = 0,plotmarker = F,norm1 = "none",log2.threshold=log2(2))$C
# or you can alse use your signature
C_ref <- 'your signature'
bulkData$Indata <- list(T = 'your bulk data',
                        C_ref = C_ref,
                        P = 'your groundtruth')
# To use CDSC correctly, please ensure that the genes of T and C_ref are the same

## Deconvolution
retult <- CDSC(data_bulk = bulkData$Indata$T, 
               data_ref = bulkData$Indata$C_ref,  
               k = dim(bulkData$Indata$C_ref)[2], 
               lambda1 = 1e-03, 
               lambda2 = 0e+00,
               lambdaC = 1000,
               Ss = SM(t(scData$Indata$T)),
               Sg = SM(scData$Indata$T))
ctlabels <- Row_label(bulkData$Indata$C_ref,retult$c,leastnum = 3)
rownames(retult$p) <- ctlabels
colnames(retult$c) <- ctlabels

## Evalution
# please ensure consistency in cell type
getPearsonRMSE(retult$p, bulkData$Indata$P)
