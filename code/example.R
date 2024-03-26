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

## Paramater
para_table <- cross_validation(bulk = bulkData$Indata$T,
                                ref = bulkData$Indata$C_ref,
                                n_folds = 5,
                                seedd = 1234)
#you can chose other strategy, "which.min(para_table$RMSE.T)" especially when performing complete deconvolution
lambda1 <- para_table$lambda1[which.max(para_table$PCC.C)]
lambda2 <- para_table$lambda2[which.max(para_table$PCC.C)]
lambdaC <- para_table$lambdaC[which.max(para_table$PCC.C)]

## Deconvolution
retult <- CDSC(data_bulk = bulkData$Indata$T, 
               data_ref = bulkData$Indata$C_ref,
               lambda1 = lambda1, 
               lambda2 = lambda2,
               lambdaC = lambdaC,
               Ss = SM(t(scData$Indata$T)),
               Sg = SM(scData$Indata$T))
ctlabels <- Row_label(bulkData$Indata$C_ref,retult$c,leastnum = 3)
rownames(retult$p) <- ctlabels
colnames(retult$c) <- ctlabels

## Evalution
# please ensure consistency in cell type
getPearsonRMSE(retult$p, bulkData$Indata$P)
