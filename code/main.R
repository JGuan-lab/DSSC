## Source
rm (list=ls ())
source("CDSC.R")
source("function_help.R")

## Read data
scData <- list(data = readRDS("XXX.rds"), full_phenoData = readRDS("XXX_phenoData.rds"))

##Simulation
bulkData <- simulation(scData)
# bulkData <- trueData

## Deconvolution
retult <- CDSC(data_bulk = bulkData$Indata$T, 
               data_ref = data_bulk$Indata$C_ref,  
               k = dim(data_bulk$Indata$C_ref)[2], 
               lambda1 = 1e-03, 
               lambda2 = 0e+00,
               lambdaC = 1000,
               Ss = SM(t(scData$Indata$T)),
               Sg = SM(scData$Indata$T))
## Evalution
ctlabels <- Row_label(data_bulk$Indata$C_ref,retult$c,leastnum = 3)
rownames(retult$p) <- ctlabels
colnames(retult$c) <- ctlabels
getPearsonRMSE(retult$p, data_bulk$Indata$P)
