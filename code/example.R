## load libraries (install via install.R)
library(dplyr)
library(Matrix)
library(limma)
library(pheatmap)
library(edgeR)
library(scater)
library(gtools)
library(reshape2)
library(data.table)
library(scales)

## Source
source("DSSC.R")
source("function_help.R")

## Prepare data
# If one generates simulated bulk data from single-cell data, please use:
scData <- list(data = readRDS("./Segerstolpe.rds"),
               full_phenoData = readRDS("./Segerstolpe_phenoData.rds"))
bulkData <- simulation(scData)

# If one analyzes true data, please use:
#bulkData <- trueData
# bulkData <- readRDS("CellLines.rds")

# The following code is available if you want to build the bulkData yourself:
# You can calculate the reference GEP matrix from scData
#C_ref <- scSimulateC(scData, leastNum = 0, plotmarker = FALSE, norm1 = "none",
#                     log2.threshold=log2(2))$C
# or you can also use your signature
#C_ref <- 'your signature'
#bulkData$Indata <- list(T = 'your bulk data',
#                        C_ref = C_ref,
#                        P = 'your groundtruth')
# To use DSSC correctly, please ensure that the genes of T and C_ref are the same

## Paramater

# Complete parameter adjustment process
# 
# para_table <- cross_validation(bulk = bulkData$Indata$T,
#                                ref = bulkData$Indata$C_ref,
#                                n_folds = 5,
#                                seedd = 1234)


# =============================Quick cross-validation example for testing DSSC =========================================
###Time difference of 56.12734 mins
para_table <- cross_validation(bulk = bulkData$Indata$T,
                                ref = bulkData$Indata$C_ref,
                                n_folds = 2,
                                seedd = 1234,lambda1 = c(0,0.001),lambda2 = c(0,0.001),lambdaC = c(1000,0))
#you can chose other strategy, "which.min(para_table$RMSE.T)" especially when performing complete deconvolution
lambda1 <- para_table$lambda1[which.max(para_table$PCC.C)]
lambda2 <- para_table$lambda2[which.max(para_table$PCC.C)]
lambdaC <- para_table$lambdaC[which.max(para_table$PCC.C)]

## Deconvolution
result <- DSSC(data_bulk = bulkData$Indata$T,
               data_ref = bulkData$Indata$C_ref,
               lambda1 = lambda1,
               lambda2 = lambda2,
               lambdaC = lambdaC,
               Ss = SM(t(bulkData$Indata$T)),
               Sg = SM(bulkData$Indata$T))

# ===================================================================================================================
# =============================Quick start for DSSC deconvolution====================================================


result <- DSSC(data_bulk = bulkData$Indata$T,
               data_ref = bulkData$Indata$C_ref,
               lambda1 = 0,
               lambda2 = 0,
               lambdaC = 1000,
               Ss = SM(t(bulkData$Indata$T)),
               Sg = SM(bulkData$Indata$T))


# =================================================================================================================



ctlabels <- Row_label(result$c,bulkData$Indata$C_ref, leastnum = 3)
rownames(result$p) <- ctlabels
colnames(result$c) <- ctlabels

## Evalution
# Please ensure consistency in cell type
getPearsonRMSE(result$p, bulkData$Indata$P)
getPearsonRMSE(result$c, bulkData$Indata$C)
