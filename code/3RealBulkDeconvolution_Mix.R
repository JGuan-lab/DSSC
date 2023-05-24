# real bulk data sets and sc data , find the deconvolution results

# ---
# title: "3RealBulkDeconvolution"
# author: "Chenqi Wang"
# date: "2022/4/11"
# output: html_document
# ---
# 

#-----hunman pancreas ///-------------

##
#------read raw data information--------
rm (list=ls ())
setwd("F:/wangchenqi/CDSC/3realBulk")

# list.files("F:/wangchenqi/CDSC/3realBulk/data/")

#-------------------------------
# CIBERSORT 提到的Melanoma数据和Wholeblood数据，可以参考CDSeq去使用模拟数据和混合数据（一共4组）

# library(GEOquery)
# BreastBlood <- list(
#   T = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/BreastBlood_GSE29832/mix.txt",row.names = 1,header = T),
#   C = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/BreastBlood_GSE29832/sig.txt",row.names = 1,header = T,sep = "\t"),
#   P = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/BreastBlood_GSE29832/coef.txt",row.names = 1,header = T,sep = "\t"),
#   data = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/BreastBlood_GSE29832/pure.txt",row.names = 1,header = T,sep = "\t"),
#   full_phenoData = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/BreastBlood_GSE29832/pure_annotations.txt",header = T,sep = "\t"),
#   ShiftGene = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/BreastBlood_GSE29832/ShiftGene.txt",sep = "\t"),
#   CT = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/BreastBlood_GSE29832/CT.txt",sep = "\t")
#   )
# GSE29832 <- getGEO("GSE29832")
# # 
# CellLines <- list(
#   T = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/CellLines_GSE11058/mix.txt"),
#   C = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/CellLines_GSE11058/sig.txt",row.names = 1,header = T,sep = "\t"),
#   P = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/CellLines_GSE11058/coef.txt",row.names = 1,header = T,sep = "\t"),
#   data = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/CellLines_GSE11058/pure.txt",row.names = 1,header = T,sep = "\t"),
#   full_phenoData = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/CellLines_GSE11058/pure_annotations.txt",header = T,sep = "\t"),
#   ShiftGene = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/CellLines_GSE11058/ShiftGene.txt",sep = "\t"),
#   CT = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/CellLines_GSE11058/CT.txt",sep = "\t"),
# 
# )
# # GSE11058 <- getGEO("GSE11058")
# 
# LiverBrainLung <- list(
#   T = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/LiverBrainLung_GSE19830/mix.txt",row.names = 1,header = T),
#   C = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/LiverBrainLung_GSE19830/sig.txt",row.names = 1,header = T,sep = "\t"),
#   P = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/LiverBrainLung_GSE19830/coef.txt",row.names = 1,header = T,sep = "\t"),
#   data = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/LiverBrainLung_GSE19830/pure.txt",row.names = 1,header = T,sep = "\t"),
#   full_phenoData = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/LiverBrainLung_GSE19830/pure_annotations.txt",header = T,sep = "\t"),
#   ShiftGene = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/LiverBrainLung_GSE19830/ShiftGene.txt",sep = "\t"),
#   CT = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/LiverBrainLung_GSE19830/CT.txt",sep = "\t")
#   
# )
# # GSE19830<- getGEO("GSE19830")
# 
# RatBrain<- list(
#   T = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/RatBrain_GSE19380/mix.txt",row.names = 1,header = T),
#   C = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/RatBrain_GSE19380/sig.txt",row.names = 1,header = T,sep = "\t"),
#   P = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/RatBrain_GSE19380/coef.txt",row.names = 1,header = T,sep = "\t"),
#   data = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/RatBrain_GSE19380/pure.txt",row.names = 1,header = T,sep = "\t"),
#   full_phenoData = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/RatBrain_GSE19380/pure_annotations.txt",header = T,sep = "\t"),
#   ShiftGene = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/RatBrain_GSE19380/ShiftGene.txt",sep = "\t"),
#   CT = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/RatBrain_GSE19380/CT.txt",sep = "\t")
#   
# )
# # 
# # GSE19380 <- getGEO("GSE19380")
# #
# Retina <- list(
#   T = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/Retina_GSE33076/mix.txt",row.names = 1,header = T),
#   C = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/Retina_GSE33076/sig.txt",row.names = 1,header = T,sep = "\t"),
#   P = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/Retina_GSE33076/coef.txt",row.names = 1,header = T,sep = "\t"),
#   data = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/Retina_GSE33076/pure.txt",row.names = 1,header = T,sep = "\t"),
#   full_phenoData = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/Retina_GSE33076/pure_annotations.txt",header = T,sep = "\t"),
#   ShiftGene = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/Retina_GSE33076/ShiftGene.txt",sep = "\t"),
#   CT = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/Retina_GSE33076/CT.txt",sep = "\t")
#   
# )
# GSE33076 <- getGEO("GSE33076")

# ExperimentalMixtures <- list(
#   data = read.table("F:/wangchenqi/CDSC/3realBulk/data/input/ExperimentalMix_GSE123604/mix.txt",header = T)
#   
# )
# library(GEOquery)
# ann<- getGEO(filename = "F:/wangchenqi/CDSC/3realBulk/data/input/ExperimentalMix_GSE123604/GSE123604_series_matrix.txt.gz", getGPL = F);ann <- ann@phenoData@data

# GSE123604 <- GEOquery::getGEO("GSE123604 ")
#------simulate--------
source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

# ----------------BreastBlood-----------------
  # 换名字
# colnames(BreastBlood$data) <- BreastBlood$full_phenoData$Class
# BreastBlood$markers <- Find_markerGene_limma(BreastBlood$data)
# BreastBlood$ShiftGene <- BreastBlood$markers$gene

# all(rownames(BreastBlood$T) == rownames(BreastBlood$C) )
# all(rownames(BreastBlood$T) == rownames(BreastBlood$data))
# 
# BreastBlood$ShiftGene <- as.matrix(BreastBlood$ShiftGene)
# BreastBlood$CT <- as.matrix(BreastBlood$CT)
# 
# length(which(rownames(BreastBlood$T) %in% BreastBlood$ShiftGene))
# 
# BreastBlood$indata$T = as.matrix(BreastBlood$T[which(rownames(BreastBlood$T) %in% BreastBlood$ShiftGene),])
# BreastBlood$indata$C = as.matrix(BreastBlood$C[which(rownames(BreastBlood$C) %in% BreastBlood$ShiftGene),])
# BreastBlood$indata$C_ref = BreastBlood$indata$C
# 
# BreastBlood$indata$P = as.matrix(BreastBlood$P)
# BreastBlood$indata$gene = BreastBlood$ShiftGene
# 
# # BreastBlood <- Deal7DataSet(BreastBlood)
# saveRDS(BreastBlood,"F:/wangchenqi/CDSC/3realBulk/bulkdata/BreastBlood.rds")
# 
# #----------------CellLines-----------------
# #换名字
# colnames(CellLines$data) <- CellLines$full_phenoData$Class
# CellLines$markers <- Find_markerGene_limma(CellLines$data)
# CellLines$ShiftGene <- CellLines$markers$gene
# 
# all(rownames(CellLines$T) == rownames(CellLines$C) )
# all(rownames(CellLines$T) == rownames(CellLines$data))
# 
# CellLines$ShiftGene <- as.matrix(CellLines$ShiftGene)
# length(which(rownames(CellLines$T) %in% CellLines$ShiftGene))
# CellLines$indata$T = as.matrix(CellLines$T[which(rownames(CellLines$T) %in% CellLines$ShiftGene),])
# CellLines$indata$C = as.matrix(CellLines$C[which(rownames(CellLines$C) %in% CellLines$ShiftGene),])
# CellLines$indata$C_ref = CellLines$indata$C
# CellLines$indata$P = as.matrix(CellLines$P)
# CellLines$indata$gene = CellLines$ShiftGene
# 
# CellLines <- Deal7DataSet(CellLines)
# saveRDS(bulkData,"F:/wangchenqi/CDSC/3realBulk/bulkdata/CellLines.rds")
# 
# #----------------LiverBrainLung-----------------
# #换名字
# all(rownames(LiverBrainLung$T) == rownames(LiverBrainLung$C) )
# all(rownames(LiverBrainLung$T) == rownames(LiverBrainLung$data))
# LiverBrainLung$ShiftGene <- as.matrix(LiverBrainLung$ShiftGene)
# length(which(rownames(LiverBrainLung$T) %in% LiverBrainLung$ShiftGene))
# LiverBrainLung$indata$T = as.matrix(LiverBrainLung$T[which(rownames(LiverBrainLung$T) %in% LiverBrainLung$ShiftGene),])
# LiverBrainLung$indata$C = as.matrix(LiverBrainLung$C[which(rownames(LiverBrainLung$C) %in% LiverBrainLung$ShiftGene),])
# LiverBrainLung$indata$C_ref = LiverBrainLung$indata$C
# LiverBrainLung$indata$P = as.matrix(LiverBrainLung$P)
# LiverBrainLung$indata$gene = LiverBrainLung$ShiftGene
# 
# # LiverBrainLung <- Deal7DataSet(LiverBrainLung)
# saveRDS(LiverBrainLung,"F:/wangchenqi/CDSC/3realBulk/bulkdata/LiverBrainLung.rds")
# 
# #----------------RatBrain-----------------
# #换名字
# all(rownames(RatBrain$T) == rownames(RatBrain$C) )
# all(rownames(RatBrain$T) == rownames(RatBrain$data))
# RatBrain$ShiftGene <- as.matrix(RatBrain$ShiftGene)
# length(which(rownames(RatBrain$T) %in% RatBrain$ShiftGene))
# RatBrain$indata$T = as.matrix(RatBrain$T[which(rownames(RatBrain$T) %in% RatBrain$ShiftGene),])
# RatBrain$indata$C = as.matrix(RatBrain$C[which(rownames(RatBrain$C) %in% RatBrain$ShiftGene),])
# RatBrain$indata$C_ref = RatBrain$indata$C
# RatBrain$indata$P = as.matrix(RatBrain$P)
# RatBrain$indata$gene = RatBrain$ShiftGene
# 
# # RatBrain <- Deal7DataSet(RatBrain)
# saveRDS(RatBrain,"F:/wangchenqi/CDSC/3realBulk/bulkdata/RatBrain.rds")
# 
# #----------------Retina-----------------
# #换名字
# all(rownames(Retina$T) == rownames(Retina$C) )
# all(rownames(Retina$T) == rownames(Retina$data))
# Retina$ShiftGene <- as.matrix(Retina$ShiftGene)
# length(which(rownames(Retina$T) %in% Retina$ShiftGene))
# Retina$indata$T = as.matrix(Retina$T[which(rownames(Retina$T) %in% Retina$ShiftGene),])
# Retina$indata$C = as.matrix(Retina$C[which(rownames(Retina$C) %in% Retina$ShiftGene),])
# Retina$indata$C_ref = Retina$indata$C
# Retina$indata$P = as.matrix(Retina$P)
# Retina$indata$gene = Retina$ShiftGene
# 
# # Retina <- Deal7DataSet(Retina)
# saveRDS(Retina,"F:/wangchenqi/CDSC/3realBulk/bulkdata/Retina.rds")

#----------------ExperimentalMixtures-----------------
#换名字
# ExperimentalMixtures$pure <- ExperimentalMixtures$data[,c(1:4,44:47)]
# cellnames <- colnames(ExperimentalMixtures$data)
# WholeBlood$C$PBMCs_5_pData = NULL
# WholeBlood$C$PBMCs_5_pData = cellnames
# WholeBlood$C$PBMCs_5_pData <- cbind(WholeBlood$C$PBMCs_5_pData, gsub("\\.[0-9]*$","",cellnames))
# WholeBlood$C$PBMCs_5_pData <- cbind(WholeBlood$C$PBMCs_5_pData, 1)
# colnames(WholeBlood$C$PBMCs_5_pData) <- c("cellID",	"cellType",	"sampleID")
# head(WholeBlood$C$PBMCs_5_pData)
# 
# WholeBlood$C$PBMCs_5_sc <- as.matrix(WholeBlood$C$PBMCs_5)
# WholeBlood$C$PBMCs_5_pData <- as.data.frame(WholeBlood$C$PBMCs_5_pData)
# rownames(WholeBlood$C$PBMCs_5_pData) <- WholeBlood$C$PBMCs_5_pData$cellID
# 
# scData <- list(data = WholeBlood$C$PBMCs_5_sc, full_phenoData = WholeBlood$C$PBMCs_5_pData)
# scDataSimulate <- scSimulateC(scData,leastNum = 0,plotmarker = F)
# 
# WholeBlood$sig$PBMCs_5_sc <- scDataSimulate$C[scDataSimulate$markerslist$gene,]
# interaction_ct = intersect(colnames(WholeBlood$sig$PBMCs_5_sc),colnames(WholeBlood$P));interaction_ct
# WholeBlood$sig$PBMCs_5_sc <- WholeBlood$sig$PBMCs_5_sc[,interaction_ct]
# WholeBlood$GroundTruth$PBMCs_5_sc <- WholeBlood$P[,interaction_ct]
# 
# WholeBlood$indata$C <- WholeBlood$sig$PBMCs_5_sc
# WholeBlood$indata$P <- t(WholeBlood$GroundTruth$PBMCs_5_sc)
# WholeBlood$indata$T <- WholeBlood$T
# 
# gene <- intersect(rownames(WholeBlood$indata$C),rownames(WholeBlood$T));length(gene)
# WholeBlood$indata$T <- WholeBlood$indata$T[gene,]
# WholeBlood$indata$T <- WholeBlood$indata$T[which(rowSums(WholeBlood$indata$T) > 0),]
# length(which(rowSums(WholeBlood$indata$T) > 0))
# gene <- gene[which(rowSums(WholeBlood$indata$T) > 0)];length(gene)
# WholeBlood$indata$C <- WholeBlood$indata$C[gene,]
# 
# # ExperimentalMixtures <- Deal7DataSet(ExperimentalMixtures)
# saveRDS(ExperimentalMixtures,"F:/wangchenqi/CDSC/3realBulk/bulkdata/ExperimentalMixtures.rds")

# > refdata = bulkData$data
# > colnames(refdata) = bulkData[["full_phenoData"]][["Class"]]
# > View(refdata)
# > MARKERS = Find_markerGene_limma(refdata, plotmarker=F)
# > View(MARKERS)
# > marklist = which(MARKERS$FDR <= 0.1 & MARKERS$P.value<0.1)
# > marklist  = MARKERS[marklist,]

#---START deconvolution-----------------------------------
rm (list=ls ())
setwd("F:/wangchenqi/CDSC/3realBulk")

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/3realBulk/bulkData"));a
a <- c("BreastBlood","CellLines","LiverBrainLung","RatBrain","Retina" );a
STRING_name = a[3]; STRING_name

getwd()
bulkData <- readRDS(paste(getwd(),"/bulkData/",STRING_name,".rds",sep=""))
result <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))
# result$all

MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,"deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",
                  "MuSiC","Bisque","SCDC", "DWLS",
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
result <- CountAllResults(result,MyMethodName)
result$all
nrow(result$all)
STRING_name
# saveRDS(result,paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))

# GetCorMatrix(result$CDSC3$dec$p,bulkData$indata$P)
# diag(GetCorMatrix(result$CDSC3$dec$p,bulkData$indata$P,1))
# diag(cor(t(result$CIBERSORT$p),t(bulkData$indata$P)))
# diag(cor(t(result$ssKL$p),t(bulkData$indata$P)))
# diag(cor(t(result$DSA$p),t(bulkData$indata$P)))

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")
# bulkData <- WholeBlood
# saveRDS(result,paste(getwd(),"/result/result_",STRING_name,".rds",sep=""))
para_lambda_123 <- readRDS(paste(getwd(),"/paramater/para_lambda_",STRING_name,"_123.rds",sep=""))
para_lambda_12 <- para_lambda_123[para_lambda_123$lambdaC == 0,]
para_lambda_12 <- para_lambda_44_8[para_lambda_44_8$lambdaC ==0,]
#关于C的ct的热图
STRING_name
map = GetCorMatrix(result$CDSC3$dec$c,bulkData$indata$C,matrix = "c");map
library(pheatmap)
plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 22, cellheight = 22#,gaps_row = c(12, 17)
                         ,fontsize = 8
                         # ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.2f" # 数字
                         # ,border_color = "black",color = colorRampPalette(MyColor)(60)
                         # ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,annotation_legend = T
                         ,main = STRING_name)
eoffice::topptx(plot_pheatmp,filename = paste(
  "F:/wangchenqi/CDSC/pictures/class3_ct_",STRING_name,".pptx",sep=""))


#-------deconvolution---
result <- list()
result$CDSC3$seedd = 44
result$CDSC3$TerCondition = 10^-8
result$CDSC3$Ss <- SM(t(bulkData$indata$T))
result$CDSC3$Sg <- SM(bulkData$indata$T)

result$CDSC3$lambda1 <- 1e-03
result$CDSC3$lambda2 <- 1e+00
result$CDSC3$lambdaC <- 1e+01

library(dplyr)
# bulkData$indata$C <- bulkData$indata$C_ref[,-1]
# bulkData$indata$P <- bulkData$indata$P[-1,]
# bulkData$indata$C_ref <- bulkData$indata$C
result$CDSC3$dec = CDSC_3(bulkData$indata$T, bulkData$indata$C_ref, dim(bulkData$indata$C_ref)[2], 
                           result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                           result$CDSC3$TerCondition,result$CDSC3$seedd,
                           result$CDSC3$Ss, result$CDSC3$Sg,all_number = 3000)
result$CDSC3$dec$p
bulkData$indata$P
cor(t(result$CDSC3$dec$p),t(bulkData$indata$P))
rownames(result$CDSC3$dec$p)
result$CDSC3$result <- calculate_result(result$CDSC3$dec$c,result$CDSC3$dec$p,
                                         bulkData$indata$T, bulkData$indata$C, bulkData$indata$C_ref, bulkData$indata$P,
                                         result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                                         result$CDSC3$dec$jump, result$CDSC3$seedd, result$CDSC3$TerCondition)
result$CDSC3$result 

result$CDSC2$lambda1 <- 1e-03
result$CDSC2$lambda2 <- 1e+00
result$CDSC2$lambdaC <- 0
result$CDSC2$dec = CDSC_2(bulkData$indata$T,  dim(bulkData$indata$C)[2], 
                           result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                           result$CDSC3$TerCondition,result$CDSC3$seedd,
                           result$CDSC3$Ss,result$CDSC3$Sg,all_number = 3000)
cor(t(result$CDSC2$dec$p),t(bulkData$indata$P))
result$CDSC2$result <- calculate_result(result$CDSC2$dec$c,result$CDSC2$dec$p,
                                         bulkData$indata$T, bulkData$indata$C, bulkData$indata$C_ref, bulkData$indata$P,
                                         result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                                         result$CDSC2$dec$jump, result$CDSC3$seedd,result$CDSC3$TerCondition)
result$CDSC2$result

#----------NNLS----------
require(NNLS)
result$NNLS$p = do.call(cbind.data.frame,lapply(apply(bulkData$indata$T,2,function(x) nnls::nnls(as.matrix(bulkData$indata$C),x)),   function(y) y$x))
result$NNLS$p = apply(result$NNLS$p,2,function(x) x/sum(x)) #explicit STO constraint
cor(t(result$NNLS$p),t(bulkData$indata$P))
rownames(result$NNLS$p) <- colnames(bulkData$indata$C)
result$NNLS$result = getPearsonRMSE(result$NNLS$p,bulkData$indata$P)
result$NNLS$result

#--------OLS------------
result$OLS$p = apply(bulkData$indata$T,2,function(x) lm(x ~ as.matrix(bulkData$indata$C))$coefficients[-1])
result$OLS$p = apply(result$OLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$OLS$p = apply(result$OLS$p,2,function(x) x/sum(x)) #explicit STO constraint
cor(t(result$OLS$p),t(bulkData$indata$P))
rownames(result$OLS$p) <- unlist(lapply(strsplit(rownames(result$OLS$p),")"),function(x) x[2]))
result$OLS$result = getPearsonRMSE(result$OLS$p,bulkData$indata$P)
result$OLS$result


#---------FARDEEP-----------
library(FARDEEP)
#result_FARDEEP = t(FARDEEP(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = t(FARDEEP::fardeep(bulkData$indata$C, bulkData$indata$T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = apply(result$FARDEEP$p,2,function(x) x/sum(x)) #explicit STO constraint
result$FARDEEP$result = getPearsonRMSE(result$FARDEEP$p,bulkData$indata$P)
result$FARDEEP$result

#----------CIBERSORT-----------------
source("F:/wangchenqi/CDSC/CIBERSORT.R")
result$CIBERSORT$p = CIBERSORT(sig_matrix = bulkData$indata$C, mixture_file = bulkData$indata$T, QN = FALSE)
result$CIBERSORT$p = t(result$CIBERSORT$p[,1:(ncol(result$CIBERSORT$p)-3)])
result$CIBERSORT$result = getPearsonRMSE(result$CIBERSORT$p,bulkData$indata$P)
result$CIBERSORT$result


#----------------"DeconRNASeq------
#nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
#datasets and reference matrix: signatures, need to be non-negative. 
#"use.scale": whether the data should be centered or scaled, default = TRUE
unloadNamespace("Seurat") #needed for PCA step
library(pcaMethods) #needed for DeconRNASeq to work
result$deconRNASeq$p = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(bulkData$indata$T), 
                                                  signatures = as.data.frame(bulkData$indata$C), proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = FALSE, fig = FALSE)$out.all)
colnames(result$deconRNASeq$p) = colnames(bulkData$indata$T)
require(Seurat)
result$deconRNASeq$result = getPearsonRMSE(result$deconRNASeq$p,bulkData$indata$P)
result$deconRNASeq$result


#-----------RLR----------------
require(MASS)
result$RLR$p = do.call(cbind.data.frame,lapply(apply(bulkData$indata$T,2,function(x) MASS::rlm(x ~ as.matrix(bulkData$indata$C), maxit=100)), function(y) y$coefficients[-1]))
result$RLR$p = apply(result$RLR$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$RLR$p = apply(result$RLR$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$RLR$p) <- unlist(lapply(strsplit(rownames(result$RLR$p),")"),function(x) x[2]))
result$RLR$result = getPearsonRMSE(result$RLR$p,bulkData$indata$P)
result$RLR$result


#-----------DCQ------------------
#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default
require(ComICS)
result$DCQ$p = t(ComICS::dcq(reference_data = bulkData$indata$C, mix_data = bulkData$indata$T, marker_set = as.data.frame(row.names(bulkData$indata$C)) , alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 10)$average)
result$DCQ$p = apply(result$DCQ$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DCQ$p = apply(result$DCQ$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DCQ$result = getPearsonRMSE(result$DCQ$p,bulkData$indata$P)
result$DCQ$result


#-----------elastic_net-----------
#standardize = TRUE by default. lambda=NULL by default 
require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
result$elastic_net$p = apply(bulkData$indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(bulkData$indata$C), y = z, alpha = 0.2, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(bulkData$indata$C), z)$lambda.1se))[1:ncol(bulkData$indata$C)+1,])
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) x/sum(x)) #explicit STO constraint
result$elastic_net$result = getPearsonRMSE(result$elastic_net$p,bulkData$indata$P)
result$elastic_net$result


#----------ridge----------------
# alpha=0
require(glmnet)
result$ridge$p = apply(bulkData$indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(bulkData$indata$C), y = z, alpha = 0, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(bulkData$indata$C), z)$lambda.1se))[1:ncol(bulkData$indata$C)+1,])
result$ridge$p = apply(result$ridge$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$ridge$p = apply(result$ridge$p,2,function(x) x/sum(x)) #explicit STO constraint
result$ridge$result = getPearsonRMSE(result$ridge$p,bulkData$indata$P)
result$ridge$result


#----------lasso-----------------
#alpha=1; shrinking some coefficients to 0. 
require(glmnet)
result$lasso$p = apply(bulkData$indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(bulkData$indata$C), y = z, alpha = 1, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(bulkData$indata$C), z)$lambda.1se))[1:ncol(bulkData$indata$C)+1,])
result$lasso$p = apply(result$lasso$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$lasso$p = apply(result$lasso$p,2,function(x) x/sum(x)) #explicit STO constraint
result$lasso$p[which(is.na(result$lasso$p) == TRUE)] <- 0 #Needed for models where glmnet drops all terms of a model and fit an intercept-only model (very unlikely but possible).
length(colsums(result$lasso$p) == 1)
result$lasso$result = getPearsonRMSE(result$lasso$p,bulkData$indata$P)
result$lasso$result

paste(getwd(),"/result/result_",STRING_name,".rds",sep="")
saveRDS(result,paste(getwd(),"/result/result_",STRING_name,".rds",sep=""))

#----------EPIC----------
require(EPIC)
markers = as.character(bulkData$ShiftGene)

C_EPIC <- list()
common_CTs <- colnames(bulkData$indata$C_ref)
C_EPIC[["sigGenes"]] <- rownames(bulkData$indata$C_ref[markers,common_CTs])
C_EPIC[["refProfiles"]] <- as.matrix(bulkData$indata$C_ref[markers,common_CTs])
C_EPIC$data <- bulkData$data
colnames(C_EPIC$data) <- bulkData$full_phenoData$Class
C_EPIC[["refProfiles.var"]] <- Get_C(C_EPIC$data)[[2]][markers,]
colnames(C_EPIC[["refProfiles.var"]]) <- common_CTs

result$EPIC$p <- t(EPIC::EPIC(bulk=as.matrix(bulkData$indata$T), reference=C_EPIC, withOtherCells=TRUE, scaleExprs=FALSE)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices
result$EPIC$p = result$EPIC$p[!rownames(result$EPIC$p) %in% "otherCells",]
result$EPIC$result = getPearsonRMSE(result$EPIC$p,bulkData$indata$P)
result$EPIC$result


source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")


#--------sc methods_————————————------------

# 
bulkData$sc$T <- bulkData$indata$T
bulkData$sc$C <- bulkData$data
# colnames(bulkData$sc$C) <- bulkData$full_phenoData$Class
bulkData$sc$phenoDataC <- NULL
bulkData$sc$phenoDataC$cellID <- bulkData$full_phenoData$GEO_ID
bulkData$sc$phenoDataC$cellType <- bulkData$full_phenoData$Class
bulkData$sc$phenoDataC$sampleID <- paste("sample",1:length(bulkData$full_phenoData$GEO_ID),sep=c(""))
all(colnames(bulkData$sc$C) == bulkData$sc$phenoDataC$cellID)

# bulkData$sc$C <- bulkData$sc$C[intersect(rownames(bulkData$sc$T),rownames(bulkData$sc$C)),]
# bulkData$sc$T <- bulkData$sc$T[intersect(rownames(bulkData$sc$T),rownames(bulkData$sc$C)),]

bulkData$sc$phenoDataC <- as.data.frame(bulkData$sc$phenoDataC)
colnames(bulkData$sc$phenoDataC) = c("cellID","cellType","SubjectName")
rownames(bulkData$sc$phenoDataC) = bulkData$sc$phenoDataC$cellID

require(xbioc)
bulkData$sc$C.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulkData$sc$C)
                                             ,phenoData = Biobase::AnnotatedDataFrame(bulkData$sc$phenoDataC))
bulkData$sc$T.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulkData$sc$T))

#---------MuSiC-----------
require(MuSiC)
source("F:/wangchenqi/CDSC/MuSiC.R")
result$MuSiC$p = t(music_prop_my(bulk.eset = bulkData$sc$T.eset, sc.eset = bulkData$sc$C.eset, 
                                     clusters = 'cellType',
                                     markers = NULL, normalize = FALSE, samples = 'SubjectName', 
                                     verbose = F)$Est.prop.weighted)
result$MuSiC$result = getPearsonRMSE(result$MuSiC$p, bulkData$indata$P)
result$MuSiC$result

#------------Bisque----------
require(BisqueRNA)
result$Bisque$p <- BisqueRNA::ReferenceBasedDecomposition(bulkData$sc$T.eset, bulkData$sc$C.eset, 
                                                             markers=NULL, use.overlap=FALSE)$bulk.props 
#use.overlap is when there's both bulk and scRNA-seq for the same set of samples
result$Bisque$result = getPearsonRMSE(result$Bisque$p,bulkData$indata$P)
result$Bisque$result

#------------deconvSeq----------
# singlecelldata = bulkData$sc$C.eset@assayData$exprs
# celltypes.sc = as.character(bulkData$sc$phenoDataC$cellType) #To avoid "Design matrix not of full rank" when removing 1 CT 
# tissuedata = bulkData$sc$T.eset @assayData$exprs
# 
# design.singlecell = model.matrix(~ -1 + as.factor(celltypes.sc))
# colnames(design.singlecell) = levels(as.factor(celltypes.sc))
# rownames(design.singlecell) = colnames(singlecelldata)
# 
# dge.singlecell = deconvSeq::getdge(singlecelldata, design.singlecell, ncpm.min = 1, nsamp.min = 4)#, method = "bin.loess"
# set.seed(1234)
# b0.singlecell = deconvSeq::getb0.rnaseq(dge.singlecell, design.singlecell, ncpm.min =1, nsamp.min = 4, sigg=NULL)#
# dge.tissue = deconvSeq::getdge(tissuedata, NULL, ncpm.min = 1, nsamp.min = 4)#, method = "bin.loess"
# result$deconvSeq$out <- deconvSeq::getx1.rnaseq(NB0 = "top_fdr",b0.singlecell, dge.tissue)
# 
# result$deconvSeq$p = t(result$deconvSeq$out$x1) #genes with adjusted p-values <0.05 after FDR correction
# result$deconvSeq$result = getPearsonRMSE(result$deconvSeq$p,bulkData$indata$P)
# result$deconvSeq$result

#--------SCDC-----
library(SCDC)
source("F:/wangchenqi/CDSC/SCDC.R")
result$SCDC$p <- t(SCDC_prop_my(bulk.eset = bulkData$sc$T.eset, sc.eset = bulkData$sc$C.eset, 
                                ct.varname = "cellType", sample = "SubjectName", 
                                ct.sub = unique(as.character(bulkData$sc$phenoDataC$cellType)), iter.max = 200,weight.basis = F)$prop.est.mvw)
result$SCDC$result = getPearsonRMSE(result$SCDC$p,bulkData$indata$P)
result$SCDC$result
#---------DWLS--------
# require(DWLS)

source('F:/wangchenqi/CDSC/DWLS.R')
getwd() 
# path=paste(getwd(),"/DWLS/results_",STRING_name,sep=""); path
# dir.create(path)
# # saveRDS(bulkData$indata$C_ref,paste(getwd(),"/DWLS/result_",STRING_name,"Sig.RData",sep=""))
# 
# if(! dir.exists(path)){ #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created
#   
#   dir.create(path)
#   Signature <- buildSignatureMatrixMAST(scdata = bulkData$sc$C.eset@assayData$exprs, 
#                                         id = as.character(bulkData$sc$phenoData$cellType), 
#                                         path = path, diff.cutoff = 0.5, pval.cutoff = 0.01)
#  
# } else {#re-load signature and remove CT column + its correspondent markers
#   
#   load(paste(path,"Sig.RData",sep="/"))
#   Signature <- Sig
#   
#   if(!is.null(elem)){#to be able to deal with full C and with removed CT
#     
#     Signature = Signature[,!colnames(Signature) %in% elem]
#     CT_to_read <- dir(path) %>% grep(paste(elem,".*RData",sep=""),.,value=TRUE)
#     load(paste(path,CT_to_read,sep="/"))
#     
#     Signature <- Signature[!rownames(Signature) %in% cluster_lrTest.table$Gene,]
#   }
#   
# }
Signature <- bulkData$indata$C_ref
result$DWLS$p <- apply(bulkData$sc$T,2, function(x){
  b = setNames(x, rownames(bulkData$sc$T))
  tr <- trimData(Signature, b)
  RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
})

rownames(result$DWLS$p) <- as.character(unique(bulkData$sc$phenoDataC$cellType))
result$DWLS$p = apply(result$DWLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DWLS$p = apply(result$DWLS$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DWLS$result = getPearsonRMSE(result$DWLS$p,bulkData$indata$P)
result$DWLS$result

result <- CountAllResults(result)
result$all
saveRDS(result, paste(getwd(),"/result/result_",STRING_name,".rds",sep=""))
#-----complete deconvolution methods————----------
#--------DSA-------------------
require(CellMix)
md = bulkData$markers1_shift
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

result$DSA = CellMix::ged(as.matrix(bulkData$indata$T), ML, method = "DSA", log = TRUE)
result$DSA$c = result$DSA@fit@W
result$DSA$p = result$DSA@fit@H
result$DSA$t = result$DSA$c%*%result$DSA$p

result$DSA$result$c = getPearsonRMSE(result$DSA$c,bulkData$indata$C)
result$DSA$result$p = getPearsonRMSE(result$DSA$p,bulkData$indata$P)
result$DSA$result$t = getPearsonRMSE(result$DSA$t,bulkData$indata$T)

result$DSA$result$all <- cbind(result$DSA$result$p,result$DSA$result$c)
result$DSA$result$all <-  cbind(result$DSA$result$all, result$DSA$result$t)
colnames(result$DSA$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$DSA$result <- result$DSA$result$all
result$DSA$result

#---------ssKL------------------
# require(CellMix)
# md = bulkData$markers1_shift
# ML = CellMix::MarkerList()
# ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

result$ssKL <- CellMix::ged(as.matrix(bulkData$indata$T), ML, method = "ssKL", sscale = FALSE, maxIter=500, log = FALSE)
result$ssKL$c = result$ssKL@fit@W
result$ssKL$p = result$ssKL@fit@H
result$ssKL$t = result$ssKL$c%*%result$ssKL$p

result$ssKL$result$c = getPearsonRMSE(result$ssKL$c,bulkData$indata$C)
result$ssKL$result$p = getPearsonRMSE(result$ssKL$p,bulkData$indata$P)
result$ssKL$result$t = getPearsonRMSE(result$ssKL$t,bulkData$indata$T)

result$ssKL$result$all <- cbind(result$ssKL$result$p,result$ssKL$result$c)
result$ssKL$result$all <-  cbind(result$ssKL$result$all,result$ssKL$result$t)
colnames(result$ssKL$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssKL$result <- result$ssKL$result$all
result$ssKL$result

#----------ssFrobenius-----------------
# require(CellMix)
# md = bulkData$markers1_shift
# ML = CellMix::MarkerList()
# ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

result$ssFrobenius <- CellMix::ged(as.matrix(bulkData$indata$T), ML, method = "ssFrobenius", sscale = TRUE, maxIter = 500, log = FALSE) #equivalent to coef(CellMix::ged(T,...)
result$ssFrobenius$c = result$ssFrobenius@fit@W
result$ssFrobenius$p = result$ssFrobenius@fit@H
result$ssFrobenius$t = result$ssFrobenius$c%*%result$ssFrobenius$p

result$ssFrobenius$result$c = getPearsonRMSE(result$ssFrobenius$c,bulkData$indata$C)
result$ssFrobenius$result$p = getPearsonRMSE(result$ssFrobenius$p,bulkData$indata$P)
result$ssFrobenius$result$t = getPearsonRMSE(result$ssFrobenius$t,bulkData$indata$T)

result$ssFrobenius$result$all <- cbind(result$ssFrobenius$result$p,result$ssFrobenius$result$c)
result$ssFrobenius$result$all <-  cbind(result$ssFrobenius$result$all,result$ssFrobenius$result$t)
colnames(result$ssFrobenius$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssFrobenius$result <- result$ssFrobenius$result$all
result$ssFrobenius$result

#-----------deconf------------
require(CellMix)
md = bulkData$markers1_shift
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

all(rownames(bulkData$indata$T)==bulkData$markers1_shift$gene)
all(rownames(bulkData$indata$T)==rownames(bulkData$shif_marker))

result$deconf <- CellMix::ged(as.matrix(bulkData$indata$T), ncol(bulkData$indata$C_ref))
result$deconf$c = result$deconf@fit@W
result$deconf$p = result$deconf@fit@H
result$deconf$t = result$deconf$c%*%result$deconf$p
ctlabels <- Row_label(result$deconf$c,bulkData$indata$C_ref,leastnum=3);ctlabels
colnames(result$deconf$c) <- ctlabels
rownames(result$deconf$p) <- ctlabels

# result$deconf <- CellMix::ged(as.matrix(bulkData$indata$T), ML, method = "deconf")
# result$deconf$c = result$deconf@fit@W
# result$deconf$p = result$deconf@fit@H
# result$deconf$t = result$deconf$c%*%result$deconf$p

result$deconf$result$c = getPearsonRMSE(result$deconf$c,bulkData$indata$C)
result$deconf$result$p = getPearsonRMSE(result$deconf$p,bulkData$indata$P)
result$deconf$result$t = getPearsonRMSE(result$deconf$t,bulkData$indata$T)

result$deconf$result$all <- cbind(result$deconf$result$p,result$deconf$result$c)
result$deconf$result$all <-  cbind(result$deconf$result$all,result$deconf$result$t)
colnames(result$deconf$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$deconf$result <- result$deconf$result$all
result$deconf$result

#--------CDSeq-------------------
# library(CDSeq)
# star_time <- Sys.time()
# result$CDSeq.noRef <- CDSeq(bulk_data =  as.matrix(bulkData$indata$T), 
#                             cell_type_number = dim(bulkData$indata$C)[2], 
#                             mcmc_iterations = 100, 
#                             cpu_number=20)
# end_time <- Sys.time()
# end_time-star_time
# result$CDSeq.noRef$result <- calculate_result(result$CDSeq.noRef$estGEP,result$CDSeq.noRef$estProp,
#                                               bulkData$indata$T, bulkData$indata$C, bulkData$indata$C, bulkData$indata$P)
# result$CDSeq.noRef$time <- end_time-star_time
# result$CDSeq.noRef$result
# 
# star_time <- Sys.time()
# result$CDSeq.haveRef<-CDSeq(bulk_data =  as.matrix(bulkData$indata$T), 
#                             cell_type_number = dim(bulkData$indata$C)[2], 
#                             mcmc_iterations = 20, 
#                             # gene_length = as.vector(gene_length), 
#                             reference_gep = as.matrix(bulkData$indata$C),  # gene expression profile of pure cell lines
#                             cpu_number = 8)
# end_time <- Sys.time()
# end_time-star_time
# result$CDSeq.haveRef$result <- calculate_result(result$CDSeq.haveRef$estGEP,result$CDSeq.haveRef$estProp,
#                                                 bulkData$indata$T, bulkData$indata$C, bulkData$indata$C, bulkData$indata$P)
# result$CDSeq.haveRef$time <- end_time-star_time
# result$CDSeq.haveRef$result
#------TOAST + NMF-------
# TOAST方法宝函数只能输入基因大于样本的数据。。
require(DeCompress)
dim(bulkData$indata$T)
result$TOAST$toast.nmf <- DeCompress::csDeCompress(Y_raw = bulkData$indata$T,
                                                   K = dim(bulkData$indata$C)[2],
                                                   nMarker = nrow(bulkData$indata$T),
                                                   FUN = nmfOut,
                                                   TotalIter = 30)
require(NMF)
result$TOAST$fin.nmf = nmf(x = bulkData$indata$T,
                           rank = dim(bulkData$indata$C)[2])
result$TOAST$ppp = t(coef(result$TOAST$fin.nmf))
result$TOAST$ppp = t(apply(result$TOAST$ppp,1,function(c) c/sum(c)))
result$TOAST$p <- t(result$TOAST$ppp)
result$TOAST$c <- basis(result$TOAST$fin.nmf)
result$TOAST$t <- result$TOAST$c%*%result$TOAST$p
nmf.res = list(prop = result$TOAST$ppp,
               sig = basis(result$TOAST$fin.nmf))

rownames(bulkData$indata$P)
labels <- Row_label(t(result$TOAST$p),t(bulkData$indata$P));labels
rownames(result$TOAST$p) <- labels
colnames(result$TOAST$c) <- labels
result$TOAST$p <- result$TOAST$p[rownames(bulkData$indata$P),]
result$TOAST$c <- result$TOAST$c[,colnames(bulkData$indata$C)]

result$TOAST$result$c = getPearsonRMSE(result$TOAST$c,bulkData$indata$C)
result$TOAST$result$p = getPearsonRMSE(result$TOAST$p,bulkData$indata$P)
result$TOAST$result$t = getPearsonRMSE(result$TOAST$t,bulkData$indata$T)

result$TOAST$result$all <- cbind(result$TOAST$result$p,result$TOAST$result$c)
result$TOAST$result$all <-  cbind(result$TOAST$result$all,result$TOAST$result$t)
colnames(result$TOAST$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$TOAST$result <- result$TOAST$result$all
result$TOAST$result

#-------Linseed---------------
require(DeCompress)
result$Linseed$Linseed.rs = DeCompress::linCor(yref = bulkData$indata$T,
                                               iters = 100,
                                               pval = 100,
                                               n.types = dim(bulkData$indata$C)[2],
                                               scree = 'drop',
                                               logTransform = F)

names(result$Linseed$Linseed.rs) = names(result$Linseed$nmf.res)
result$Linseed$c <- result$Linseed$Linseed.rs[[2]]
result$Linseed$p <- result$Linseed$Linseed.rs[[1]]
result$Linseed$t <- result$Linseed$c%*%result$Linseed$p
rownames(bulkData$indata$P)
labels <- Row_label(t(result$Linseed$p),t(bulkData$indata$P));labels
rownames(result$Linseed$p) <- labels
colnames(result$Linseed$c) <- labels
result$Linseed$p <- result$Linseed$p[rownames(bulkData$indata$P),]
result$Linseed$c <- result$Linseed$c[,colnames(bulkData$indata$C)]

result$Linseed$result$c = getPearsonRMSE(result$Linseed$c,bulkData$indata$C)
result$Linseed$result$p = getPearsonRMSE(result$Linseed$p,bulkData$indata$P)
result$Linseed$result$t = getPearsonRMSE(result$Linseed$t,bulkData$indata$T)

result$Linseed$result$all <- cbind(result$Linseed$result$p,result$Linseed$result$c)
result$Linseed$result$all <-  cbind(result$Linseed$result$all,result$Linseed$result$t)
colnames(result$Linseed$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$Linseed$result <- result$Linseed$result$all
result$Linseed$result

#-----CellDistinguisher-----
require(DeCompress)
result$CellDistinguisher <- CellDistinguisher::gecd_CellDistinguisher(
  bulkData$indata$T,
  genesymb=rownames(bulkData$indata$T),
  numCellClasses = dim(bulkData$indata$C)[2],
  minDistinguisherAlternatives=1,
  maxDistinguisherAlternatives=100,
  minAlternativesLengthsNormalized=0.5,
  expressionQuantileForFilter=0.999,
  expressionConcentrationRatio=0.333,
  verbose=0)
result$CellDist.deconv <-
  tryCatch(CellDistinguisher::gecd_DeconvolutionByDistinguishers(
    as.matrix(bulkData$indata$T),
    result$CellDistinguisher$bestDistinguishers,
    nonNegativeOnly = FALSE,
    convexSolution = FALSE,
    verbose = 0),
    error = function(e) return(list(sampleCompositions =
                                      matrix(rep(1/dim(bulkData$indata$C)[2],
                                                 dim(bulkData$indata$C)[2]*ncol(bulkData$indata$T)),
                                             ncol=dim(bulkData$indata$C)[2]))))
if (length(result$CellDist.deconv) > 1){
  result$CellDistinguisher$p <- result$CellDist.deconv$sampleComposition
  result$CellDistinguisher$c <- result$CellDist.deconv$cellSubclassSignatures
}else {
  result$CellDistinguisher$p <- nmf.res$prop
  result$CellDistinguisher$c <- nmf.res$sig
}
result$CellDistinguisher$t <- result$CellDistinguisher$c%*%result$CellDistinguisher$p
rownames(bulkData$indata$P)
labels <- Row_label(t(result$CellDistinguisher$p),t(bulkData$indata$P))
rownames(result$CellDistinguisher$p) <- labels
colnames(result$CellDistinguisher$c) <- labels
result$CellDistinguisher$p <- result$CellDistinguisher$p[rownames(bulkData$indata$P),]
result$CellDistinguisher$c <- result$CellDistinguisher$c[,colnames(bulkData$indata$C)]
result$CellDistinguisher$result$c = getPearsonRMSE(result$CellDistinguisher$c,bulkData$indata$C)
result$CellDistinguisher$result$p = getPearsonRMSE(result$CellDistinguisher$p,bulkData$indata$P)
result$CellDistinguisher$result$t = getPearsonRMSE(result$CellDistinguisher$t,bulkData$indata$T)
result$CellDistinguisher$result$all <- cbind(result$CellDistinguisher$result$p,result$CellDistinguisher$result$c)
result$CellDistinguisher$result$all <-  cbind(result$CellDistinguisher$result$all,result$CellDistinguisher$result$t)
colnames(result$CellDistinguisher$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$CellDistinguisher$result <- result$CellDistinguisher$result$all
result$CellDistinguisher$result

#-------------DeCompress--------
# result$DeCompress$decompress.res = bestDeconvolution( bulkData$indata$T,
#                                                       n.types = dim(bulkData$indata$C)[2],
#                                                       scree = 'cumvar',
#                                                       logTransform = F,
#                                                       known.props = NULL,
#                                                       methods = c('TOAST',
#                                                                   'Linseed',
#                                                                   'CellDistinguisher'))

#--------------all————————-----------------
result <- CountAllResults(result,MyMethodName)
result$all
STRING_name
paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep="")
saveRDS(result,paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))

#----------find Paramater ///----------
lambda1 <- c(0,10^-4,10^-3,10^-2,10^-1,1,10,100)
lambda2 <- c(0,10^-4,10^-3,10^-2,10^-1,1,10,100)
lambdaC <- c(0,10^-2,10^-1,10^0,10^1,10^2,10^3,10^4,10^5)
bulkData$seedd = 44
bulkData$indata$TerCondition = 10^-8

#------------cyclic to do deconvolution----------
pb <- txtProgressBar(style = 3)
star_time <- Sys.time()
bulkData$Ss <- SM(t(bulkData$indata$T))
bulkData$Sg <- SM(bulkData$indata$T)
result_para_c = NULL
result_para_p = NULL
pearson_para_c = NULL
number_iter = NULL
num = 1

para_lambda_44_8 = NULL
library(dplyr)
for (dir_i in 1:length(lambda1)){
  for (dir_j in 1:length(lambda2)){
    for (dir_k in 1:length(lambdaC)){
      result_CDSC = CDSC_3(bulkData$indata$T, bulkData$indata$C, dim(bulkData$indata$C)[2], 
                           lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k],
                           bulkData$indata$TerCondition,bulkData$seedd,
                           bulkData$Ss,bulkData$Sg,all_number=3000)
      result_para_c[[num]] <- result_CDSC[[1]]
      result_para_p[[num]] <- result_CDSC[[2]]
      number_iter[num] <- result_CDSC[[3]]
      
      pearson_para_c[[num]] <- cor(result_para_c[[num]],bulkData$indata$C)
      result1 = NULL
      if(!all(is.na(pearson_para_c[[num]]) == FALSE) ){
        break
      }
      result1 <- calculate_result(result_para_c[[num]],result_para_p[[num]],
                                  bulkData$indata$T,bulkData$indata$C,bulkData$indata$C_ref,bulkData$indata$P,
                                  lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k],number_iter[num],
                                  bulkData$seedd,bulkData$indata$TerCondition)
      # result
      # rownames(result) = 1:length(rownames(result))
      para_lambda_44_8 =  rbind(para_lambda_44_8,result1)
      # cat("\nlambda1=",lambda1[dir_i],",lambda2=",lambda2[dir_j],
      #     ",lambdac=",lambdaC[dir_k],",迭代次数",number_iter[num])
      num = num+1
      setTxtProgressBar(pb, num/(length(lambda1)*length(lambda2)*length(lambdaC)))
    }
  }
}
end_time <- Sys.time()
close(pb)
print(star_time)
print(end_time)
end_time-star_time

STRING_name
saveRDS(para_lambda_44_8,paste(getwd(),"/paramater/para_lambda_",STRING_name,"_123.rds",sep=""))


#-------------hot plot---------
rm (list=ls ())
setwd("F:/wangchenqi/CDSC/3realBulk")
source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/3realBulk/bulkData"));a
a <- c("BreastBlood","CellLines","LiverBrainLung","RatBrain","Retina" )
MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,"deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",
                  "MuSiC","Bisque","SCDC", "DWLS",
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
bulkData = NULL
result = NULL 
for(i in 1:length(a)){
  STRING_name = a[i]; STRING_name
  # bulkData[[i]] <- readRDS(paste(getwd(),"/bulkData/",STRING_name,".rds",sep=""))
  result[[i]] <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))
}
map = NULL
for(i in 1:length(a)){
  # result[[i]] <- CountAllResults(result[[i]],MyMethodName)
  # map <- cbind(map,result[[i]]$all$RMSE_to_P)
  map <- cbind(map,result[[i]]$all$Peason_to_P)
}
rownames(map) <- MyMethodName;
colnames(map) <- a
map[which(map == 0)] <- NA;map


library(pheatmap)#c("navy", "white", "firebrick3") color
plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 42, cellheight = 12,gaps_row = c(12,16)
                         ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.3f"
                         ,border_color = "black",color = colorRampPalette(c( "white", "firebrick3"))(30)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "Pearson of P")

plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 42, cellheight = 12,gaps_row = c(12,16)
                         ,display_numbers = TRUE, number_format = "%.3f"
                         ,border_color = "black",color = colorRampPalette(c( "white", "navy"))(70)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "RMSE of P")

plot_pheatmp <- pheatmap(map[c(1,17:24),],cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 42, cellheight = 12#,gaps_row = c(11,16)
                         ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.3f"
                         ,border_color = "black",color = colorRampPalette(c("white", "firebrick3"))(30)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "Pearson of C")

plot_pheatmp <- pheatmap(map[c(1,17:24),],cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 42, cellheight = 12#,gaps_row = c(10)
                         ,display_numbers = TRUE, number_format = "%.1f"
                         ,border_color = "black",color = colorRampPalette(c( "white", "navy"))(70)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "RMSE of C")
apply(map, 1, mean)
plot_pheatmp;
library(ggplot2)
eoffice::topptx(plot_pheatmp,filename = 
                    "F:/wangchenqi/CDSC/pictures/class3_Pearson_P.pptx")

# 泛化能力的箱线图--p----
mapBox <- NULL
mapBox =  data.frame(pearson = map[1,], method = rownames(map)[1],row.names = NULL)
nn1 <- nrow(mapBox)
for (i in 1:length(rownames(map))) {
  mapBox <- rbind(mapBox, data.frame(pearson = map[i,], method = rownames(map)[i],row.names = NULL))
}

mapBox <- mapBox[-(1:nn1),]
library(ggplot2)
plot_class1 =  ggplot(mapBox, aes(factor(method,levels=MyMethodName),pearson,fill=method)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson of P")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_class1;
eoffice::topptx(plot_class1,filename = 
                  "F:/wangchenqi/CDSC/pictures/class3_fanhua_p.pptx")

# 泛化能力的箱线图----c-----
mapBox <- NULL
methodsNames2 <-  c("CDSC3","CDSC2", "DSA","ssKL" ,"ssFrobenius","deconf",
                    "TOAST","Linseed","CellDistinguisher")
map2 <- map[methodsNames2,]
mapBox =  data.frame(pearson = map2[1,], method = rownames(map2)[1],row.names = NULL)
nn1 <- nrow(mapBox)
for (i in 1:length(rownames(map2))) {
  mapBox <- rbind(mapBox, data.frame(pearson = map2[i,], method = rownames(map2)[i],row.names = NULL))
}

mapBox <- mapBox[-(1:nn1),]
library(ggplot2)
plot_class1 =  ggplot(mapBox, aes(factor(method,levels=methodsNames2),pearson,fill=method)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson of C")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_class1;
eoffice::topptx(plot_class1,filename = 
                  "F:/wangchenqi/CDSC/pictures/class3_fanhua_c.pptx")

#------柱状图--p------
rm (list=ls ())
setwd("F:/wangchenqi/CDSC/3realBulk")

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/3realBulk/bulkData"));a
a <- c("CellLines","LiverBrainLung","RatBrain" );a

# result$all

bulkData = NULL
result = NULL 

for(i in 1:length(a)){
  STRING_name = a[i]; STRING_name
  bulkData[[i]] <- readRDS(paste(getwd(),"/bulkData/",STRING_name,".rds",sep=""))
  result[[i]] <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))
}
rownames(result[[i]]$CDSC3$dec$p)

for (i in 1:length(a)) {
  ctlabels <- Row_label(result[[i]]$CDSC3$dec$c,bulkData[[i]]$indata$C,leastnum=3)
  rownames(result[[i]]$CDSC3$dec$p) <- ctlabels
  colnames(result[[i]]$CDSC3$dec$c) <- ctlabels
  ctlabels_ <- rownames(bulkData[[i]]$indata$P)
  result[[i]]$CDSC3$p <- result[[i]]$CDSC3$dec$p[ctlabels_, ]
  result[[i]]$CDSC3$c <- result[[i]]$CDSC3$dec$c[ ,ctlabels_]
}
# for Retina

map = NULL
for(i in 1:length(a)){
  
  map <- rbind(map,data.frame(group = a[i],
                              pearson = diag(cor(result[[i]]$CDSC3$p,bulkData[[i]]$indata$P))))
}

map = NULL
for(i in 1:length(a)){
  map <- rbind(map,data.frame(group = a[i],
                              pearson = diag(cor(result[[i]]$CDSC3$p,bulkData[[i]]$indata$P))))
}
map[is.na(map)] <-1
library(dplyr)
map_mean <- map %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarize(
    count = n(),
    mean = mean(pearson),
    sd = sd(pearson)
  )
library(ggplot2)
p1 <- ggplot()+ 
  geom_bar(data=map_mean,mapping=aes(x=factor(group,levels = a),y=mean,fill=group), # fill填充
           position="dodge", # 柱状图格式
           stat="identity", # 数据格式
           width = 0.7,
           show.legend = F)+  # 柱状图尺寸
  scale_fill_manual(values = c("#80AFBF","#80AFBF","#80AFBF","#80AFBF","#80AFBF","#80AFBF"))+ # 柱状图颜色, "#DA635D","#B1938B"
  # geom_signif(data=plot_data2,mapping=aes(x=group,y=SOD), # 不同组别的显著性
  #             comparisons = list(c("C", "HT"), # 哪些组进行比较
  #                                c("HI", "HT")),
  #             annotation=c("**"), # 显著性差异做标记
  #             map_signif_level=T, # T为显著性，F为p value
  #             tip_length=c(0.04,0.04,0.05,0.05), # 修改显著性那个线的长短
  #             y_position = c(4100,3000), # 设置显著性线的位置高度
  #             size=1, # 修改线的粗细
  #             textsize = 10, # 修改*标记的大小
  #             test = "t.test")+ # 检验的类型
  geom_errorbar(data=map_mean,mapping=aes(x = group,ymin = mean-sd, ymax = mean+sd), # 误差线添加
                width = 0.1, #误差线的宽度
                color = 'black', #颜色
                size=0.8)+ #粗细
  scale_y_continuous(limits =c(0, 1.1) ,expand = c(0,0))+ # y轴的范围
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  
  labs(title="",x="Dataset",y="Pearson of samples in P")+ # 添加标题，x轴，y轴内容
  theme(plot.title = element_text(size = 20,
                                  colour = "red",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", # 修改图例的标题
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 修改X轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 1, # 位置
                                   hjust = 1, 
                                   angle = 45), #角度
        axis.text.y = element_text(size = 13,  # 修改y轴上字体大小，
                                   # family = "myFont", # 类型
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0) #角度
  ) 
# emf(file = "SOD.emf") # 打开一个矢量图画布，这种格式的图片放在word里不会失真
print(p1) # 打印图片
eoffice::topptx(p1,filename = "F:/wangchenqi/CDSC/pictures/class3_bar_samples+.pptx")

# ---- sample情况的散点图----
library(ggplot2)
plot_class2 =  ggplot(map, aes(factor(group,levels=a),pearson,fill=group)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_class2
eoffice::topptx(plot_class2,filename = "F:/wangchenqi/CDSC/pictures/class3_boxplot_samples.pptx")

