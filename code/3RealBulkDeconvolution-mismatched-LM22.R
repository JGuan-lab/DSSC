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
# Melanoma <- list(
#   T = read.table("F:/wangchenqi/CDSC/3realBulk/data/Melanoma_bulk.txt",row.names = 1,header = T),
#   C = read.table("F:/wangchenqi/CDSC/3realBulk/data/Melanoma_signature.txt",row.names = 1,header = T,sep = "\t"),
#   P = read.table("F:/wangchenqi/CDSC/3realBulk/data/Melanoma_groundtruth.txt",row.names = 1,header = T,sep = "\t")
# )
WholeBlood <- list(
  T = read.table("F:/wangchenqi/CDSC/3realBulk/data/WholeBlood_bulk.txt",row.names = 1,header = T),
  # C = read.table("F:/wangchenqi/CDSC/3realBulk/data/LM22_source_GEPs.txt",row.names = 1,header = T,sep = "\t"),
  P = read.table("F:/wangchenqi/CDSC/3realBulk/data/WholeBlood_groundtruth.txt",row.names = 1,header = T,sep = "\t")
)
#------simulate--------
source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")


#----------WholeBlood-------
# this is a complete data form CIBERSORT, which contains T,C,P
# LM22 signature
# WholeBlood$C$nsclc <- read.table("F:/wangchenqi/CDSC/3realBulk/data/NSCLC_PBMCs_scRNAseq_sigmatrix.txt",row.names = 1,header = T,sep = "\t")
# WholeBlood$C$nsclc_sc <- read.table("F:/wangchenqi/CDSC/3realBulk/data/NSCLC_PBMCs_scRNAseq_refsample.txt",row.names = 1,header = T,sep = "\t")
# WholeBlood$C$FL <- read.table("F:/wangchenqi/CDSC/3realBulk/data/FL_scRNAseq_matrix.txt",row.names = 1,header = T,sep = "\t")
colnames(WholeBlood$P)
# "Neutrophils" 
# "Lymphocytes" 
# "Monocytes"   
# "T.cells"     
# "T.cells.CD8" 
# "T.cells.CD4" 
# "B.cells"     
# "NK.cells"
# 8 kinds of celltypes

# --LM22-----------
WholeBlood$C$LM22 <- read.table("F:/wangchenqi/CDSC/3realBulk/data/LM22.txt",row.names = 1,header = T,sep = "\t")
library(readxl)
WholeBlood$C$LM22_gene_CT <- read_excel("F:/wangchenqi/CDSC/3realBulk/data/LM22_geneToCT.xlsx",col_names = T)

colnames(WholeBlood$C$LM22)
# [1] "B.cells.naive"                "B.cells.memory"               "Plasma.cells"                 "T.cells.CD8"
# [5] "T.cells.CD4.naive"            "T.cells.CD4.memory.resting"   "T.cells.CD4.memory.activated" "T.cells.follicular.helper"
# [9] "Tregs"                        "T.cells.gamma.delta"          "NK.cells.resting"             "NK.cells.activated"
# [13] "Monocytes"                    "Macrophages.M0"               "Macrophages.M1"               "Macrophages.M2"
# [17] "Dendritic.cells.resting"      "Dendritic.cells.activated"    "Mast.cells.resting"           "Mast.cells.activated"
# [21] "Eosinophils"                  "Neutrophils"
WholeBlood$C$LM22_gene_CT <- as.data.frame(WholeBlood$C$LM22_gene_CT)
rownames(WholeBlood$C$LM22_gene_CT) <- WholeBlood$C$LM22_gene_CT[,1]
WholeBlood$C$LM22_gene_CT <- WholeBlood$C$LM22_gene_CT[,-1]
colnames(WholeBlood$C$LM22_gene_CT) <- colnames(WholeBlood$C$LM22)

WholeBlood$sig$lm22 <- WholeBlood$C$LM22[,c("B.cells.naive", "B.cells.memory",
                                            "T.cells.CD8",
                                            "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
                                            "NK.cells.resting", "NK.cells.activated",
                                            "Monocytes")]
WholeBlood$sig$lm22_gene_CT <- WholeBlood$C$LM22_gene_CT[,c("B.cells.naive", "B.cells.memory",
                                                            "T.cells.CD8",
                                                            "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated",
                                                            "NK.cells.resting", "NK.cells.activated",
                                                            "Monocytes")]
# colnames(WholeBlood$sig$lm22) <- c("B.cells","T.cells.CD8","T.cells.CD4","NK.cells","Monocytes")
# colnames(WholeBlood$sig$lm22_gene_CT) <- c("B.cells","T.cells.CD8","T.cells.CD4","NK.cells","Monocytes")

WholeBlood$sig$lm22_gene_CT <- WholeBlood$sig$lm22_gene_CT[which(rowSums(WholeBlood$sig$lm22_gene_CT)>0),]
celltypeName = colnames(WholeBlood$sig$lm22)

gene_ct <- WholeBlood$sig$lm22_gene_CT;genes <- rownames(gene_ct)
ct<-list();RightGene <- NULL
length(genes)
for(i in celltypeName){
  ct[[i]] <- genes[which(gene_ct[[i]] ==1)]
  genes <- genes[-which(gene_ct[[i]] ==1)]
  gene_ct = gene_ct[-which(gene_ct[[i]] ==1),]
  RightGene = c(RightGene,ct[[i]])
}
length(genes);length(RightGene)

WholeBlood$sig$lm22_gene_CT <- WholeBlood$sig$lm22_gene_CT[RightGene,]
WholeBlood$sig$lm22 <- WholeBlood$sig$lm22[RightGene,]
# #
colnames(WholeBlood$sig$lm22)
# interaction_ct = intersect(colnames(WholeBlood$sig$lm22),colnames(WholeBlood$P))
# interaction_ct
# WholeBlood$sig$lm22 <- WholeBlood$sig$lm22[,interaction_ct]
# WholeBlood$GroundTruth$lm22 <- WholeBlood$P[,interaction_ct]

WholeBlood$indata$C <- WholeBlood$sig$lm22
WholeBlood$indata$P <- t(WholeBlood$P)
WholeBlood$indata$T <- WholeBlood$T

gene <- intersect(rownames(WholeBlood$indata$C),rownames(WholeBlood$T))
length(gene)
WholeBlood$indata$T <- WholeBlood$indata$T[gene,]
WholeBlood$indata$C <- WholeBlood$indata$C[gene,]
# #
markers <- WholeBlood$sig$lm22_gene_CT
Alist <- list();markerslist = NULL
for (i in celltypeName) {
  markers[ct[[i]],][[i]] <- i
  Alist[[i]] <- markers[ct[[i]],][[i]]
  markerslist <- c(markerslist,Alist[[i]])
}

markerslist=as.data.frame(markerslist);rownames(markerslist) <- rownames(WholeBlood$sig$lm22_gene_CT)
WholeBlood$scDataSimulate$markerslist <- markerslist
WholeBlood$scDataSimulate$markerslist$gene <- rownames(WholeBlood$scDataSimulate$markerslist)
WholeBlood$scDataSimulate$markerslist$CT <- WholeBlood$scDataSimulate$markerslist$markerslist

WholeBlood$indata$T <- as.matrix(WholeBlood$indata$T)
WholeBlood$indata$C <- as.matrix(WholeBlood$indata$C)
WholeBlood$indata$P <- as.matrix(SumEqual_1(WholeBlood$indata$P))
WholeBlood$indata$C_ref <- WholeBlood$indata$C

WholeBlood$T ['BLK',]
WholeBlood$indata$T ['BLK',]
# md = bulkData$scDataSimulate$markerslist[rownames(bulkData$indata$T),]
# ML = CellMix::MarkerList()
# ML@.Data <- tapply(as.character(rownames(md)),as.character(md$CT),list)
saveRDS(WholeBlood,"F:/wangchenqi/CDSC/3realBulk/bulkdata/WholeBlood_lm22.rds")

### --3'PBMCs-----------
WholeBlood$C$PBMCs_3 <- read.table("F:/wangchenqi/CDSC/3realBulk/data/3PBMCs_scRNAseq_matrix.txt",row.names = 1,header = T,sep = "\t")

cellnames <- colnames(WholeBlood$C$PBMCs_3)
WholeBlood$C$PBMCs_3_pData = NULL
WholeBlood$C$PBMCs_3_pData = cellnames
WholeBlood$C$PBMCs_3_pData <- cbind(WholeBlood$C$PBMCs_3_pData, gsub("\\.[0-9]*$","",cellnames))
WholeBlood$C$PBMCs_3_pData <- cbind(WholeBlood$C$PBMCs_3_pData, 1)
colnames(WholeBlood$C$PBMCs_3_pData) <- c("cellID",	"cellType",	"sampleID")

WholeBlood$C$PBMCs_3 <- as.matrix(WholeBlood$C$PBMCs_3)
WholeBlood$C$PBMCs_3_pData <- as.data.frame(WholeBlood$C$PBMCs_3_pData)
rownames(WholeBlood$C$PBMCs_3_pData) <- WholeBlood$C$PBMCs_3_pData$cellID
# find marker genes
scData <- list(data = WholeBlood$C$PBMCs_3, full_phenoData = WholeBlood$C$PBMCs_3_pData)
WholeBlood$scDataSimulate <- scSimulateC(scData,
                                        leastNum = 0,
                                        plotmarker = F,
                                        norm1 = "none",
                                        log2.threshold=log2(2))
table(WholeBlood$scDataSimulate$markerslist$CT)
nrow(WholeBlood$scDataSimulate$markerslist)

# shift CT
CellType <- intersect(colnames(WholeBlood$scDataSimulate$C),colnames(WholeBlood$P));CellType
WholeBlood$sig$PBMCs_3 <- WholeBlood$scDataSimulate$C[,CellType]
WholeBlood$GroundTruth$PBMCs_3 <- WholeBlood$P[,CellType]
WholeBlood$scDataSimulate$markerslist <- WholeBlood$scDataSimulate$markerslist[
  which(WholeBlood$scDataSimulate$markerslist$CT %in% CellType),]
# shift genes
WholeBlood$indata$C <- WholeBlood$sig$PBMCs_3
WholeBlood$indata$P <- t(WholeBlood$GroundTruth$PBMCs_3)
WholeBlood$indata$T <- WholeBlood$T

gene <- intersect(rownames(WholeBlood$indata$C),rownames(WholeBlood$T));gene <- intersect(gene,WholeBlood$scDataSimulate$markerslist$gene);length(gene)
WholeBlood$indata$T <- WholeBlood$indata$T[gene,]
WholeBlood$indata$T <- WholeBlood$indata$T[which(rowSums(WholeBlood$indata$T) > 0),]
length(which(rowSums(WholeBlood$indata$T) > 0))
gene <- gene[which(rowSums(WholeBlood$indata$T) > 0)];length(gene)
WholeBlood$indata$C <- WholeBlood$indata$C[gene,]

WholeBlood$indata$T <- as.matrix(WholeBlood$indata$T)
WholeBlood$indata$C <- as.matrix(WholeBlood$indata$C)
WholeBlood$indata$P <- as.matrix(SumEqual_1(WholeBlood$indata$P))
WholeBlood$indata$C_ref <- WholeBlood$indata$C

WholeBlood$T ['BLK',]
WholeBlood$indata$T ['BLK',]
saveRDS(WholeBlood,"F:/wangchenqi/CDSC/3realBulk/bulkdata/WholeBlood_3pbmcs.rds")

### --5'PBMCs-------------
WholeBlood$C$PBMCs_5 <- read.table("F:/wangchenqi/CDSC/3realBulk/data/5PBMCs_scRNAseq_matrix.txt",row.names = 1,header = T,sep = "\t")

cellnames <- colnames(WholeBlood$C$PBMCs_5)
WholeBlood$C$PBMCs_5_pData = NULL
WholeBlood$C$PBMCs_5_pData = cellnames
WholeBlood$C$PBMCs_5_pData <- cbind(WholeBlood$C$PBMCs_5_pData, gsub("\\.[0-9]*$","",cellnames))
WholeBlood$C$PBMCs_5_pData <- cbind(WholeBlood$C$PBMCs_5_pData, 1)
colnames(WholeBlood$C$PBMCs_5_pData) <- c("cellID",	"cellType",	"sampleID")

WholeBlood$C$PBMCs_5 <- as.matrix(WholeBlood$C$PBMCs_5)
WholeBlood$C$PBMCs_5_pData <- as.data.frame(WholeBlood$C$PBMCs_5_pData)
rownames(WholeBlood$C$PBMCs_5_pData) <- WholeBlood$C$PBMCs_5_pData$cellID
# find marker genes
scData <- list(data = WholeBlood$C$PBMCs_5, 
               full_phenoData = WholeBlood$C$PBMCs_5_pData)
WholeBlood$scDataSimulate <- scSimulateC(scData,
                                         leastNum = 0,
                                         plotmarker = F,
                                         norm1 = "none",
                                         log2.threshold=log2(2))
table(WholeBlood$scDataSimulate$markerslist$CT)
# shift CT
CellType <- intersect(colnames(WholeBlood$scDataSimulate$C),colnames(WholeBlood$P));CellType
WholeBlood$sig$PBMCs_5 <- WholeBlood$scDataSimulate$C[,CellType]
WholeBlood$GroundTruth$PBMCs_5 <- WholeBlood$P[,CellType]
WholeBlood$scDataSimulate$markerslist <- WholeBlood$scDataSimulate$markerslist[
  which(WholeBlood$scDataSimulate$markerslist$CT %in% CellType),]
# shift genes
WholeBlood$indata$C <- WholeBlood$sig$PBMCs_5
WholeBlood$indata$P <- t(WholeBlood$GroundTruth$PBMCs_5)
WholeBlood$indata$T <- WholeBlood$T

gene <- intersect(rownames(WholeBlood$indata$C),rownames(WholeBlood$T));gene <- intersect(gene,WholeBlood$scDataSimulate$markerslist$gene);length(gene)
WholeBlood$indata$T <- WholeBlood$indata$T[gene,]
WholeBlood$indata$T <- WholeBlood$indata$T[which(rowSums(WholeBlood$indata$T) > 0),]
length(which(rowSums(WholeBlood$indata$T) > 0))
gene <- gene[which(rowSums(WholeBlood$indata$T) > 0)];length(gene)
WholeBlood$indata$C <- WholeBlood$indata$C[gene,]

WholeBlood$indata$T <- as.matrix(WholeBlood$indata$T)
WholeBlood$indata$C <- as.matrix(WholeBlood$indata$C)
WholeBlood$indata$P <- as.matrix(SumEqual_1(WholeBlood$indata$P))
WholeBlood$indata$C_ref <- WholeBlood$indata$C

WholeBlood$T ['BLK',]
WholeBlood$indata$T ['BLK',]

saveRDS(WholeBlood,"F:/wangchenqi/CDSC/3realBulk/bulkdata/WholeBlood_5pbmcs.rds")


#--------------------------------------
rm (list=ls ())
setwd("F:/wangchenqi/CDSC/3realBulk")
source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

#--------------------------------------
a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/3realBulk/bulkData"));a
a = c("WholeBlood_lm22","WholeBlood_3pbmcs","WholeBlood_5pbmcs")
STRING_name = a[1]; STRING_name
# STRING_name <- c("WholeBlood")
getwd()
bulkData <- readRDS(paste(getwd(),"/bulkData/",STRING_name,".rds",sep=""))
result <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))
result$lasso$p
# MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,     
#                      "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" , "EPIC",      
#                      "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
# 
# result <- CountAllResults(result,MyMethodName)
result$all
nrow(result$all)
# is.na(result$Linseed);is.na(result$celldistinguisher)
# result$Linseed <- result$Linseed
# result$Linseed <- NULL
# result$CellDistinguisher <- result$celldistinguisher
# result$celldistinguisher <- NULL
all(rownames(bulkData$indata$T) == rownames(bulkData$indata$C))
all(colnames(bulkData$indata$T) == colnames(bulkData$indata$P))
all(colnames(bulkData$indata$C) == rownames(bulkData$indata$P))

dim(bulkData$T)
dim(bulkData$P)

# result <- CountAllResults(result,MyMethodName)
# result$all
# nrow(result$all)
# STRING_name
# saveRDS(result,paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))
#关于C的ct的热图
STRING_name
# map = GetCorMatrix(result$CDSC3$dec$p,bulkData$indata$P,matrix = "p");map
map = GetCorMatrixLM22(result$CDSC3$dec$c,
                       bulkData$indata$C,
                       bulkData$indata$C_ref, matrix = "c");map
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
library(ggplot2)
# eoffice::topptx(plot_pheatmp,filename = paste(
# "F:/wangchenqi/CDSC/pictures/class4_samples_",STRING_name,".pptx",sep=""))
eoffice::topptx(plot_pheatmp,filename = paste(
  "F:/wangchenqi/CDSC/pictures/class4_NO_ct_",STRING_name,".pptx",sep=""))

# saveRDS(result,paste(getwd(),"/result/result_",STRING_name,".rds",sep=""))
# para_lambda_123 <- readRDS(paste(getwd(),"/paramater/para_lambda_",STRING_name,"_123.rds",sep=""))
para_lambda_12 <- para_lambda_44_8[which(para_lambda_44_8$lambdaC == 0),]
# result$all
result$CDSC3$result
result$CDSC2$result

# para_lambda_baron_Muraro_12 <- readRDS(
#   paste(getwd(),"/paramater/para_lambda_",STRING_name,"_12.rds",sep=""))
# para_lambda_baron_Muraro_123 <- readRDS(
#   paste(getwd(),"/paramater/para_lambda_",STRING_name,"_123.rds",sep=""))
# AllResult <- result$all;result$all

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")
# WholeBlood <- bulkData
bulkData <- WholeBlood
#-------CDSC3------
result <- list()
result$CDSC3$seedd = 44
result$CDSC3$TerCondition = 10^-8
result$CDSC3$Ss <- SM(t(bulkData$indata$T))
result$CDSC3$Sg <- SM(bulkData$indata$T)

result$CDSC3$lambda1 <- 1e-04
result$CDSC3$lambda2 <- 1e-02
result$CDSC3$lambdaC <- 1e+02

library(dplyr)

result$CDSC3$dec = CDSC_3(bulkData$indata$T, bulkData$indata$C_ref, dim(bulkData$indata$C_ref)[2], 
                          result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                          result$CDSC3$TerCondition,result$CDSC3$seedd,
                          result$CDSC3$Ss,result$CDSC3$Sg,all_number = 3000)

result$CDSC3$result <- calculate_pbmcs_result(result$CDSC3$dec$c,result$CDSC3$dec$p,
                                        bulkData$indata$T,bulkData$indata$C,bulkData$indata$C_ref,bulkData$indata$P,
                                        result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                                        number_iter = result$CDSC3$dec$jump,
                                        seedd = result$CDSC3$seedd,
                                        TerCondition = result$CDSC3$TerCondition)
result$CDSC3$result

#-------CDSC2------
result$CDSC2$lambda1 <- 1e-02
result$CDSC2$lambda2 <- 1e-03
result$CDSC2$lambdaC <- 0
result$CDSC2$dec = CDSC_2(bulkData$indata$T,  dim(bulkData$indata$C_ref)[2], 
                          result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                          result$CDSC3$TerCondition,result$CDSC3$seedd,
                          result$CDSC3$Ss,result$CDSC3$Sg,all_number = 3000)
colSums(result$CDSC2$dec$p)
result$CDSC2$result <- calculate_pbmcs_result(result$CDSC2$dec$c,result$CDSC2$dec$p,
                                        bulkData$indata$T, bulkData$indata$C, bulkData$indata$C_ref, bulkData$indata$P,
                                        result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                                        number_iter = result$CDSC2$dec$jump,
                                        seedd = result$CDSC3$seedd,
                                        TerCondition = result$CDSC3$TerCondition)
result$CDSC2$result
# bulkData <- WholeBlood
# result <- NULL
#----------NNLS----------
require(nnls)
result$NNLS$p = do.call(cbind.data.frame,lapply(apply(bulkData$indata$T,2,function(x) nnls::nnls(as.matrix(bulkData$indata$C_ref),x)),   function(y) y$x))
result$NNLS$p = apply(result$NNLS$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$NNLS$p) <- colnames(bulkData$indata$C_ref)
result$NNLS$p <- MergeCellType(result$NNLS$p)
result$NNLS$result = getPearsonRMSE(result$NNLS$p,bulkData$indata$P)
result$NNLS$result

#--------OLS------------
result$OLS$p = apply(bulkData$indata$T,2,function(x) lm(x ~ as.matrix(bulkData$indata$C_ref))$coefficients[-1])
result$OLS$p = apply(result$OLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$OLS$p = apply(result$OLS$p,2,function(x) x/sum(x)) #explicit STO constraint
diag(cor(result$OLS$p,bulkData$indata$P))
diag(cor(t(result$OLS$p),t(bulkData$indata$P)))
rownames(result$OLS$p) <- unlist(lapply(strsplit(rownames(result$OLS$p),")"),function(x) x[2]))
result$OLS$p <- MergeCellType(result$OLS$p)
colSums(result$OLS$p)
result$OLS$result = getPearsonRMSE(result$OLS$p,bulkData$indata$P)
result$OLS$result


#---------FARDEEP-----------
library(FARDEEP)
#result_FARDEEP = t(FARDEEP(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = t(FARDEEP::fardeep(bulkData$indata$C_ref, bulkData$indata$T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = apply(result$FARDEEP$p,2,function(x) x/sum(x)) #explicit STO constraint
diag(cor(result$FARDEEP$p,bulkData$indata$P))
diag(cor(t(result$FARDEEP$p),t(bulkData$indata$P)))
result$FARDEEP$p <- MergeCellType(result$FARDEEP$p)
colSums(result$OLS$p)
result$FARDEEP$result = getPearsonRMSE(result$FARDEEP$p,bulkData$indata$P)
result$FARDEEP$result


#----------CIBERSORT-----------------
source("F:/wangchenqi/CDSC/CIBERSORT.R")
result$CIBERSORT$p = CIBERSORT(sig_matrix = bulkData$indata$C_ref, mixture_file = bulkData$indata$T, QN = FALSE)
result$CIBERSORT$p = t(result$CIBERSORT$p[,1:(ncol(result$CIBERSORT$p)-3)])
diag(cor(result$CIBERSORT$p,bulkData$indata$P))
diag(cor(t(result$CIBERSORT$p),t(bulkData$indata$P)))
result$CIBERSORT$p <- MergeCellType(result$CIBERSORT$p)
colSums(result$CIBERSORT$p)
result$CIBERSORT$result = getPearsonRMSE(result$CIBERSORT$p,bulkData$indata$P)
result$CIBERSORT$result


#----------------"DeconRNASeq------
#nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
#datasets and reference matrix: signatures, need to be non-negative. 
#"use.scale": whether the data should be centered or scaled, default = TRUE
unloadNamespace("Seurat") #needed for PCA step
library(pcaMethods) #needed for DeconRNASeq to work
result$deconRNASeq$p = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(bulkData$indata$T), 
                                                  signatures = as.data.frame(bulkData$indata$C_ref), 
                                                  proportions = NULL, 
                                                  checksig = FALSE, 
                                                  known.prop = FALSE, 
                                                  use.scale = FALSE, 
                                                  fig = FALSE)$out.all)
colnames(result$deconRNASeq$p) = colnames(bulkData$indata$T)
# require(Seurat)
result$deconRNASeq$p <- MergeCellType(result$deconRNASeq$p)
colSums(result$deconRNASeq$p)
result$deconRNASeq$result = getPearsonRMSE(result$deconRNASeq$p,bulkData$indata$P)
result$deconRNASeq$result


#-----------RLR----------------
require(MASS)
result$RLR$p = do.call(cbind.data.frame,lapply(apply(bulkData$indata$T,2,function(x) 
  MASS::rlm(x ~ as.matrix(bulkData$indata$C_ref), maxit=100)), function(y) y$coefficients[-1]))
result$RLR$p = apply(result$RLR$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$RLR$p = apply(result$RLR$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$RLR$p) <- unlist(lapply(strsplit(rownames(result$RLR$p),")"),function(x) x[2]))
result$RLR$p <- MergeCellType(result$RLR$p)
colSums(result$RLR$p )
result$RLR$result = getPearsonRMSE(result$RLR$p,bulkData$indata$P)
result$RLR$result


#-----------DCQ------------------
#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default
require(ComICS)
result$DCQ$p = t(ComICS::dcq(reference_data = bulkData$indata$C_ref, 
                             mix_data = bulkData$indata$T, 
                             marker_set = as.data.frame(row.names(bulkData$indata$C_ref)) , 
                             alpha_used = 0.05, 
                             lambda_min = 0.2, 
                             number_of_repeats = 10)$average)
result$DCQ$p = apply(result$DCQ$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DCQ$p = apply(result$DCQ$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DCQ$p <- MergeCellType(result$DCQ$p)
colSums(result$DCQ$p)
result$DCQ$result = getPearsonRMSE(result$DCQ$p,bulkData$indata$P)
result$DCQ$result


#-----------elastic_net-----------
#standardize = TRUE by default. lambda=NULL by default 
require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
result$elastic_net$p = apply(bulkData$indata$T, 2, function(z) 
  coef(glmnet::glmnet(x = as.matrix(bulkData$indata$C_ref), 
                      y = z, 
                      alpha = 0.2, 
                      standardize = TRUE, 
                      lambda = glmnet::cv.glmnet(as.matrix(bulkData$indata$C_ref), z)$lambda.1se))[1:ncol(bulkData$indata$C_ref)+1,])
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) x/sum(x)) #explicit STO constraint
result$elastic_net$p <- MergeCellType(result$elastic_net$p)
colSums(result$elastic_net$p)
result$elastic_net$result = getPearsonRMSE(result$elastic_net$p,bulkData$indata$P)
result$elastic_net$result


#----------ridge----------------
# alpha=0
require(glmnet)
result$ridge$p = apply(bulkData$indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(bulkData$indata$C_ref), y = z, alpha = 0, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(bulkData$indata$C_ref), z)$lambda.1se))[1:ncol(bulkData$indata$C_ref)+1,])
result$ridge$p = apply(result$ridge$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$ridge$p = apply(result$ridge$p,2,function(x) x/sum(x)) #explicit STO constraint
result$ridge$p <- MergeCellType(result$ridge$p)
colSums(result$ridge$p)
result$ridge$result = getPearsonRMSE(result$ridge$p,bulkData$indata$P)
result$ridge$result


#----------lasso-----------------
#alpha=1; shrinking some coefficients to 0. 
require(glmnet)
result$lasso$p = apply(bulkData$indata$T, 2, function(z) 
  coef(glmnet::glmnet(x = as.matrix(bulkData$indata$C_ref), 
                      y = z, 
                      alpha = 0.5, 
                      standardize = TRUE, 
                      lambda = glmnet::cv.glmnet(as.matrix(bulkData$indata$C_ref), z)$lambda.1se))[1:ncol(bulkData$indata$C_ref)+1,])
result$lasso$p = apply(result$lasso$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$lasso$p = apply(result$lasso$p,2,function(x) x/sum(x)) #explicit STO constraint
result$lasso$p[which(is.na(result$lasso$p) == TRUE)] <- 0 #Needed for models where glmnet drops all terms of a model and fit an intercept-only model (very unlikely but possible).
length(colSums(result$lasso$p) == 1)
result$lasso$p <- MergeCellType(result$lasso$p)
colSums(result$lasso$p)
result$lasso$result = getPearsonRMSE(result$lasso$p,bulkData$indata$P)
result$lasso$result

#----------EPIC----------
require(EPIC)
#for PBMCs
bulkData$marker_distrib = bulkData$scDataSimulate$markerslist[bulkData$scDataSimulate$markerslist$gene %in% rownames(bulkData$indata$C_ref),]
markers = as.character(bulkData$marker_distrib$gene)

C_EPIC <- list()
common_CTs <- colnames(bulkData$indata$C_ref)
C_EPIC[["sigGenes"]] <- rownames(bulkData$indata$C_ref[markers,common_CTs])
C_EPIC[["refProfiles"]] <- as.matrix(bulkData$indata$C_ref[markers,common_CTs])
C_EPIC[["refProfiles.var"]] <- bulkData$scDataSimulate$refProfiles.var[markers,common_CTs]

#for LM22
# markers = as.character(rownames(bulkData$indata$T))
# 
# C_EPIC <- list()
# common_CTs <- intersect(colnames(bulkData$indata$C_ref), rownames(bulkData$indata$P))
# C_EPIC[["sigGenes"]] <- rownames(bulkData$indata$C_ref[markers,common_CTs])
# C_EPIC[["refProfiles"]] <- as.matrix(bulkData$sig$lm22_gene_CT[markers,common_CTs])
# C_EPIC[["refProfiles.var"]] <- bulkData$sig$lm22_gene_CT[markers,common_CTs]

result$EPIC$p <- t(EPIC::EPIC(bulk=as.matrix(bulkData$indata$T), reference=C_EPIC, withOtherCells=TRUE, scaleExprs=FALSE)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices
result$EPIC$p = result$EPIC$p[!rownames(result$EPIC$p) %in% "otherCells",]
result$EPIC$p <- MergeCellType(result$EPIC$p)
colSums(result$EPIC$p )
result$EPIC$result = getPearsonRMSE(result$EPIC$p,bulkData$indata$P)
result$EPIC$result


source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")



#-----complete deconvolution methods————----------
#--------DSA-------------------
# a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/3realBulk/bulkData"));a
# STRING_name = a[12]; STRING_name
# # STRING_name <- c("WholeBlood")
# getwd()
# bulkData <- readRDS(paste(getwd(),"/bulkData/",STRING_name,".rds",sep=""))
# result <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))
# result$all

require(CellMix)
md = bulkData$scDataSimulate$markerslist[rownames(bulkData$indata$T),]
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(rownames(md)),as.character(md$CT),list)

#FOR LM22
# md = bulkData$sig$lm22_gene_CT[rownames(bulkData$indata$T),]
# dim(md)
# length(which(rowSums(md) > 0));md = md[which(rowSums(md) > 0),]
# ctnames <- colnames(md)
# marker = NULL
# for (i in ctnames) {
#   mdmatrix <- matrix(data = i,nrow = length(which(md[[i]] == 1)))
#   rownames(mdmatrix) <- rownames(md)[which(md[[i]] == 1)]
#   marker = rbind(marker,mdmatrix)
#   md = md[-which(md[[i]] == 1),]
# }
# colnames(marker) = c("CT");marker = as.data.frame(marker)
# require(CellMix)
# md= marker
# ML = CellMix::MarkerList()
# ML@.Data <- tapply(as.character(rownames(md)),as.character(md$CT),list)

result$DSA = CellMix::ged(as.matrix(bulkData$indata$T), ML, method = "DSA", log = TRUE)
result$DSA$c = result$DSA@fit@W
result$DSA$p = result$DSA@fit@H
result$DSA$t = result$DSA$c%*%result$DSA$p

result$DSA$c <- MergeCellType(result$DSA$c,'c')
result$DSA$p <- MergeCellType(result$DSA$p,'p')

result$DSA$result$c = getPearsonRMSE(result$DSA$c,bulkData$indata$C)
result$DSA$result$p = getPearsonRMSE(result$DSA$p,bulkData$indata$P)
result$DSA$result$t = getPearsonRMSE(result$DSA$t,bulkData$indata$T)

result$DSA$result$all <- cbind(result$DSA$result$p,result$DSA$result$c)
result$DSA$result$all <-  cbind(result$DSA$result$all, result$DSA$result$t)
colnames(result$DSA$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$DSA$result <- result$DSA$result$all
result$DSA$result

#---------ssKL------------------

result$ssKL <- CellMix::ged(as.matrix(bulkData$indata$T), ML, method = "ssKL", sscale = FALSE, maxIter=500, log = FALSE)
result$ssKL$c = result$ssKL@fit@W
result$ssKL$p = result$ssKL@fit@H
result$ssKL$t = result$ssKL$c%*%result$ssKL$p

result$ssKL$c <- MergeCellType(result$ssKL$c,'c')
result$ssKL$p <- MergeCellType(result$ssKL$p,'p')

result$ssKL$result$c = getPearsonRMSE(result$ssKL$c,bulkData$indata$C)
result$ssKL$result$p = getPearsonRMSE(result$ssKL$p,bulkData$indata$P)
result$ssKL$result$t = getPearsonRMSE(result$ssKL$t,bulkData$indata$T)

result$ssKL$result$all <- cbind(result$ssKL$result$p,result$ssKL$result$c)
result$ssKL$result$all <-  cbind(result$ssKL$result$all,result$ssKL$result$t)
colnames(result$ssKL$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssKL$result <- result$ssKL$result$all
result$ssKL$result

#----------ssFrobenius-----------------

result$ssFrobenius <- CellMix::ged(as.matrix(bulkData$indata$T), ML, method = "ssFrobenius", sscale = TRUE, maxIter = 500, log = FALSE) #equivalent to coef(CellMix::ged(T,...)
result$ssFrobenius$c = result$ssFrobenius@fit@W
result$ssFrobenius$p = result$ssFrobenius@fit@H
result$ssFrobenius$t = result$ssFrobenius$c%*%result$ssFrobenius$p

result$ssFrobenius$c <- MergeCellType(result$ssFrobenius$c,'c')
result$ssFrobenius$p <- MergeCellType(result$ssFrobenius$p,'p')

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

result$deconf$c <- MergeCellType(result$deconf$c,'c')
result$deconf$p <- MergeCellType(result$deconf$p,'p')

result$deconf$result$c = getPearsonRMSE(result$deconf$c,bulkData$indata$C)
result$deconf$result$p = getPearsonRMSE(result$deconf$p,bulkData$indata$P)
result$deconf$result$t = getPearsonRMSE(result$deconf$t,bulkData$indata$T)

result$deconf$result$all <- cbind(result$deconf$result$p,result$deconf$result$c)
result$deconf$result$all <-  cbind(result$deconf$result$all,result$deconf$result$t)
colnames(result$deconf$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$deconf$result <- result$deconf$result$all
result$deconf$result


#------TOAST + NMF-------
# TOAST方法宝函数只能输入基因大于样本的数据。。
require(DeCompress)
source("F:/wangchenqi/CDSC/DeCompress.R")# csDeCompress_my

dim(bulkData$indata$T)
result$TOAST$toast.nmf <- csDeCompress_my(Y_raw = bulkData$indata$T,
                                                   K = dim(bulkData$indata$C_ref)[2],
                                                   nMarker = nrow(bulkData$indata$T),
                                                   FUN = nmfOut,
                                                   TotalIter = 30)
require(NMF)
result$TOAST$fin.nmf = nmf(x = bulkData$indata$T,
                           rank = dim(bulkData$indata$C_ref)[2])
result$TOAST$ppp = t(coef(result$TOAST$fin.nmf))
result$TOAST$ppp = t(apply(result$TOAST$ppp,1,function(c) c/sum(c)))

result$TOAST$p <- t(result$TOAST$ppp)
result$TOAST$c <- basis(result$TOAST$fin.nmf)
result$TOAST$t <- result$TOAST$c%*%result$TOAST$p
nmf.res = list(prop = result$TOAST$ppp,
               sig = basis(result$TOAST$fin.nmf))

rownames(bulkData$indata$P)
# labels <- Row_label(t(result$TOAST$p),t(bulkData$indata$P));labels
labels <- Row_label(result$TOAST$c,bulkData$indata$C_ref);labels
rownames(result$TOAST$p) <- labels
colnames(result$TOAST$c) <- labels

result$TOAST$c <- MergeCellType(result$TOAST$c,'c')
result$TOAST$p <- MergeCellType(result$TOAST$p,'p')

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
                                               n.types = dim(bulkData$indata$C_ref)[2],
                                               scree = 'drop',
                                               logTransform = F)

names(result$Linseed$Linseed.rs) = names(result$Linseed$nmf.res)
result$Linseed$c <- result$Linseed$Linseed.rs[[2]]
result$Linseed$p <- result$Linseed$Linseed.rs[[1]]
result$Linseed$t <- result$Linseed$c%*%result$Linseed$p
rownames(bulkData$indata$P)
labels <- Row_label(result$Linseed$c, bulkData$indata$C_ref);labels
rownames(result$Linseed$p) <- labels
colnames(result$Linseed$c) <- labels

result$Linseed$c <- MergeCellType(result$Linseed$c,'c')
result$Linseed$p <- MergeCellType(result$Linseed$p,'p')

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
  numCellClasses = dim(bulkData$indata$C_ref)[2],
  minDistinguisherAlternatives=1,
  maxDistinguisherAlternatives=100,
  minAlternativesLengthsNormalized=0.5,
  expressionQuantileForScale = 0.75,
  expressionQuantileForFilter=0.999,
  expressionConcentrationRatio=0.333,
  probesWithGenesOnly = FALSE,
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
result$CellDistinguisher$p <- result$CellDist.deconv$sampleComposition
result$CellDistinguisher$c <- result$CellDist.deconv$cellSubclassSignatures
result$CellDistinguisher$t <- result$CellDistinguisher$c%*%result$CellDistinguisher$p
rownames(bulkData$indata$P)
labels <- Row_label(result$CellDistinguisher$c, bulkData$indata$C_ref);labels
rownames(result$CellDistinguisher$p) <- labels
colnames(result$CellDistinguisher$c) <- labels

result$CellDistinguisher$c <- MergeCellType(result$CellDistinguisher$c,'c')
result$CellDistinguisher$p <- MergeCellType(result$CellDistinguisher$p,'p')

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
MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,
                  "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",
                  "MuSiC","Bisque","SCDC", "DWLS",
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
result <- CountAllResults(result,MyMethodName)
result$all;dim(result$all)
STRING_name
paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep="")
saveRDS(result,paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

#----------find Paramater ///----------
lambda1 <- c(0,10^-4,10^-3,0.005,10^-2,0.05,10^-1)
lambda2 <- c(0,10^-4,10^-3,0.005,10^-2,0.05,10^-1)
lambdaC <- c(0,10^-1,10^0,5,10^1,50,10^2,10^3)
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
for (dir_i in 4:length(lambda1)){
  for (dir_j in 1:length(lambda2)){
    for (dir_k in 1:length(lambdaC)){
      result_CDSC = CDSC_3(bulkData$indata$T, bulkData$indata$C_ref, dim(bulkData$indata$C_ref)[2], 
                           lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k],
                           bulkData$indata$TerCondition,bulkData$seedd,
                           bulkData$Ss,bulkData$Sg,all_number = 3000)
      result_para_c[[num]] <- result_CDSC[[1]]
      result_para_p[[num]] <- result_CDSC[[2]]
      number_iter[num] <- result_CDSC[[3]]
      
      pearson_para_c[[num]] <- cor(result_para_c[[num]],bulkData$indata$C)
      result1 = NULL
      if(!all(is.na(pearson_para_c[[num]]) == FALSE) ){
        break
      }
      result1 <- calculate_pbmcs_result(result_para_c[[num]],result_para_p[[num]],
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
saveRDS(para_lambda_44_8,
        paste(getwd(),"/paramater/para_lambda_",STRING_name,"_noCPM+.rds",sep=""))


# #-----热图----
# rm (list=ls ())
# setwd("F:/wangchenqi/CDSC/3realBulk")
# source("F:/wangchenqi/CDSC/CDSC.R")
# source("F:/wangchenqi/CDSC/CDSC_expand.R")
# 
# a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/3realBulk/bulkData"));a
# a = c("WholeBlood_lm22","WholeBlood_3pbmcs","WholeBlood_5pbmcs")
# 
# MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,     
#                      "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,'EPIC',       
#                      "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
# 
# bulkData = NULL
# result = NULL 
# for(i in 1:length(a)){
#   STRING_name = a[i]; STRING_name
#   # bulkData[[i]] <- readRDS(paste(getwd(),"/bulkData/",STRING_name,".rds",sep=""))
#   result[[i]] <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,".rds",sep=""))
# }
# map = NULL
# for(i in 1:length(a)){
#   # result[[i]] <- CountAllResults(result[[i]],MyMethodName)
#   map <- cbind(map,result[[i]]$all$Peason_to_P)
# }
# rownames(map) <- MyMethodName;colnames(map) <- a
# # map[which(map == 0)] <- NA;map
# rownames(map) <- MyMethodName;map
# library(pheatmap)#c("navy", "white", "firebrick3") color
# plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
#                          ,cellwidth = 42, cellheight = 12,gaps_row = c(11)
#                          ,breaks = c(seq(0.4,1,by = 0.02))
#                          ,display_numbers = TRUE, number_format = "%.3f"
#                          ,border_color = "black",color = colorRampPalette(c( "white", "firebrick3"))(30)
#                          ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
#                          ,main = "Pearson to P")
# 
# plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
#                          ,cellwidth = 42, cellheight = 12,gaps_row = c(11)
#                          ,display_numbers = TRUE, number_format = "%.3f"
#                          ,border_color = "black",color = colorRampPalette(c( "white", "navy"))(70)
#                          ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
#                          ,main = "RMSE to P")
# 
# plot_pheatmp <- pheatmap(map[c(1,12:19),],cluster_row = FALSE,cluster_cols = FALSE
#                          ,cellwidth = 30, cellheight = 12#,gaps_row = c(12, 17)
#                          ,breaks = c(seq(0.4,1,by = 0.02))
#                          ,display_numbers = TRUE, number_format = "%.3f"
#                          ,border_color = "black",color = colorRampPalette(c("white", "firebrick3"))(30)
#                          ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
#                          ,main = "Pearson to C")
# 
# plot_pheatmp <- pheatmap(map[c(1,12:19),],cluster_row = FALSE,cluster_cols = FALSE
#                          ,cellwidth = 42, cellheight = 12#,gaps_row = c(10)
#                          ,display_numbers = TRUE, number_format = "%.0f"
#                          ,border_color = "black",color = colorRampPalette(c( "white", "navy"))(70)
#                          ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
#                          ,main = "RMSE to C")
# apply(map, 1, mean)
# plot_pheatmp;
# library(ggplot2)
# eoffice::topptx(plot_pheatmp,filename = 
#                   "F:/wangchenqi/CDSC/pictures/class4_RMSE_P.pptx")
# 
# # 泛化能力的箱线图--p----
# mapBox <- NULL
# mapBox =  data.frame(pearson = map[1,], method = rownames(map)[1],row.names = NULL)
# nn1 <- nrow(mapBox)
# for (i in 1:length(rownames(map))) {
#   mapBox <- rbind(mapBox, data.frame(pearson = map[i,], method = rownames(map)[i],row.names = NULL))
# }
# 
# mapBox <- mapBox[-(1:nn1),]
# library(ggplot2)#coef=1e30,
# plot_class1 =  ggplot(mapBox, aes(factor(method,levels=MyMethodName),pearson,fill=method)) +  # background
#   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot
#   geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot
#   geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#   # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#   # scale_colour_manual(values = c("black"))
#   labs(x="Methods",y = "Pearson")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# # ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# plot_class1;
# eoffice::topptx(plot_class1,filename = 
#                   "F:/wangchenqi/CDSC/pictures/class4_fanhua_p.pptx")
# 
# # 泛化能力的箱线图----c-----
# mapBox <- NULL
# MyMethodName2 <-  c("CDSC3","CDSC2", "DSA","ssKL" ,"ssFrobenius",
#                     "TOAST","Linseed","CellDistinguisher")
# map2 <- map[MyMethodName2,]
# mapBox =  data.frame(pearson = map2[1,], method = rownames(map2)[1],row.names = NULL)
# nn1 <- nrow(mapBox)
# for (i in 1:length(rownames(map2))) {
#   mapBox <- rbind(mapBox, data.frame(pearson = map2[i,], method = rownames(map2)[i],row.names = NULL))
# }
# 
# mapBox <- mapBox[-(1:nn1),]
# library(ggplot2)
# plot_class1 =  ggplot(mapBox, aes(factor(method,levels=MyMethodName2),pearson,fill=method)) +  # background
#   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot
#   geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot
#   geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#   # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#   # scale_colour_manual(values = c("black"))
#   labs(x="Methods",y = "Pearson")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# # ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# plot_class1;
# eoffice::topptx(plot_class1,filename = "F:/wangchenqi/CDSC/pictures/class4_fanhua_c.pptx")
# 
# #------------ALL RESULT-----
# rm (list=ls ())
# setwd("F:/wangchenqi/CDSC/3realBulk")
# source("F:/wangchenqi/CDSC/CDSC.R")
# source("F:/wangchenqi/CDSC/CDSC_expand.R")
# library(ggplot2)
# STRING_name <- c("WholeBlood")
# WH <- list()
# # WH$result_nsclc <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_nsclc.rds",sep=""))
# # WH$result_nsclc_sc <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_nsclc_sc.rds",sep=""))
# WH$result_lm22 <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_lm22.rds",sep=""))
# WH$result_3pbmcs <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_3pbmcs.rds",sep=""))
# WH$result_5pbmcs <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_5pbmcs.rds",sep=""))
# # WH$result_FL <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_FL.rds",sep=""))
# RS <- list()
# # RS$result_nsclc <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_nsclc.rds",sep=""))
# # RS$result_nsclc_sc <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_nsclc_sc.rds",sep=""))
# RS$result_lm22 <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_lm22.rds",sep=""))
# RS$result_3pbmcs <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_3pbmcs.rds",sep=""))
# RS$result_5pbmcs <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_5pbmcs.rds",sep=""))
# # RS$result_FL <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_FL.rds",sep=""))
# #-----------PLOT-----
# Refnames <- names(WH);Refnames
# # [1] "result_nsclc_sc" "result_lm22"     "result_3pbmcs"   "result_5pbmcs"   "result_FL" 
# Refnames <- c("result_lm22","result_3pbmcs","result_5pbmcs")
# names(RS)
# # names(Aresult)
# #-----第一种-----每个参考的P整体指标（整个矩阵一起算相似性）共五个值-----
# OverAllResult <- NULL
# OverAllResult <- WH[[Refnames[1]]]$all$Peason_to_P
# for (i in 1:length(names(RS))) {
#   OverAllResult <- cbind(OverAllResult, WH[[Refnames[i]]]$all$Peason_to_P)
# }
# OverAllResult <- OverAllResult[,-1]
# rownames(OverAllResult) <- rownames(WH[[Refnames[1]]]$all); colnames(OverAllResult) <- Refnames
# OverAllResult
# 
# OverAllResultC <- NULL
# OverAllResultC <- WH[[Refnames[1]]]$all$Peason_to_C
# for (i in 1:length(names(RS))) {
#   OverAllResultC <- cbind(OverAllResultC, WH[[Refnames[i]]]$all$Peason_to_C)
# }
# OverAllResultC <- OverAllResultC[,-1]
# rownames(OverAllResultC) <- rownames(WH[[Refnames[1]]]$all);colnames(OverAllResultC) <- Refnames
# OverAllResultC
# ### 箱线图
# MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,     
#                   "deconRNASeq","RLR","DCQ","elastic_net","ridge", "lasso" ,       
#                   "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
# MethodName <- names(WH[[Refnames[1]]]);MethodName=setdiff(MethodName, c("all","CellDist.deconv"));MethodName
# MatrixCref <- intersect(c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT","deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso","EPIC"),MyMethodName)
# ScCref <- intersect(c("MuSiC","Bisque","deconvSeq","SCDC","DWLS"),MyMethodName)
# NoCref <- intersect(c("CDSC2","DSA","ssKL","ssFrobenius","deconf","CDSeq","TOAST","Linseed","CellDistinguisher"),MyMethodName)
# 
# Oneref <- union(MatrixCref,NoCref)
# 
# myOneref <- intersect(Oneref,MethodName);#myOneref=setdiff(myOneref,c("CDSC2"))
# myTworef <- intersect(NoCref,MethodName);myTworef=union(c("CDSC3"),myTworef);#myTworef=setdiff(myTworef,c("CDSC2"))
# myMatrixCref <- intersect(MatrixCref,MethodName)
# myScCref <- intersect(ScCref,MethodName)
# myNoCref <- intersect(NoCref,MethodName)
# # 先给CDSC3和CDSC2方法赋予正确的labels，但是不同数据的labels顺序时不一样的。
# Refnames
# for (i in 1:length(Refnames)) {
#   ctlabels <- Row_label(WH[[Refnames[i]]]$CDSC3$dec$c,RS[[Refnames[i]]]$indata$C,leastnum=3)
#   rownames(WH[[Refnames[i]]]$CDSC3$dec$p) <- ctlabels
#   colnames(WH[[Refnames[i]]]$CDSC3$dec$c) <- ctlabels
#   ctlabels_ <- rownames(RS[[Refnames[i]]]$indata$P)
#   WH[[Refnames[i]]]$CDSC3$p <- WH[[Refnames[i]]]$CDSC3$dec$p[ctlabels_, ]
#   WH[[Refnames[i]]]$CDSC3$c <- WH[[Refnames[i]]]$CDSC3$dec$c[ ,ctlabels_]
# }
# for (i in 1:length(Refnames)) {
#   ctlabels <- Row_label(WH[[Refnames[i]]]$CDSC2$dec$c,RS[[Refnames[i]]]$indata$C,leastnum=3)
#   rownames(WH[[Refnames[i]]]$CDSC2$dec$p) <- ctlabels
#   colnames(WH[[Refnames[i]]]$CDSC2$dec$c) <- ctlabels
#   ctlabels_ <- rownames(RS[[Refnames[i]]]$indata$P)
#   WH[[Refnames[i]]]$CDSC2$p <- WH[[Refnames[i]]]$CDSC2$dec$p[ctlabels_, ]
#   WH[[Refnames[i]]]$CDSC2$c <- WH[[Refnames[i]]]$CDSC2$dec$c[ ,ctlabels_]
# }
# # 赋予一致的细胞类型顺序
# Refnames 
# ctlabels_ <- c("NK.cells" ,   "Monocytes"  , "B.cells" ,    "T.cells.CD4" ,"T.cells.CD8")
# RS_ok<-list()
# for(i_method in MethodName){
#   for (i_data in Refnames){
#     WH[[i_data]][[i_method]]$p <- WH[[i_data]][[i_method]]$p[ctlabels_, ]
#     WH[[i_data]][[i_method]]$c <- WH[[i_data]][[i_method]]$c[ ,ctlabels_]      
#     RS[[i_data]]$indata$P <- RS[[i_data]]$indata$P[ctlabels_, ]
#     RS[[i_data]]$indata$C <- RS[[i_data]]$indata$C[,ctlabels_]
#   }
# }
# # RS_ok$indata$P <- RS[[Refnames[1]]]$indata$P[ctlabels_, ]
# 
# Refnames 
# MethodName
# library(RColorBrewer)
# ## 部分反卷积-P-sample
# 
# PartSampleResult <- NULL
# plot_1 <- NULL
# for (i in Refnames) {
#   PartSampleResult[[i]] <- NULL
#   PartSampleResult[[i]] <- data.frame(pearson = diag(cor(WH[[i]]$CDSC3$p,RS[[i]]$indata$P)),
#                            method = c("CDSC3"),row.names = NULL)
#   nn1 <- nrow(PartSampleResult[[i]])
#   for (i_data in i) {
#     for (i_method in myOneref) {
#         PartSampleResult[[i]] <- rbind(PartSampleResult[[i]], 
#                                 data.frame(pearson = diag(cor(WH[[i_data]][[i_method]]$p,RS[[i_data]]$indata$P)),
#                                                      method = i_method,row.names = NULL))
#     }
#   }
#   PartSampleResult[[i]] <- PartSampleResult[[i]][-(1:nn1),]
#   names_1 <-myOneref 
#   names_1 [which(names_1 == c("CDSC3"))] <- c("CDSC3")
#   names_1 [which(names_1 == c("CDSC2"))] <- c("CDSC2")
#   PartSampleResult[[i]][which(PartSampleResult[[i]]$method == c("CDSC3")),2] <- c("CDSC3")
#   PartSampleResult[[i]][which(PartSampleResult[[i]]$method == c("CDSC2")),2] <- c("CDSC2")
#   plot_1[[i]] = ggplot(PartSampleResult[[i]], aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
#     # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
#     geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#BEBADA") + # Boxplot 
#     geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#     # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#     # sacle_color_brewer(palette='set1')+
#     labs(x="Methods",y = "Pearson")+
#     theme_classic()+
#     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white")) +   #panel.background = element_rect(fill = '#d8dfea')) # 底色
#     ggtitle(strsplit(i,"_")[[1]][2]) +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# }
# plot_1[["result_lm22"]]
# eoffice::topptx(plot_1[["result_lm22"]],filename =
#                   "F:/wangchenqi/CDSC/pictures/class4_no1-1.pptx")
# plot_1[["result_3pbmcs"]]
# eoffice::topptx(plot_1[["result_3pbmcs"]],filename =
#                   "F:/wangchenqi/CDSC/pictures/class4_no1-2.pptx")
# plot_1[["result_5pbmcs"]]
# eoffice::topptx(plot_1[["result_5pbmcs"]],filename =
#                   "F:/wangchenqi/CDSC/pictures/class4_no1-3.pptx")
# 
# #------画在一起----------
# PartSampleResultALL <- PartSampleResult[["result_lm22"]]
# PartSampleResultALL <- rbind(PartSampleResultALL,PartSampleResult[["result_3pbmcs"]])
# PartSampleResultALL <- rbind(PartSampleResultALL,PartSampleResult[["result_5pbmcs"]])
# 
# plot_1_all = ggplot(PartSampleResultALL, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
#   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
#   geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
#   geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#   # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#   # scale_colour_manual(values = c("black"))
#   labs(x="Methods",y = "Pearson")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white")) +   #panel.background = element_rect(fill = '#d8dfea')) # 底色
#   ggtitle("All") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# plot_1_all
# eoffice::topptx(plot_1_all,filename = 
#                   "F:/wangchenqi/CDSC/pictures/class4_no1-all_samples.pptx")
# 
# 
# colorRampPalette(brewer.pal(9,'Reds')[c(1,2,4,6,9)])(5)
# facebook = c("#3b5998","#6d84b4", "#afbdd4", "#d8dfea")
# google = c("#5380E4", "#E12A3C", "#FFBF03","#00B723")
# etsy = c("#F14000", "#67B6C3", "#F0DA47", "#EBEBE6", "#D0D0CB")
# twitter = c("#55ACEE", "#292f33", "#8899a6", "#e1e8ed")
# # 部分反卷积-P-celltype
# # PartCTResult <- NULL
# # PartCTResult <- data.frame(pearson = diag(cor(t(WH[[Refnames[1]]]$CDSC3$p),t(RS[[Refnames[1]]]$indata$P))),
# #                        method = c("CDSC3"),row.names = NULL)
# # nn2 <- nrow(PartCTResult)
# # for (i_data in Refnames) {
# #   for (i_method in myOneref) {
# #       PartCTResult <- rbind(PartCTResult,
# #                             data.frame(pearson = diag(cor(t( WH[[i_data]][[i_method]]$p),t(RS[[i_data]]$indata$P))),
# #                                              method = i_method,row.names = NULL))
# #   }
# # }
# # PartCTResult <- PartCTResult[-(1:nn2),]
# # PartCTResult[which(PartCTResult$method == c("CDSC3")),2] <- c("CDSC")
# #   ggplot(PartCTResult, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# 
# # ## 完全反卷积-P-sample
# # CompSampleResult <- NULL
# # CompSampleResult <- data.frame(pearson = diag(cor(WH[[Refnames[1]]]$CDSC3$p,RS[[Refnames[1]]]$indata$P)),
# #                            method = c("CDSC3"),row.names = NULL)
# # nn3 <- nrow(CompSampleResult)
# # for (i_data in Refnames) {
# #   for (i_method in myNoCref) {
# #     CompSampleResult <- rbind(CompSampleResult, 
# #                               data.frame(pearson = diag(cor(WH[[i_data]][[i_method]]$p,RS[[i_data]]$indata$P)),
# #                                                      method = i_method,row.names = NULL))
# #   }
# # }
# # CompSampleResult <- CompSampleResult[-(1:nn3),]
# # ggplot(CompSampleResult, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot 
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# # 
# # ## 完全反卷积-P-celltype
# # CompCTResult <- NULL
# # CompCTResult <- data.frame(pearson = diag(cor(t(WH[[Refnames[1]]]$CDSC3$p),t(RS[[Refnames[1]]]$indata$P))),
# #                            method = c("CDSC3"),row.names = NULL)
# # nn4 <- nrow(CompCTResult)
# # for (i_data in Refnames) {
# #   for (i_method in myTworef) {
# #     CompCTResult <- rbind(CompCTResult, 
# #                           data.frame(pearson = diag(cor(t(WH[[i_data]][[i_method]]$p),t(RS[[i_data]]$indata$P))),
# #                                                      method = i_method,row.names = NULL))    
# #   }
# # }
# # CompCTResult <- CompCTResult[-(1:nn4),]
# # ggplot(CompCTResult, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot 
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# 
# ## 完全反卷积-C-celltype
# CompCTResultC <- NULL
# plot_2 <- NULL
# for (i in Refnames) {
#   CompCTResultC[[i]] <- NULL
#     CompCTResultC[[i]] <- data.frame(pearson = diag(cor(WH[[i]]$CDSC3$c,RS[[i]]$indata$C)),
#                                       method = c("CDSC3"),row.names = NULL)
#   nn1 <- nrow(CompCTResultC[[i]])
#   for (i_data in i) {
#     for (i_method in myTworef) {
#       CompCTResultC[[i]] <- rbind(CompCTResultC[[i]], 
#                                      data.frame(pearson = diag(cor(WH[[i_data]][[i_method]]$c,RS[[i_data]]$indata$C)),
#                                                 method = i_method,row.names = NULL))
#     }
#   }
#   CompCTResultC[[i]] <- CompCTResultC[[i]][-(1:nn1),]
#   names_2 <- myTworef 
#   names_2 [which(names_2 == c("CDSC3"))] <- c("CDSC3")
#   names_2 [which(names_2 == c("CDSC2"))] <- c("CDSC2")
#   CompCTResultC[[i]][which(CompCTResultC[[i]]$method == c("CDSC3")),2] <- c("CDSC3")
#   CompCTResultC[[i]][which(CompCTResultC[[i]]$method == c("CDSC2")),2] <- c("CDSC2")
#   
#   plot_2[[i]] = ggplot(CompCTResultC[[i]], aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
#     # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
#     geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#91D1C27F") + # Boxplot 
#     geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#     # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#     # scale_colour_manual(values = c("black"))
#     labs(x="Methods",y = "Pearson")+
#     theme_classic()+
#     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white")) +   #panel.background = element_rect(fill = '#d8dfea')) # 底色
#     ggtitle(strsplit(i,"_")[[1]][2]) +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# }
# plot_2[["result_lm22"]]
# eoffice::topptx(plot_2[["result_lm22"]],filename = "F:/wangchenqi/CDSC/pictures/class4_no1-4.pptx")
# plot_2[["result_3pbmcs"]]
# eoffice::topptx(plot_2[["result_3pbmcs"]],filename = "F:/wangchenqi/CDSC/pictures/class4_no1-5.pptx")
# plot_2[["result_5pbmcs"]]
# eoffice::topptx(plot_2[["result_5pbmcs"]],filename = "F:/wangchenqi/CDSC/pictures/class4_no1-6.pptx")
# 
# CompCTResultCALL <- CompCTResultC[["result_lm22"]]
# CompCTResultCALL <- rbind(CompCTResultCALL,CompCTResultC[["result_3pbmcs"]])
# CompCTResultCALL <- rbind(CompCTResultCALL,CompCTResultC[["result_5pbmcs"]])
# 
# plot_2_all = ggplot(CompCTResultCALL, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
#   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
#   geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
#   geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#   # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#   # scale_colour_manual(values = c("black"))
#   labs(x="Methods",y = "Pearson")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white")) +   #panel.background = element_rect(fill = '#d8dfea')) # 底色
#   ggtitle("All") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# plot_2_all
# eoffice::topptx(plot_2_all,filename = "F:/wangchenqi/CDSC/pictures/class4_no1-all_ct.pptx")
# 
# #-----第二种-----先求平均，再求评价指标-----
# Refnames 
# MethodName 
# # FL参考数据只有3种细胞类型，少于真实的细胞类型数量，所以不予计算
# 
# WH_OK <- list()
# RS_OK <- list()
# # MethodName <- names(WH[[Refnames[1]]]);MethodName=setdiff(MethodName,c("all","lasso","CellDist.deconv"));MethodName
# MethodName <- names(WH[[Refnames[1]]]);MethodName=setdiff(MethodName, c("all","CellDist.deconv"));MethodName
# 
# for (i_method in MethodName) {
#   WH_OK[[i_method]]$p <- SumEqual_1(WH[[Refnames[1]]][[i_method]]$p)
#   for (i in 2:length(Refnames)){
#   WH_OK[[i_method]]$p <- WH_OK[[i_method]]$p + SumEqual_1(WH[[Refnames[i]]][[i_method]]$p)
#   }
#   WH_OK[[i_method]]$p = WH_OK[[i_method]]$p/5
#   WH_OK[[i_method]]$p <- WH_OK[[i_method]]$p[ctlabels_, ]
# }
# ## 由于不同参考的marker基因数，不一样，所以，无法进行此项实验
# # for (i_method in myNoCref) {
# #   WH_OK[[i_method]]$c <- SumEqual_1(WH[[Refnames[1]]][[i_method]]$c)
# #   for (i in 2:length(Refnames)){
# #     WH_OK[[i_method]]$c <- WH_OK[[i_method]]$c + SumEqual_1(WH[[Refnames[i]]][[i_method]]$c)
# #   }
# #   WH_OK[[i_method]]$c = WH_OK[[i_method]]$c/5
# #   WH_OK[[i_method]]$c <- WH_OK[[i_method]]$c[,ctlabels_ ]
# # 
# # }
# RS_OK$indata$P <- RS[[Refnames[1]]]$indata$P[ctlabels_, ]
# ##  算平均后的总体指标
# for(i_method in MethodName){
#   WH_OK[[i_method]]$result_p = getPearsonRMSE(WH_OK[[i_method]]$p,RS_OK$indata$P)
# }
# WH_OK$all_result <- NULL
# WH_OK$all_result <- WH_OK$CDSC3$result_p
# for (i_method in MethodName) {
#   WH_OK$all_result <- rbind(WH_OK$all_result,WH_OK[[i_method]]$result_p)
# };WH_OK$all_result <- WH_OK$all_result[-1,]
# rownames(WH_OK$all_result) <- MethodName
# WH_OK$all_result
# 
# ##  部分反卷积-P-sample
# EqualPartSampleResult <- NULL
# EqualPartSampleResult <- data.frame(pearson = diag(cor(WH_OK[[MethodName[1]]]$p,RS_OK$indata$P)),
#                                 method = MethodName[1],row.names = NULL)
# nn1 <- nrow(EqualPartSampleResult)
# for (i_method in myOneref) {
#   EqualPartSampleResult <- rbind(EqualPartSampleResult, data.frame(pearson = diag(cor(WH_OK[[i_method]]$p,RS_OK$indata$P)),
#                                                  method = i_method,row.names = NULL))
# }
# EqualPartSampleResult <-EqualPartSampleResult[-(1:nn1),]
# names_1 <-myOneref 
# names_1 [which(names_1 == c("CDSC3"))] <- c("CDSC3")
# names_1 [which(names_1 == c("CDSC2"))] <- c("CDSC2")
# EqualPartSampleResult[which(EqualPartSampleResult$method == c("CDSC3")),2] <- c("CDSC3")
# EqualPartSampleResult[which(EqualPartSampleResult$method == c("CDSC2")),2] <- c("CDSC2")
# plot_3 =  ggplot(EqualPartSampleResult, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
#   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
#   geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#4DBBD57F") + # Boxplot 
#   geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#   # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#   # scale_colour_manual(values = c("black"))
#   labs(x="Methods",y = "Pearson")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# # ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# plot_3
# eoffice::topptx(plot_3,filename = "F:/wangchenqi/CDSC/pictures/class4_no2-1.pptx")
# 
# 
# ##  部分反卷积-P-celltype
# # EqualPartCTResult <- NULL
# # EqualPartCTResult <- data.frame(pearson = diag(cor(t(WH_OK[[MethodName[1]]]$p),t(RS_OK$indata$P))),
# #                                 method = MethodName[1],row.names = NULL)
# # nn2 <- nrow(EqualPartCTResult)
# # for (i_method in myMatrixCref) {
# #   EqualPartCTResult <- rbind(EqualPartCTResult, data.frame(pearson = diag(cor(t(WH_OK[[i_method]]$p),t(RS_OK$indata$P))),
# #                                                            method = i_method,row.names = NULL))
# # }
# # EqualPartCTResult <-EqualPartCTResult[-(1:nn2),]
# # ggplot(EqualPartCTResult, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot 
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# 
# ##  完全反卷积-P-sample
# # EqualCompSampleResult <- NULL
# # EqualCompSampleResult <- data.frame(pearson = diag(cor(WH_OK[[MethodName[1]]]$p,RS_OK$indata$P)),
# #                                     method = MethodName[1],row.names = NULL)
# # nn3 <- nrow(EqualCompSampleResult)
# # for (i_method in myNoCref) {
# #   EqualCompSampleResult <- rbind(EqualCompSampleResult, data.frame(pearson = diag(cor(WH_OK[[i_method]]$p,RS_OK$indata$P)),
# #                                                                    method = i_method,row.names = NULL))
# # }
# # EqualCompSampleResult <-EqualCompSampleResult[-(1:nn3),]
# # ggplot(EqualCompSampleResult, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot 
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# # 
# # ##  完全反卷积-P-celltype
# # EqualCompCTResult <- NULL
# # EqualCompCTResult <- data.frame(pearson = diag(cor(t(WH_OK[[MethodName[1]]]$p),t(RS_OK$indata$P))),
# #                                 method = MethodName[1],row.names = NULL)
# # nn4 <- nrow(EqualCompCTResult)
# # for (i_method in myNoCref) {
# #   EqualCompCTResult <- rbind(EqualCompCTResult, data.frame(pearson = diag(cor(t(WH_OK[[i_method]]$p),t(RS_OK$indata$P))),
# #                                                            method = i_method,row.names = NULL))
# # }
# # EqualCompCTResult <-EqualCompCTResult[-(1:nn4),]
# # ggplot(EqualCompCTResult, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot 
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# 
# 
# #-----第三种-----先算指标；再平均-----
# # 各个参考计算出P矩阵的整体相似性指标之后，再对指标平均
# ThreeOverAllResult <- as.data.frame(round(as.numeric(rowSums(OverAllResult)/ncol(OverAllResult)),8))
# rownames(ThreeOverAllResult) <- rownames(OverAllResult)
# ThreeOverAllResult
# ThreeOverAllResultC <- as.data.frame(round(as.numeric(rowSums(OverAllResultC)/ncol(OverAllResult)),8))
# rownames(ThreeOverAllResultC) <- rownames(OverAllResultC)
# ThreeOverAllResultC
# 
# ## 部分反卷积-P-samples
# TwoPartSampleResult <- NULL
# TwoPartSampleResult$all <- NULL
# for (i_method in myOneref) {
#   TwoPartSampleResult[[i_method]] <- data.frame(pearson = diag(cor(WH[[Refnames[1]]]$CDSC3$p,RS[[Refnames[1]]]$indata$P)),
#                                                 row.names = NULL)
#   for (i_data in Refnames) {
#     TwoPartSampleResult[[i_method]] <- cbind(TwoPartSampleResult[[i_method]],
#                               data.frame(Pearson = diag(cor(WH[[i_data]][[i_method]]$p,RS[[i_data]]$indata$P)),
#                                          row.names = NULL))
#   }
#   TwoPartSampleResult[[i_method]] <- TwoPartSampleResult[[i_method]][,-1]
#   colnames(TwoPartSampleResult[[i_method]]) <- Refnames
#   TwoPartSampleResult[[i_method]] <- as.data.frame(round(as.numeric(rowSums(TwoPartSampleResult[[i_method]])/length(Refnames)), 8))
#   TwoPartSampleResult$all <- rbind(TwoPartSampleResult$all,
#                                     data.frame(pearson = TwoPartSampleResult[[i_method]],
#                                                method = i_method,row.names = NULL))
# }
# colnames(TwoPartSampleResult$all) <- c("pearson","method")
# TwoPartSampleResult$all[which(TwoPartSampleResult$all$method == c("CDSC3")),2] <- c("CDSC3")
# TwoPartSampleResult$all[which(TwoPartSampleResult$all$method == c("CDSC2")),2] <- c("CDSC2")
# 
# plot_4 = ggplot(TwoPartSampleResult$all, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
#   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
#   geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#00A0877F") + # Boxplot 
#   geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#   # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#   # scale_colour_manual(values = c("black"))
#   labs(x="Methods",y = "Pearson")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# # ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# plot_4;
# eoffice::topptx(plot_4,filename = "F:/wangchenqi/CDSC/pictures/class4_no3-1.pptx")
# 
# ## 部分反卷积-P-celltyps
# # TwoPartCTResult <- NULL
# # TwoPartCTResult$all <- NULL
# # for (i_method in myMatrixCref) {
# #   TwoPartCTResult[[i_method]] <- data.frame(pearson = diag(cor(t(WH[[Refnames[1]]]$CDSC3$p),t(RS[[Refnames[1]]]$indata$P))),
# #                                                 row.names = NULL)
# #   for (i_data in Refnames) {
# #     TwoPartCTResult[[i_method]] <- cbind(TwoPartCTResult[[i_method]],
# #                                              data.frame(Pearson = diag(cor(t(WH[[i_data]][[i_method]]$p),t(RS[[i_data]]$indata$P))),
# #                                                         row.names = NULL))
# #   }
# #   TwoPartCTResult[[i_method]] <- TwoPartCTResult[[i_method]][,-1]
# #   colnames(TwoPartCTResult[[i_method]]) <- Refnames
# #   TwoPartCTResult[[i_method]] <- as.data.frame(round(as.numeric(rowSums(TwoPartCTResult[[i_method]])/length(Refnames)), 8))
# #   TwoPartCTResult$all <- rbind(TwoPartCTResult$all,
# #                                    data.frame(pearson = TwoPartCTResult[[i_method]],
# #                                               method = i_method,row.names = NULL))
# # }
# # colnames(TwoPartCTResult$all) <- c("pearson","method")
# # ggplot(TwoPartCTResult$all, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot 
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# # 
# # ## 完全反卷积-P-samples
# # TwoCompSampleResult <- NULL
# # TwoCompSampleResult$all <- NULL
# # for (i_method in myNoCref) {
# #   TwoCompSampleResult[[i_method]] <- data.frame(pearson = diag(cor(WH[[Refnames[1]]]$CDSC2$p,RS[[Refnames[1]]]$indata$P)),
# #                                                 row.names = NULL)
# #   for (i_data in Refnames) {
# #     TwoCompSampleResult[[i_method]] <- cbind(TwoCompSampleResult[[i_method]],
# #                                              data.frame(Pearson = diag(cor(WH[[i_data]][[i_method]]$p,RS[[i_data]]$indata$P)),
# #                                                         row.names = NULL))
# #   }
# #   TwoCompSampleResult[[i_method]] <- TwoCompSampleResult[[i_method]][,-1]
# #   colnames(TwoCompSampleResult[[i_method]]) <- Refnames
# #   TwoCompSampleResult[[i_method]] <- as.data.frame(round(as.numeric(rowSums(TwoCompSampleResult[[i_method]])/length(Refnames)), 8))
# #   TwoCompSampleResult$all <- rbind(TwoCompSampleResult$all,
# #                                    data.frame(pearson = TwoCompSampleResult[[i_method]],
# #                                               method = i_method,row.names = NULL))
# # }
# # colnames(TwoCompSampleResult$all) <- c("pearson","method")
# # ggplot(TwoCompSampleResult$all, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot 
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# # 
# # ## 完全反卷积-P-celltyps
# # TwoCompCTResult <- NULL
# # TwoCompCTResult$all <- NULL
# # for (i_method in myNoCref) {
# #   TwoCompCTResult[[i_method]] <- data.frame(pearson = diag(cor(t(WH[[Refnames[1]]]$CDSC2$p),t(RS[[Refnames[1]]]$indata$P))),
# #                                             row.names = NULL)
# #   for (i_data in Refnames) {
# #     TwoCompCTResult[[i_method]] <- cbind(TwoCompCTResult[[i_method]],
# #                                          data.frame(Pearson = diag(cor(t(WH[[i_data]][[i_method]]$p),t(RS[[i_data]]$indata$P))),
# #                                                     row.names = NULL))
# #   }
# #   TwoCompCTResult[[i_method]] <- TwoCompCTResult[[i_method]][,-1]
# #   colnames(TwoCompCTResult[[i_method]]) <- Refnames
# #   TwoCompCTResult[[i_method]] <- as.data.frame(round(as.numeric(rowSums(TwoCompCTResult[[i_method]])/length(Refnames)), 8))
# #   TwoCompCTResult$all <- rbind(TwoCompCTResult$all,
# #                                data.frame(pearson = TwoCompCTResult[[i_method]],
# #                                           method = i_method,row.names = NULL))
# # }
# # colnames(TwoCompCTResult$all) <- c("pearson","method")
# # ggplot(TwoCompCTResult$all, aes(method,pearson,fill=method)) +  # background
# #   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
# #   geom_boxplot(notch = F,width=0.5,outlier.shape = NA) + # Boxplot 
# #   geom_jitter(shape=16, position=position_jitter(0.1)) + #plot Scatter diagram
# #   stat_summary(fun = mean, geom = "point", shape = 23, size=4) # add Mean
# 
# ## 完全反卷积-C-celltype
# TwoCompCTResultC <- NULL
# TwoCompCTResultC$all <- NULL
# for (i_method in myTworef) {
#   TwoCompCTResultC[[i_method]] <- data.frame(pearson = diag(cor(WH[[Refnames[1]]]$CDSC2$c,RS[[Refnames[1]]]$indata$C)),
#                                                 row.names = NULL)
#   for (i_data in Refnames) {
#     TwoCompCTResultC[[i_method]] <- cbind(TwoCompCTResultC[[i_method]],
#                                              data.frame(Pearson = diag(cor(WH[[i_data]][[i_method]]$c,RS[[i_data]]$indata$C)),
#                                                         row.names = NULL))
#   }
#   TwoCompCTResultC[[i_method]] <- TwoCompCTResultC[[i_method]][,-1]
#   colnames(TwoCompCTResultC[[i_method]]) <- Refnames
#   TwoCompCTResultC[[i_method]] <- as.data.frame(round(as.numeric(rowSums(TwoCompCTResultC[[i_method]])/length(Refnames)), 8))
#   TwoCompCTResultC$all <- rbind(TwoCompCTResultC$all,
#                                    data.frame(pearson = TwoCompCTResultC[[i_method]],
#                                               method = i_method,row.names = NULL))
# }
# colnames(TwoCompCTResultC$all) <- c("pearson","method")
# TwoCompCTResultC$all[which(TwoCompCTResultC$all$method == c("CDSC3")),2] <- c("CDSC3")
# TwoCompCTResultC$all[which(TwoCompCTResultC$all$method == c("CDSC2")),2] <- c("CDSC2")
# 
# plot_5 = ggplot(TwoCompCTResultC$all, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
#   # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
#   geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#00A0877F") + # Boxplot 
#   geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
#   # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
#   # scale_colour_manual(values = c("black"))
#   labs(x="Methods",y = "Pearson")+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
#   # ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
# plot_5
# eoffice::topptx(plot_5,filename = "F:/wangchenqi/CDSC/pictures/class4_no3-2+.pptx")
# ### 热图
# library(pheatmap)
# pheatmap(WholeBlood$sig$lm22,cluster_row = FALSE,cluster_cols = FALSE
#          ,cellwidth = 10, cellheight = 3,gaps_row = c(12, 17)
#          ,border_color = "black",color = colorRampPalette(c("navy", "white", "firebrick3"))(70)
#          ,main = "RMSE of P")
# 
