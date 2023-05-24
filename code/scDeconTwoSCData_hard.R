# 2 scSimulated data sets, find the deconvolution results

# ---
# title: "2scDeconvolution"
# author: "Chenqi Wang"
# date: "2022/3/12"
# output: html_document
# ---
# 

#-----hunman pancreas ///-------------

##
#------read raw data information--------
# rm (list=ls ())
# setwd("F:/wangchenqi/CDSC/2scDeconvolution")


# 人类肾脏数据有：Baron, Muraro, Segerstolpe, xin, Enge, wang

#-----human pancreas--------------------------
Baron <- list(
  data = readRDS("F:/wangchenqi/CDSC/dataset/Baron/process/arrange/sc_Baron.rds"),
  full_phenoData = readRDS("F:/wangchenqi/CDSC/dataset/Baron/process/arrange/Baron_phenoData.rds")
)
Muraro <- list(
  data = readRDS("F:/wangchenqi/CDSC/dataset/Muraro/process/arrange/sc_Muraro.rds"),
  full_phenoData = readRDS("F:/wangchenqi/CDSC/dataset/Muraro/process/arrange/Muraro_phenoData.rds")
)
Segerstolpe <- list(
  data = readRDS("F:/wangchenqi/CDSC/dataset/Segerstolpe/process/arrange/sc_Segerstolpe.rds"),
  full_phenoData = readRDS("F:/wangchenqi/CDSC/dataset/Segerstolpe/process/arrange/Segerstolpe_phenoData.rds")
)
# xin <- list(
#   data = readRDS("F:/wangchenqi/CDSC/dataset/xin/process/arrange/sc_xin.rds"),
#   full_phenoData = readRDS("F:/wangchenqi/CDSC/dataset/xin/process/arrange/xin_phenoData.rds")
# )
# Enge <- list(
#   data = readRDS("F:/wangchenqi/CDSC/dataset/Enge/process/arrange/sc_Enge.rds"),
#   full_phenoData = readRDS("F:/wangchenqi/CDSC/dataset/Enge/process/arrange/Enge_phenoData.rds")
# )

#---Mouse Retina----------------------------
# Macosko <- list(
#   data = readRDS("F:/wangchenqi/CDSC/dataset/Macosko/process/arrange/sc_Macosko.rds"),
#   full_phenoData = readRDS("F:/wangchenqi/CDSC/dataset/Macosko/process/arrange/Macosko_phenoData.rds")
# )
# Shekhar <- list(
#   data = readRDS("F:/wangchenqi/CDSC/dataset/shekhar/process/arrange/sc_shekhar.rds"),
#   full_phenoData = readRDS("F:/wangchenqi/CDSC/dataset/shekhar/process/arrange/shekhar_phenoData.rds")
# )
# #------simulate--------
# source("F:/wangchenqi/CDSC/CDSC.R")
# source("F:/wangchenqi/CDSC/CDSC_expand.R")

# #-----
# Baron
# Segerstolpe
# Muraro
# #----
# shekhar
# Macosko
# 
# a = c("Baron_Muraro","Muraro_Baron","Baron_Segerstolpe","Segerstolpe_Baron","Segerstolpe_Muraro",
#       "Muraro_Segerstolpe","Macosko_Shekhar","Shekhar_Macosko" );a

ComDate <- Combine_2scdata_hard( Shekhar,Macosko )
# saveRDS(ComDate,"F:/wangchenqi/CDSC/2scDeconvolution/data2sc_hard/Shekhar_Macosko.rds")
# saveRDS(ComDate,paste(getwd(),"/data2sc/",STRING_name,".rds",sep=""))

# #----START----
rm (list=ls ())
setwd("F:/wangchenqi/CDSC/2scDeconvolution")

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/2scDeconvolution/data2sc"));a
a = c("Baron_Muraro","Muraro_Baron","Baron_Segerstolpe","Segerstolpe_Baron","Segerstolpe_Muraro","Muraro_Segerstolpe","Macosko_Shekhar","Shekhar_Macosko" )
# a = c("Baron_Muraro","Muraro_Baron","Segerstolpe_Muraro","Muraro_Segerstolpe")
STRING_name = a[1]; STRING_name
getwd()
ComDate<- readRDS(paste(getwd(),"/data2sc_hard/",STRING_name,".rds",sep=""))
result <- readRDS(paste(getwd(),"/result_hard/result_",STRING_name,".rds",sep=""))

# all(colnames(ComDate$indata$P) == colnames(ComDate$indata$T))
# 
# MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,
#                   "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",
#                   "MuSiC","Bisque","SCDC", "DWLS",
#                   "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
# result <- CountAllResults(result,MyMethodName)
result$all
nrow(result$all)

# saveRDS(result,paste(getwd(),"/result_hard/result_",STRING_name,".rds",sep=""))
# 
# for (i in MyMethodName) {
#   print(i )
#   print(dim(result[[i]]$p))
# }
# result$CellDistinguisher <- result$CellDistinguisher
# result$CellDistinguisher <- NULL
# result <- CountAllResults(result,MyMethodName)
# result$all
# nrow(result$all)
# STRING_name
# saveRDS(result,paste(getwd(),"/result_hard/result_",STRING_name,".rds",sep=""))

# GetCorMatrix(result$CDSC3$dec$p,ComDate$indata$P)
# diag(GetCorMatrix(result$CDSC3$dec$p,ComDate$indata$P))
# diag(cor(t(result$CIBERSORT$p),t(ComDate$indata$P)))
# diag(cor(t(result$EPIC$p),t(ComDate$indata$P)))
# diag(cor(t(result$ssFrobenius$p),t(ComDate$indata$P)))

#关于C的ct的热图
STRING_name
map = GetCorMatrix(result$CDSC3$dec$c,ComDate$indata$C,ComDate$indata$C_ref,matrix = "c");map
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
eoffice::topptx(plot_pheatmp,
                filename = paste("F:/wangchenqi/CDSC/pictures/class2_hard_NO_ct_",STRING_name,".pptx",sep=""))

para_lambda_Baron_Muraro_123 <- readRDS(
  paste(getwd(),"/paramater/para_lambda_",STRING_name,"_123.rds",sep=""))


# diag(cor(t(result$CDSC3$dec$p),t(ComDate$indata$P)))
# diag(cor(t(result$EPIC$p),t(ComDate$indata$P)))
#-------CDSC---------------------
result <- list()
result$CDSC3$seedd = 44
result$CDSC3$TerCondition = 10^-8
result$CDSC3$Ss <- SM(t(ComDate$indata$T))
result$CDSC3$Sg <- SM(ComDate$indata$T)

result$CDSC3$lambda1 <- 1e-01
result$CDSC3$lambda2 <- 1e+01
result$CDSC3$lambdaC <- 1e+01

library(dplyr)

result$CDSC3$dec = CDSC_3(ComDate$indata$T, ComDate$indata$C_ref, dim(ComDate$indata$C_ref)[2],
                           result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                           result$CDSC3$TerCondition,result$CDSC3$seedd,
                           result$CDSC3$Ss, result$CDSC3$Sg,all_number = 3000)

result$CDSC3$result <- calculate_result(result$CDSC3$dec$c,result$CDSC3$dec$p,
                                         ComDate$indata$T, ComDate$indata$C, ComDate$indata$C_ref, ComDate$indata$P,
                                         result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                                         result$CDSC3$dec$jump, result$CDSC3$seedd, result$CDSC3$TerCondition)

result$CDSC3$result

result$CDSC2$lambda1 <- 1e-03
result$CDSC2$lambda2 <- 1e-02
result$CDSC2$lambdaC <- 0
result$CDSC2$dec = CDSC_2(ComDate$indata$T,  dim(ComDate$indata$C_ref)[2],
                           result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                           result$CDSC3$TerCondition,result$CDSC3$seedd,
                           result$CDSC3$Ss, result$CDSC3$Sg,all_number = 3000)

result$CDSC2$result <- calculate_result(result$CDSC2$dec$c,result$CDSC2$dec$p,
                                         ComDate$indata$T,ComDate$indata$C,ComDate$indata$C_ref, ComDate$indata$P,
                                         result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                                         result$CDSC2$dec$jump, result$CDSC3$seedd, result$CDSC3$TerCondition)
result$CDSC2$result


#----------NNLS----------
require(NNLS)
result$NNLS$p = do.call(cbind.data.frame,lapply(apply(ComDate$indata$T,2,function(x) nnls::nnls(as.matrix(ComDate$indata$C_ref),x)),   function(y) y$x))
result$NNLS$p = apply(result$NNLS$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$NNLS$p) <- colnames(ComDate$indata$C_ref)
result$NNLS$result = getPearsonRMSE(result$NNLS$p,ComDate$indata$P)
result$NNLS$result

#--------OLS------------
result$OLS$p = apply(ComDate$indata$T,2,function(x) lm(x ~ as.matrix(ComDate$indata$C_ref))$coefficients[-1])
result$OLS$p = apply(result$OLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$OLS$p = apply(result$OLS$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$OLS$p) <- unlist(lapply(strsplit(rownames(result$OLS$p),")"),function(x) x[2]))
result$OLS$result = getPearsonRMSE(result$OLS$p,ComDate$indata$P)
result$OLS$result


#---------FARDEEP-----------
library(FARDEEP)
#result_FARDEEP = t(FARDEEP(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = t(FARDEEP::fardeep(ComDate$indata$C_ref, ComDate$indata$T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = apply(result$FARDEEP$p,2,function(x) x/sum(x)) #explicit STO constraint
result$FARDEEP$result = getPearsonRMSE(result$FARDEEP$p,ComDate$indata$P)
result$FARDEEP$result


#----------CIBERSORT-----------------
source("F:/wangchenqi/CDSC/CIBERSORT.R")
result$CIBERSORT$p = CIBERSORT(sig_matrix =ComDate$indata$C_ref, mixture_file = ComDate$indata$T, QN = FALSE)
result$CIBERSORT$p = t(result$CIBERSORT$p[,1:(ncol(result$CIBERSORT$p)-3)])
result$CIBERSORT$result = getPearsonRMSE(result$CIBERSORT$p,ComDate$indata$P)
result$CIBERSORT$result


#----------------"DeconRNASeq------
#nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
#datasets and reference matrix: signatures, need to be non-negative. 
#"use.scale": whether the data should be centered or scaled, default = TRUE
unloadNamespace("DWLS")
unloadNamespace("Seurat") #needed for PCA step
library(pcaMethods) #needed for DeconRNASeq to work
result$deconRNASeq$p = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(ComDate$indata$T), 
                                                  signatures = as.data.frame(ComDate$indata$C_ref), 
                                                  proportions = NULL, 
                                                  checksig = FALSE, 
                                                  known.prop = FALSE, 
                                                  use.scale = TRUE, 
                                                  fig = FALSE)$out.all)
colnames(result$deconRNASeq$p) = colnames(ComDate$indata$T)
require(Seurat)
result$deconRNASeq$result = getPearsonRMSE(result$deconRNASeq$p,ComDate$indata$P)
result$deconRNASeq$result


#-----------RLR----------------
require(MASS)
result$RLR$p = do.call(cbind.data.frame,lapply(apply(ComDate$indata$T,2,function(x) MASS::rlm(x ~ as.matrix(ComDate$indata$C_ref), maxit=100)), function(y) y$coefficients[-1]))
result$RLR$p = apply(result$RLR$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$RLR$p = apply(result$RLR$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$RLR$p) <- unlist(lapply(strsplit(rownames(result$RLR$p),")"),function(x) x[2]))
result$RLR$result = getPearsonRMSE(result$RLR$p,ComDate$indata$P)
result$RLR$result


#-----------DCQ------------------
#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default
require(ComICS)
result$DCQ$p = t(ComICS::dcq(reference_data =ComDate$indata$C_ref, 
                             mix_data = ComDate$indata$T, 
                             marker_set = as.data.frame(row.names(ComDate$indata$C_ref)) ,
                             alpha_used = 0.99, 
                             lambda_min = 0.1, 
                             number_of_repeats = 10)$average)
result$DCQ$p = apply(result$DCQ$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DCQ$p = apply(result$DCQ$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DCQ$result = getPearsonRMSE(result$DCQ$p,ComDate$indata$P)
result$DCQ$result


#-----------elastic_net-----------
#standardize = TRUE by default. lambda=NULL by default 
require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
result$elastic_net$p = apply(ComDate$indata$T, 2, function(z) 
  coef(glmnet::glmnet(x = as.matrix(ComDate$indata$C_ref), 
                      y = z, 
                      alpha = 0.5, 
                      standardize = TRUE, 
                      lambda = glmnet::cv.glmnet(as.matrix(ComDate$indata$C_ref), z)$lambda.1se))[1:ncol(ComDate$indata$C_ref)+1,])
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) x/sum(x)) #explicit STO constraint
result$elastic_net$result = getPearsonRMSE(result$elastic_net$p,ComDate$indata$P)
result$elastic_net$result


#----------ridge----------------
# alpha=0
require(glmnet)
result$ridge$p = apply(ComDate$indata$T, 2, function(z) 
  coef(glmnet::glmnet(x = as.matrix(ComDate$indata$C_ref), 
                      y = z, 
                      alpha = 0.5, 
                      standardize = TRUE, 
                      lambda = glmnet::cv.glmnet(as.matrix(ComDate$indata$C_ref), z)$lambda.1se))[1:ncol(ComDate$indata$C_ref)+1,])
result$ridge$p = apply(result$ridge$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$ridge$p = apply(result$ridge$p,2,function(x) x/sum(x)) #explicit STO constraint
result$ridge$result = getPearsonRMSE(result$ridge$p,ComDate$indata$P)
result$ridge$result


#----------lasso-----------------
#alpha=1; shrinking some coefficients to 0. 
require(glmnet)
result$lasso$p = apply(ComDate$indata$T, 2, function(z) 
  coef(glmnet::glmnet(x = as.matrix(ComDate$indata$C_ref), 
                      y = z, 
                      alpha = 0.5, 
                      standardize = TRUE, 
                      lambda = glmnet::cv.glmnet(as.matrix(ComDate$indata$C_ref), z)$lambda.1se))[1:ncol(ComDate$indata$C_ref)+1,])
result$lasso$p = apply(result$lasso$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$lasso$p = apply(result$lasso$p,2,function(x) x/sum(x)) #explicit STO constraint
result$lasso$p[which(is.na(result$lasso$p) == TRUE)] <- 0 #Needed for models where glmnet drops all terms of a model and fit an intercept-only model (very unlikely but possible).
length(colSums(result$lasso$p) == 1)
result$lasso$result = getPearsonRMSE(result$lasso$p,ComDate$indata$P)
result$lasso$result


# saveRDS(result,paste(getwd(),"/result_hard/result_",STRING_name,".rds",sep=""))

#----------EPIC----------
require(EPIC)
markers = as.character(ComDate$markerslist$gene)
C_EPIC <- list()
common_CTs <- intersect(colnames(ComDate$indata$C_ref),colnames(ComDate$indata$refProfiles.var))
C_EPIC[["sigGenes"]] <- rownames(ComDate$indata$C_ref[markers,common_CTs])
C_EPIC[["refProfiles"]] <- as.matrix(ComDate$indata$C_ref[markers,common_CTs])
C_EPIC[["refProfiles.var"]] <- ComDate$refProfiles.var[markers,common_CTs]

result$EPIC$p <- t(EPIC::EPIC(bulk=as.matrix(ComDate$indata$T), 
                              reference=C_EPIC, 
                              withOtherCells=TRUE, 
                              scaleExprs=T)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices
result$EPIC$result = getPearsonRMSE(result$EPIC$p,ComDate$indata$P)
result$EPIC$result

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")


#--------sc methods_————————————------------

# 
ComDate$sc$T <- ComDate$indata$T
ComDate$sc$C <- ComDate$Train$data
colnames(ComDate$sc$C) <- ComDate$Train$original_cell_names
ComDate$sc$phenoDataC <- ComDate$Train$pData
all(colnames(ComDate$sc$C) == ComDate$sc$phenoDataC$cellID)

ComDate$sc$C <- ComDate$sc$C[intersect(rownames(ComDate$sc$T),rownames(ComDate$sc$C)),]
ComDate$sc$T <- ComDate$sc$T[intersect(rownames(ComDate$sc$T),rownames(ComDate$sc$C)),]


#Bisque requires "SubjectName" in phenoDataC
if(length(grep("[N-n]ame",colnames(ComDate$sc$phenoDataC))) > 0){
  ComDate$sc$sample_column = grep("[N-n]ame",colnames(ComDate$sc$phenoDataC))
} else {
  ComDate$sc$sample_column = grep("[S-s]ample|[S-s]ubject",colnames(ComDate$sc$phenoDataC))
}

colnames(ComDate$sc$phenoDataC)[ComDate$sc$sample_column] = "SubjectName"
rownames(ComDate$sc$phenoDataC) = ComDate$sc$phenoDataC$cellID

require(xbioc)
ComDate$sc$C.eset <- Biobase::ExpressionSet(assayData = as.matrix(ComDate$sc$C)
                                            ,phenoData = Biobase::AnnotatedDataFrame(ComDate$sc$phenoDataC))
ComDate$sc$T.eset <- Biobase::ExpressionSet(assayData = as.matrix(ComDate$sc$T))

#---------MuSiC-----------
require(MuSiC)
source("F:/wangchenqi/CDSC/MuSiC.R")
result$MuSiC$p = t(music_prop_my(bulk.eset = ComDate$sc$T.eset, sc.eset = ComDate$sc$C.eset, 
                                     clusters = 'cellType',
                                     markers = NULL, normalize = FALSE, samples = 'SubjectName', 
                                     verbose = F)$Est.prop.weighted)
result$MuSiC$p <- result$MuSiC$p[rownames(ComDate$indata$P),]
result$MuSiC$result = getPearsonRMSE(result$MuSiC$p, ComDate$indata$P)
result$MuSiC$result

#------------Bisque----------
require(Bisque)
result$Bisque$p <- BisqueRNA::ReferenceBasedDecomposition(ComDate$sc$T.eset, 
                                                          ComDate$sc$C.eset, 
                                                          markers=NULL, 
                                                          use.overlap=FALSE)$bulk.props 
#use.overlap is when there's both bulk and scRNA-seq for the same set of samples
result$Bisque$result = getPearsonRMSE(result$Bisque$p,ComDate$indata$P)
result$Bisque$result


#--------SCDC-----
library(SCDC)
source("F:/wangchenqi/CDSC/SCDC.R")
result$SCDC$p <- t(SCDC::SCDC_prop(bulk.eset = ComDate$sc$T.eset, sc.eset = ComDate$sc$C.eset, 
                                   ct.varname = "cellType", sample = "SubjectName", 
                                   ct.sub = unique(as.character(ComDate$sc$phenoDataC$cellType)), iter.max = 200)$prop.est.mvw)
result$SCDC$result = getPearsonRMSE(result$SCDC$p,ComDate$indata$P)
result$SCDC$result
#---------DWLS--------
require(DWLS)
source('F:/wangchenqi/CDSC/DWLS.R')
getwd()
path=paste(getwd(),"/DWLS/results_hard_",STRING_name,sep=""); path; dir.exists(path)
if(! dir.exists(path)){ #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created

  dir.create(path)
  Signature <- buildSignatureMatrixMAST(scdata = ComDate$sc$C,
                                        id = as.character(ComDate$sc$phenoData$cellType),
                                        path = path, diff.cutoff = 0.5, pval.cutoff = 0.01)

} else {#re-load signature and remove CT column + its correspondent markers

  load(paste(path,"Sig.RData",sep="/"))
  Signature <- Sig

  if(!is.null(elem)){#to be able to deal with full C and with removed CT

    Signature = Signature[,!colnames(Signature) %in% elem]
    CT_to_read <- dir(path) %>% grep(paste(elem,".*RData",sep=""),.,value=TRUE)
    load(paste(path,CT_to_read,sep="/"))

    Signature <- Signature[!rownames(Signature) %in% cluster_lrTest.table$Gene,]
  }

}

result$DWLS$p <- apply(ComDate$sc$T,2, function(x){
  b = setNames(x, rownames(ComDate$sc$T))
  tr <- trimData(Signature, b)
  RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
})

rownames(result$DWLS$p) <- as.character(unique(ComDate$sc$phenoDataC$cellType))
result$DWLS$p = apply(result$DWLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DWLS$p = apply(result$DWLS$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DWLS$result = getPearsonRMSE(result$DWLS$p,ComDate$indata$P)
result$DWLS$result

# result <- CountAllResults(result)
# result$all
# saveRDS(result, paste(getwd(),"/result/result_",STRING_name,".rds",sep=""))
# saveRDS(result, paste(getwd(),"/result_hard/result_",STRING_name,".rds",sep=""))

#-----complete deconvolution methods————----------
#--------DSA-------------------
require(CellMix)
md = ComDate$markerslist
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

result$DSA = CellMix::ged(as.matrix(ComDate$indata$T), ML, method = "DSA", log = TRUE)
result$DSA$c = result$DSA@fit@W
result$DSA$p = result$DSA@fit@H
result$DSA$t = result$DSA$c%*%result$DSA$p

result$DSA$result$c = getPearsonRMSE(result$DSA$c,ComDate$indata$C)
result$DSA$result$p = getPearsonRMSE(result$DSA$p,ComDate$indata$P)
result$DSA$result$t = getPearsonRMSE(result$DSA$t,ComDate$indata$T)

result$DSA$result$all <- cbind(result$DSA$result$p,result$DSA$result$c)
result$DSA$result$all <-  cbind(result$DSA$result$all, result$DSA$result$t)
colnames(result$DSA$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$DSA$result <- result$DSA$result$all
result$DSA$result

#---------ssKL------------------
require(CellMix)
md = ComDate$markerslist
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

result$ssKL <- CellMix::ged(as.matrix(ComDate$indata$T), ML, method = "ssKL", sscale = FALSE, maxIter=500, log = FALSE)
result$ssKL$c = result$ssKL@fit@W
result$ssKL$p = result$ssKL@fit@H
result$ssKL$t = result$ssKL$c%*%result$ssKL$p

result$ssKL$result$c = getPearsonRMSE(result$ssKL$c,ComDate$indata$C)
result$ssKL$result$p = getPearsonRMSE(result$ssKL$p,ComDate$indata$P)
result$ssKL$result$t = getPearsonRMSE(result$ssKL$t,ComDate$indata$T)

result$ssKL$result$all <- cbind(result$ssKL$result$p,result$ssKL$result$c)
result$ssKL$result$all <-  cbind(result$ssKL$result$all,result$ssKL$result$t)
colnames(result$ssKL$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssKL$result <- result$ssKL$result$all
result$ssKL$result

#----------ssFrobenius-----------------
require(CellMix)
md = ComDate$markerslist
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

result$ssFrobenius <- CellMix::ged(as.matrix(ComDate$indata$T), ML, method = "ssFrobenius", sscale = TRUE, maxIter = 500, log = FALSE) #equivalent to coef(CellMix::ged(T,...)
result$ssFrobenius$c = result$ssFrobenius@fit@W
result$ssFrobenius$p = result$ssFrobenius@fit@H
result$ssFrobenius$t = result$ssFrobenius$c%*%result$ssFrobenius$p

result$ssFrobenius$result$c = getPearsonRMSE(result$ssFrobenius$c,ComDate$indata$C)
result$ssFrobenius$result$p = getPearsonRMSE(result$ssFrobenius$p,ComDate$indata$P)
result$ssFrobenius$result$t = getPearsonRMSE(result$ssFrobenius$t,ComDate$indata$T)

result$ssFrobenius$result$all <- cbind(result$ssFrobenius$result$p,result$ssFrobenius$result$c)
result$ssFrobenius$result$all <-  cbind(result$ssFrobenius$result$all,result$ssFrobenius$result$t)
colnames(result$ssFrobenius$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssFrobenius$result <- result$ssFrobenius$result$all
result$ssFrobenius$result

#-----------deconf------------
require(CellMix)
md = ComDate$markerslist
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

all(rownames(ComDate$indata$T)==ComDate$markers1_shift$gene)
all(rownames(ComDate$indata$T)==rownames(ComDate$shif_marker))

result$deconf <- CellMix::ged(as.matrix(ComDate$indata$T), ncol(ComDate$indata$C_ref))
result$deconf$c = result$deconf@fit@W
result$deconf$p = result$deconf@fit@H
result$deconf$t = result$deconf$c%*%result$deconf$p
ctlabels <- Row_label(result$deconf$c,ComDate$indata$C_ref,leastnum=3);ctlabels
colnames(result$deconf$c) <- ctlabels
rownames(result$deconf$p) <- ctlabels
# result$deconf <- CellMix::ged(as.matrix(ComDate$indata$T), ML, method = "deconf")
# result$deconf$c = result$deconf@fit@W
# result$deconf$p = result$deconf@fit@H
# result$deconf$t = result$deconf$c%*%result$deconf$p

result$deconf$result$c = getPearsonRMSE(result$deconf$c,ComDate$indata$C)
result$deconf$result$p = getPearsonRMSE(result$deconf$p,ComDate$indata$P)
result$deconf$result$t = getPearsonRMSE(result$deconf$t,ComDate$indata$T)

result$deconf$result$all <- cbind(result$deconf$result$p,result$deconf$result$c)
result$deconf$result$all <-  cbind(result$deconf$result$all,result$deconf$result$t)
colnames(result$deconf$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$deconf$result <- result$deconf$result$all
result$deconf$result

#--------CDSeq-------------------
# library(CDSeq)
# star_time <- Sys.time()
# result$CDSeq.noRef <- CDSeq(bulk_data =  as.matrix(ComDate$indata$T), 
#                             cell_type_number = dim(ComDate$indata$C_ref)[2], 
#                             mcmc_iterations = 100, 
#                             cpu_number=20)
# end_time <- Sys.time()
# end_time-star_time
# result$CDSeq.noRef$result <- calculate_result(result$CDSeq.noRef$estGEP,result$CDSeq.noRef$estProp,
#                                               ComDate$indata$T,ComDate$indata$C,ComDate$indata$C_ref, ComDate$indata$P)
# result$CDSeq.noRef$time <- end_time-star_time
# result$CDSeq.noRef$result
# 
# star_time <- Sys.time()
# result$CDSeq.haveRef<-CDSeq(bulk_data =  as.matrix(ComDate$indata$T), 
#                             cell_type_number = dim(ComDate$indata$C_ref)[2], 
#                             mcmc_iterations = 20, 
#                             # gene_length = as.vector(gene_length), 
#                             reference_gep = as.matrix(ComDate$indata$C_ref),  # gene expression profile of pure cell lines
#                             cpu_number = 8)
# end_time <- Sys.time()
# end_time-star_time
# result$CDSeq.haveRef$result <- calculate_result(result$CDSeq.haveRef$estGEP,result$CDSeq.haveRef$estProp,
#                                                 ComDate$indata$T,ComDate$indata$C,ComDate$indata$C_ref, ComDate$indata$P)
# result$CDSeq.haveRef$time <- end_time-star_time
# result$CDSeq.haveRef$result
#------TOAST + NMF-------
# TOAST方法宝函数只能输入基因大于样本的数据。。
require(DeCompress)
source("F:/wangchenqi/CDSC/DeCompress.R")# csDeCompress_my

dim(ComDate$indata$T)
result$TOAST$toast.nmf <- DeCompress::csDeCompress(Y_raw = ComDate$indata$T,
                                                   K = dim(ComDate$indata$C_ref)[2],
                                                   nMarker = nrow(ComDate$indata$T),
                                                   FUN = nmfOut,
                                                   TotalIter = 30)
require(NMF)
result$TOAST$fin.nmf = nmf(x = ComDate$indata$T,
                           rank = dim(ComDate$indata$C_ref)[2])
result$TOAST$ppp = t(coef(result$TOAST$fin.nmf))
result$TOAST$ppp = t(apply(result$TOAST$ppp,1,function(c) c/sum(c)))
result$TOAST$p <- t(result$TOAST$ppp)
result$TOAST$c <- basis(result$TOAST$fin.nmf)
result$TOAST$t <- result$TOAST$c%*%result$TOAST$p
nmf.res = list(prop = result$TOAST$ppp,
               sig = basis(result$TOAST$fin.nmf))
all(rownames(result$TOAST$c) == rownames(ComDate$indata$C_ref))
# result$Linseed$c <- result$CellDistinguisher$c[rownames(ComDate$indata$C_ref),]


rownames(ComDate$indata$P)
rownames(result$TOAST$p)
result$TOAST$result <- calculate_result(result$TOAST$c,result$TOAST$p,
                   ComDate$indata$T, ComDate$indata$C, ComDate$indata$C_ref, ComDate$indata$P)
result$TOAST$result

#-------Linseed---------------
require(DeCompress)
result$Linseed$Linseed.rs = DeCompress::linCor(yref = ComDate$indata$T,
                                               iters = 100,
                                               pval = 100,
                                               n.types = dim(ComDate$indata$C_ref)[2],
                                               scree = 'drop',
                                               logTransform = F)

# names(result$Linseed$Linseed.rs) = names(result$Linseed$nmf.res)
result$Linseed$c <- result$Linseed$Linseed.rs[[2]]
result$Linseed$p <- result$Linseed$Linseed.rs[[1]]
# result$Linseed$t <- result$Linseed$c%*%result$Linseed$p

all(rownames(result$Linseed$c) == rownames(ComDate$indata$C_ref))
all(rownames(result$Linseed$c) == rownames(ComDate$indata$T))

result$Linseed$c <- result$Linseed$c[rownames(ComDate$indata$C_ref),]
# result$Linseed$t <- result$Linseed$t[rownames(ComDate$indata$C_ref),]


rownames(ComDate$indata$P)
rownames(result$Linseed$p)

result$Linseed$result <- calculate_result(result$Linseed$c,result$Linseed$p,
                                        ComDate$indata$T, ComDate$indata$C, ComDate$indata$C_ref, ComDate$indata$P)
result$Linseed$result

#-----CellDistinguisher-----
require(DeCompress)
result$CellDistinguisher <- CellDistinguisher::gecd_CellDistinguisher(
  ComDate$indata$T,
  genesymb=rownames(ComDate$indata$T),
  numCellClasses = dim(ComDate$indata$C_ref)[2],
  minDistinguisherAlternatives=1,
  maxDistinguisherAlternatives=100,
  minAlternativesLengthsNormalized=0.5,
  expressionQuantileForScale = 0.75,
  expressionQuantileForFilter=0.999,
  expressionConcentrationRatio=0.333,
  probesWithGenesOnly = F,
  verbose=1)
result$CellDist.deconv <-
  tryCatch(CellDistinguisher::gecd_DeconvolutionByDistinguishers(
    as.matrix(ComDate$indata$T),
    result$CellDistinguisher$bestDistinguishers,
    nonNegativeOnly = T,
    convexSolution = T,
    verbose = 0),
    error = function(e) return(list(sampleCompositions =
                                      matrix(rep(1/dim(ComDate$indata$C_ref)[2],
                                                 dim(ComDate$indata$C_ref)[2]*ncol(ComDate$indata$T)),
                                             ncol=dim(ComDate$indata$C_ref)[2]))))
if (length(result$CellDist.deconv) > 1){
  result$CellDistinguisher$p <- result$CellDist.deconv$sampleComposition
  result$CellDistinguisher$c <- result$CellDist.deconv$cellSubclassSignatures
}else {
  result$CellDistinguisher$p <- t(nmf.res$prop)
  result$CellDistinguisher$c <- nmf.res$sig
}
all(rownames(result$CellDistinguisher$c) == rownames(ComDate$indata$C_ref))
result$Linseed$c <- result$CellDistinguisher$c[rownames(ComDate$indata$C_ref),]

# result$CellDistinguisher$t <- result$CellDistinguisher$c%*%result$CellDistinguisher$p
rownames(ComDate$indata$P)
rownames(result$CellDistinguisher$p)
result$CellDistinguisher$result <- calculate_result(result$CellDistinguisher$c,result$CellDistinguisher$p,
                                          ComDate$indata$T, ComDate$indata$C, ComDate$indata$C_ref, ComDate$indata$P)
result$CellDistinguisher$result
#-------------DeCompress--------
# result$DeCompress$decompress.res = bestDeconvolution( ComDate$indata$T,
# n.types = dim(ComDate$indata$C_ref)[2],
# scree = 'cumvar',
# logTransform = F,
# known.props = NULL,
# methods = c('TOAST',
#             'Linseed',
#             'CellDistinguisher'))

#--------------all————————-----------------
MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,
                  "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",
                  "MuSiC","Bisque","SCDC", "DWLS",
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
result <- CountAllResults(result,MyMethodName)
result$all;dim(result$all)
STRING_name

saveRDS(result,paste(getwd(),"/result_hard/result_",STRING_name,".rds",sep=""))

#----------find Paramater ///----------
lambda1 <- c(0,10^-3,10^-2,10^-1,1,10,100)
lambda2 <- c(0,10^-3,10^-2,10^-1,1,10,100)
lambdaC <- c(0,10^-1,10^0,10^1,10^2,10^3,10^4,10^5)

#------------cyclic to do deconvolution----------
pb <- txtProgressBar(style = 3)
star_time <- Sys.time()
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
      result_CDSC = CDSC_3(ComDate$indata$T,ComDate$indata$C_ref, dim(ComDate$indata$C_ref)[2], 
                           lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k],
                           result$CDSC3$TerCondition,result$CDSC3$seedd,
                           result$CDSC3$Ss,result$CDSC3$Sg,all_number = 3000)
      result_para_c[[num]] <- result_CDSC[[1]]
      result_para_p[[num]] <- result_CDSC[[2]]
      number_iter[num] <- result_CDSC[[3]]
      
      pearson_para_c[[num]] <- cor(result_para_c[[num]],ComDate$indata$C_ref)
      result1 = NULL
      if(!all(is.na(pearson_para_c[[num]]) == FALSE) ){
        break
      }
      result1 <- calculate_result(result_para_c[[num]],result_para_p[[num]],
                                  ComDate$indata$T,ComDate$indata$C,ComDate$indata$C_ref,ComDate$indata$P,
                                  lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k],number_iter[num],
                                  result$CDSC3$seedd,result$CDSC3$TerCondition)
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
# para_lambda_44_8_12 <- para_lambda_44_8

saveRDS(para_lambda_44_8,
        paste(getwd(),"/paramater_hard/para_lambda_",STRING_name,".rds",sep=""))


# 热图
rm (list=ls ())
setwd("F:/wangchenqi/CDSC/2scDeconvolution")

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/2scDeconvolution/data2sc_hard"));a
a = c("Baron_Muraro","Muraro_Baron","Baron_Segerstolpe","Segerstolpe_Baron",
      "Segerstolpe_Muraro","Muraro_Segerstolpe","Macosko_Shekhar","Shekhar_Macosko" )
ComDate = NULL
result = NULL
for(i in 1:length(a)){
  STRING_name = a[i]; STRING_name
  # ComData[[i]] <- readRDS(paste(getwd(),"/data2sc/",STRING_name,".rds",sep=""))
  result[[i]] <- readRDS(paste(getwd(),"/result_hard/result_",STRING_name,".rds",sep=""))
}
MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,
                  "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",
                  "MuSiC","Bisque","SCDC", "DWLS",
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
names(result[[1]])
map = NULL
for(i in 1:length(a)){
  result[[i]] <- CountAllResults(result[[i]],MyMethodName)
  result[[i]]$all;dim(result[[i]]$all)
  if(nrow(result[[i]]$all)!=length(MyMethodName)){
    print("ERROR")
    break;}
  map <- cbind(map,result[[i]]$all$RMSE_to_C)
}
rownames(map) <- MyMethodName;colnames(map) <- a
map;dim(map)

library(pheatmap)#c("navy", "white", "firebrick3") color
plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
             ,cellwidth = 42, cellheight = 12,gaps_row = c(12,16)
             ,breaks = c(seq(0.4,1,by = 0.02))
             ,display_numbers = TRUE, number_format = "%.3f"
             ,border_color = "black",color = colorRampPalette(c( "white", "firebrick3"))(30)
             ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
             ,main = "Pearson to P")

plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 42, cellheight = 12,gaps_row =c(12,16)
                         ,display_numbers = TRUE, number_format = "%.3f"
                         ,border_color = "black",color = colorRampPalette(c( "white", "navy"))(70)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "RMSE to P")

plot_pheatmp <- pheatmap(map[c(1,17:24),],cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 30, cellheight = 12#,gaps_row = c(12, 17)
                         ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.3f"
                         ,border_color = "black",color = colorRampPalette(c("white", "firebrick3"))(30)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "Pearson to C")
plot_pheatmp <- pheatmap(map[c(1,17:24),],cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 30, cellheight = 12#,gaps_row = c(12, 17)
                         # ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.0f"
                         ,border_color = "black",color = colorRampPalette(c("white", "navy"))(30)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "RMSE to C")
apply(map, 1, mean)
plot_pheatmp
library(ggplot2)
eoffice::topptx(plot_pheatmp,filename = 
                  "F:/wangchenqi/CDSC/pictures/class2_hard_NO_RMSE_of_C.pptx")

# 泛化能力的箱线图--p----
mapBox <- NULL
mapBox =  data.frame(pearson = map[1,], method = rownames(map)[1],row.names = NULL)
nn1 <- nrow(mapBox)
for (i in 1:length(rownames(map))) {
  mapBox <- rbind(mapBox, data.frame(pearson = map[i,], method = rownames(map)[i],row.names = NULL))
}

mapBox <- mapBox[-(1:nn1),]
library(ggplot2)#coef=1e30,
plot_class1 =  ggplot(mapBox, aes(factor(method,levels=MyMethodName),pearson,fill=method)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_class1;
eoffice::topptx(plot_class1,
                filename = "F:/wangchenqi/CDSC/pictures/class2_hard_NO_fanhua_p.pptx")

# 泛化能力的箱线图----c-----
mapBox <- NULL
methodsNames2 <-  c("CDSC3","CDSC2", "DSA","ssKL" ,"ssFrobenius",
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
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_class1;
eoffice::topptx(plot_class1,filename = 
                  "F:/wangchenqi/CDSC/pictures/class2_hard_NO_fanhua_c.pptx")

# 柱状图--p
rm (list=ls ())
setwd("F:/wangchenqi/CDSC/2scDeconvolution")

source("F:/wangchenqi/CDSC/CDSC.R")
source("F:/wangchenqi/CDSC/CDSC_expand.R")

a = gsub('.rds','',list.files("F:/wangchenqi/CDSC/2scDeconvolution/data2sc_hard"));a
a = c("Baron_Muraro","Muraro_Baron","Baron_Segerstolpe","Segerstolpe_Baron",
      "Segerstolpe_Muraro","Muraro_Segerstolpe","Macosko_Shekhar","Shekhar_Macosko" )

ComDate = NULL
result = NULL

for(i in 1:length(a)){
  STRING_name = a[i]; STRING_name
  ComDate[[i]] <- readRDS(paste(getwd(),"/data2sc_hard/",STRING_name,".rds",sep=""))
  result[[i]] <- readRDS(paste(getwd(),"/result_hard/result_",STRING_name,".rds",sep=""))
}
rownames(result[[i]]$CDSC3$dec$p)

for (i in 1:length(a)) {
  ctlabels <- Row_label(result[[i]]$CDSC3$dec$c,ComDate[[i]]$indata$C_ref,leastnum=3);ctlabels
  rownames(result[[i]]$CDSC3$dec$p) <- ctlabels
  colnames(result[[i]]$CDSC3$dec$c) <- ctlabels
  ctlabels_ <- intersect(rownames(ComDate[[i]]$indata$P),as.character(ctlabels));ctlabels_
  result[[i]]$CDSC3$p <- result[[i]]$CDSC3$dec$p[ctlabels_, ]
  result[[i]]$CDSC3$c <- result[[i]]$CDSC3$dec$c[ ,ctlabels_]
  ComDate[[i]]$indata$C = as.matrix(ComDate[[i]]$indata$C[rownames(result[[i]]$CDSC3$c),colnames(result[[i]]$CDSC3$c)])
  ComDate[[i]]$indata$P = as.matrix(ComDate[[i]]$indata$P[rownames(result[[i]]$CDSC3$p),colnames(result[[i]]$CDSC3$p)])
}

map = NULL
for(i in 1:length(a)){
  map <- rbind(map,data.frame(group = a[i],
                              pearson = diag(cor(result[[i]]$CDSC3$p,ComDate[[i]]$indata$P))))
}
map = map[which(is.na(map$pearson) == F),]
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
  scale_fill_manual(values = c("#80AFBF","#80AFBF","#80AFBF","#80AFBF","#80AFBF","#80AFBF","#80AFBF","#80AFBF"))+ # 柱状图颜色, "#DA635D","#B1938B"
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

  labs(title="",x="METHODS",y="Pearson of samples in P")+ # 添加标题，x轴，y轴内容
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
eoffice::topptx(p1,
                filename = "F:/wangchenqi/CDSC/pictures/class2_hard_NO_bar_samples.pptx")


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
eoffice::topptx(plot_class2,
                filename = "F:/wangchenqi/CDSC/pictures/class2_hard_NO_boxplot_samples.pptx")
