rm (list=ls ())

source("CDSC.R")
source("function_help.R")

#------data preparation----
Segerstolpe <- list(data = readRDS("Segerstolpe.rds"), full_phenoData = readRDS("Segerstolpe_phenoData.rds"))
sort(table(Segerstolpe$full_phenoData$cellType))
# Segerstolpe$simualte <- scSimulate(Segerstolpe, leastNum=50, plotmarker = F)
Segerstolpe$simulate1 <- scSimulateSplit(Segerstolpe,
                                         leastNum=50, plotmarker = F,
                                         norm1 = "CPM",log2.threshold = 1)
nrow(Segerstolpe$simulate1$markerslist)
table(Segerstolpe$simulate1$markerslist$CT)
Segerstolpe <- scSimulateShift(Segerstolpe,"all",standardization=TRUE)
Segerstolpe$simulate1$T["REG1A",1:5]
Segerstolpe$Indata$T["REG1A",1:5]
saveRDS(Segerstolpe,"XXX.rds")
#--------sc------------

a = gsub('.rds','',list.files("1scSimulate/dataSimulate"));a
a = c("Nestorowa","Manno","Darmanis","Camp","Segerstolpe")
STRING_name = a[3]; STRING_name
getwd()
scData <- readRDS(paste(getwd(),"/dataSimulate/",STRING_name,".rds",sep=""))
result <- readRDS(paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))
# result0 <- result

result$all;dim(result$all)

MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,
                  "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",
                  "MuSiC","Bisque","SCDC", "DWLS",
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
result <- CountAllResults(result,MyMethodName)
result$all;dim(result$all)
# paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep="")
# saveRDS(result,paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))

# paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep="")
# saveRDS(result,paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))

# GetCorMatrix(result$CDSC3$dec$p,scData$Indata$P,matrix = "p")
# diag(GetCorMatrix(result$CDSC3$dec$p,scData$Indata$P,matrix = "p"))
# diag(cor(t(result$EPIC$p),t(scData$Indata$P)))
# diag(cor(t(result$MuSiC$p),t(scData$Indata$P)))
# diag(cor(t(result$linseed$p),t(scData$Indata$P)))

#pheatMap of C's ct
library(pheatmap)
STRING_name
map = GetCorMatrix(result$CDSC3$dec$c,scData$Indata$C,matrix = "c");map
plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 15, cellheight =15#,gaps_row = c(12, 17)
                         ,fontsize = 6
                         # ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.2f" # 鏁板瓧
                         # ,border_color = "black",color = colorRampPalette(MyColor)(60)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,annotation_legend = T
                         ,main = STRING_name)

eoffice::topptx(plot_pheatmp,filename = 
                  paste("pictures/class1_ct_NO_",
                        STRING_name,".pptx",sep=""))



dim(scData$data)
dim(scData$Indata$C)
result$all
# 
# paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep="")
# saveRDS(result,paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))

para_lambda_123 <- readRDS(paste(getwd(),"/paramater/para_lambda_",STRING_name,"_123.rds",sep=""))
para_lambda_12 <- readRDS(paste(getwd(),"/paramater/para_lambda_",STRING_name,"_12.rds",sep=""))
result$all
result$CDSC3$result
result$CDSC2$result

#-------CDSC-------------
result <- list()
result$CDSC3$seedd = 44
result$CDSC3$TerCondition = 10^-8
result$CDSC3$Ss <- SM(t(scData$Indata$T))
result$CDSC3$Sg <- SM(scData$Indata$T)
# 
result$CDSC3$lambda1 <- 1e-03
result$CDSC3$lambda2 <- 0e+00
result$CDSC3$lambdaC <- 10000

library(dplyr)



# result
result$CDSC3$dec = CDSC(scData$Indata$T, 
                          scData$Indata$C_ref, 
                          dim(scData$Indata$C_ref)[2], 
                          result$CDSC3$lambda1, 
                          result$CDSC3$lambda2, 
                          result$CDSC3$lambdaC,
                          result$CDSC3$TerCondition,
                          result$CDSC3$seedd,
                          result$CDSC3$Ss,
                          result$CDSC3$Sg,
                          all_number = 3000)

result$CDSC3$result <- calculate_result(result$CDSC3$dec$c,result$CDSC3$dec$p,
                            scData$Indata$T,scData$Indata$C,scData$Indata$C_ref,scData$Indata$P,
                            result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                            number_iter = result$CDSC3$dec$jump,
                            seedd = result$CDSC3$seedd,
                            TerCondition = result$CDSC3$TerCondition)
result$CDSC3$result

result$CDSC2$lambda1 <- 1e-02
result$CDSC2$lambda2 <- 1e+01
result$CDSC2$lambdaC <- 0
result$CDSC2$dec = CDSC(scData$Indata$T,  NULL, dim(scData$Indata$C_ref)[2], 
                           result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                           result$CDSC3$TerCondition,result$CDSC3$seedd,
                           result$CDSC3$Ss,result$CDSC3$Sg,all_number = 3000)

result$CDSC2$result <- calculate_result(result$CDSC2$dec$c,result$CDSC2$dec$p,
                                         scData$Indata$T, scData$Indata$C, scData$Indata$C_ref, scData$Indata$P,
                                         result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                                        number_iter = result$CDSC2$dec$jump,
                                        seedd = result$CDSC3$seedd,
                                        TerCondition = result$CDSC3$TerCondition)
result$CDSC2$result
# scData = Camp
#----------NNLS----------
require(nnls)
result$NNLS$p = do.call(cbind.data.frame,lapply(apply(scData$Indata$T,2,function(x) nnls::nnls(as.matrix(scData$Indata$C_ref),x)),   function(y) y$x))
result$NNLS$p = apply(result$NNLS$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$NNLS$p) <- colnames(scData$Indata$C_ref)
result$NNLS$result = getPearsonRMSE(result$NNLS$p,scData$Indata$P)
result$NNLS$result
diag(cor(t(result$NNLS$p),t(scData$Indata$P)))

#--------OLS------------
result$OLS$p = apply(scData$Indata$T,2,function(x) lm(x ~ as.matrix(scData$Indata$C_ref))$coefficients[-1])
result$OLS$p = apply(result$OLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$OLS$p = apply(result$OLS$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$OLS$p) <- unlist(lapply(strsplit(rownames(result$OLS$p),")"),function(x) x[2]))
result$OLS$result = getPearsonRMSE(result$OLS$p,scData$Indata$P)
result$OLS$result


#---------FARDEEP-----------
library(FARDEEP)
#result_FARDEEP = t(FARDEEP(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = t(FARDEEP::fardeep(scData$Indata$C_ref, scData$Indata$T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = apply(result$FARDEEP$p,2,function(x) x/sum(x)) #explicit STO constraint
result$FARDEEP$result = getPearsonRMSE(result$FARDEEP$p,scData$Indata$P)
result$FARDEEP$result


#----------CIBERSORT-----------------
source("CIBERSORT.R")
result$CIBERSORT$p = CIBERSORT(sig_matrix = scData$Indata$C_ref, mixture_file = scData$Indata$T, QN = FALSE)
result$CIBERSORT$p = t(result$CIBERSORT$p[,1:(ncol(result$CIBERSORT$p)-3)])
result$CIBERSORT$result = getPearsonRMSE(result$CIBERSORT$p,scData$Indata$P)
result$CIBERSORT$result
diag(cor(t(result$CIBERSORT$p),t(scData$Indata$P)))

#----------------DeconRNASeq------
#nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
#datasets and reference matrix: signatures, need to be non-negative. 
#"use.scale": whether the data should be centered or scaled, default = TRUE
unloadNamespace("Seurat") #needed for PCA step
library(pcaMethods) #needed for DeconRNASeq to work
result$deconRNASeq$p = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(scData$Indata$T), 
                                                  signatures = as.data.frame(scData$Indata$C_ref), 
                                                  proportions = NULL, checksig = FALSE, 
                                                  known.prop = FALSE, use.scale = TRUE,
                                                  fig = FALSE)$out.all)
colnames(result$deconRNASeq$p) = colnames(scData$Indata$T)
result$deconRNASeq$p[,1:10]
require(Seurat)
result$deconRNASeq$result = getPearsonRMSE(result$deconRNASeq$p,scData$Indata$P)
result$deconRNASeq$result


#-----------RLR----------------
require(MASS)
result$RLR$p = do.call(cbind.data.frame,lapply(apply(scData$Indata$T,2,function(x) MASS::rlm(x ~ as.matrix(scData$Indata$C_ref), maxit=100)), function(y) y$coefficients[-1]))
result$RLR$p = apply(result$RLR$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$RLR$p = apply(result$RLR$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$RLR$p) <- unlist(lapply(strsplit(rownames(result$RLR$p),")"),function(x) x[2]))
result$RLR$result = getPearsonRMSE(result$RLR$p,scData$Indata$P)
result$RLR$result


#-----------DCQ------------------
#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default
require(ComICS)
result$DCQ$p = t(ComICS::dcq(reference_data = scData$Indata$C_ref, 
                             mix_data = scData$Indata$T, 
                             marker_set = as.data.frame(row.names(scData$Indata$C_ref)) , 
                             alpha_used = 0.99, 
                             lambda_min = 0.1, 
                             number_of_repeats = 10)$average)
result$DCQ$p = apply(result$DCQ$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DCQ$p = apply(result$DCQ$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DCQ$result = getPearsonRMSE(result$DCQ$p,scData$Indata$P)
result$DCQ$result
# cor(t(result$DCQ$p),t(scData$Indata$P))

#-----------elastic_net-----------
#standardize = TRUE by default. lambda=NULL by default 
require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
result$elastic_net$p = apply(scData$Indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(scData$Indata$C_ref), y = z, alpha = 0.2, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(scData$Indata$C_ref), z)$lambda.1se))[1:ncol(scData$Indata$C_ref)+1,])
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) x/sum(x)) #explicit STO constraint
result$elastic_net$result = getPearsonRMSE(result$elastic_net$p,scData$Indata$P)
result$elastic_net$result


#----------ridge----------------
# alpha=0
require(glmnet)
result$ridge$p = apply(scData$Indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(scData$Indata$C_ref), y = z, alpha = 0, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(scData$Indata$C_ref), z)$lambda.1se))[1:ncol(scData$Indata$C_ref)+1,])
result$ridge$p = apply(result$ridge$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$ridge$p = apply(result$ridge$p,2,function(x) x/sum(x)) #explicit STO constraint
result$ridge$result = getPearsonRMSE(result$ridge$p,scData$Indata$P)
result$ridge$result


#----------lasso-----------------
#alpha=1; shrinking some coefficients to 0. 
require(glmnet)
result$lasso$p = apply(scData$Indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(scData$Indata$C_ref), y = z, alpha = 1, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(scData$Indata$C_ref), z)$lambda.1se))[1:ncol(scData$Indata$C_ref)+1,])
result$lasso$p = apply(result$lasso$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$lasso$p = apply(result$lasso$p,2,function(x) x/sum(x)) #explicit STO constraint
result$lasso$p[which(is.na(result$lasso$p) == TRUE)] <- 0 #Needed for models where glmnet drops all terms of a model and fit an intercept-only model (very unlikely but possible).
length(colSums(result$lasso$p) == 1)
result$lasso$result = getPearsonRMSE(result$lasso$p,scData$Indata$P)
result$lasso$result


# saveRDS(result,paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))

#
# source("CDSC.R")
# source("CDSC_expand.R")

#----------EPIC----------
require(EPIC)
C_EPIC <- list()
C_EPIC$marker = scData$simulate1$markerslist[scData$simulate1$markerslist$gene %in% rownames(scData$Indata$C_ref),]
C_EPIC$markers = as.character(C_EPIC$marker$gene)
C_EPIC$CTs <- intersect(colnames(scData$Indata$C_ref),colnames(scData$Indata$refProfiles.var))
C_EPIC[["sigGenes"]] <- rownames(scData$Indata$C_ref[C_EPIC$markers,C_EPIC$CTs])
C_EPIC[["refProfiles"]] <- as.matrix(scData$Indata$C_ref[C_EPIC$markers,C_EPIC$CTs])
C_EPIC[["refProfiles.var"]] <- scData$Indata$refProfiles.var[C_EPIC$markers,C_EPIC$CTs]

result$EPIC$p <- t(EPIC::EPIC(bulk=as.matrix(scData$Indata$T), 
                              reference=C_EPIC, withOtherCells=TRUE, 
                              scaleExprs=TRUE)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices
result$EPIC$p = result$EPIC$p[!rownames(result$EPIC$p) %in% "otherCells",]
result$EPIC$result = getPearsonRMSE(result$EPIC$p,scData$Indata$P)
result$EPIC$result

#--------sc methods------------
scData$sc$keep <- which(rownames(scData$simulate1$split$train) %in% scData$simulate1$shif_marker$gene)
# scData$sc$keep <- which(rownames(scData$simulate1$split$train) %in% rownames(scData$simulate1$split$train))
scData$sc$T <- scData$simulate1$T[scData$sc$keep,]

scData$sc$C <- scData$simulate1$split$train[scData$sc$keep,]
colnames(scData$sc$C) <- scData$simulate1$original_cell_names[scData$simulate1$split$training]
scData$sc$C <- scData$sc$C[rownames(scData$sc$T),]
scData$sc$phenoDataC <- scData$simulate1$split$pData_train

#Bisque requires "SubjectName" in phenoDataC
if(length(grep("[N-n]ame",colnames(scData$sc$phenoDataC))) > 0){
  scData$sc$sample_column = grep("[N-n]ame",colnames(scData$sc$phenoDataC))
} else {
  scData$sc$sample_column = grep("[S-s]ample|[S-s]ubject",colnames(scData$sc$phenoDataC))
}

colnames(scData$sc$phenoDataC)[scData$sc$sample_column] = "SubjectName"
rownames(scData$sc$phenoDataC) = scData$sc$phenoDataC$cellID

require(xbioc)
colnames(scData$sc$C) <- scData$simulate1$original_cell_names[scData$simulate1$split$training]
scData$sc$C.eset <- Biobase::ExpressionSet(assayData = as.matrix(scData$sc$C)
                                           ,phenoData = Biobase::AnnotatedDataFrame(scData$sc$phenoDataC))
scData$sc$T.eset <- Biobase::ExpressionSet(assayData = as.matrix(scData$sc$T))

#---------MuSiC-----------
library(MuSiC)
source("MuSiC.R")
result$MuSiC$p = t(music_prop_my(bulk.eset = scData$sc$T.eset, 
                                     sc.eset = scData$sc$C.eset,
                                     markers = NULL,
                                     clusters = "cellType",
                                     samples = 'SubjectName',
                                     select.ct = unique(as.character(scData$sc$phenoDataC$cellType)),
                                     verbose = F,
                                     normalize = F
                                      
                                     )$Est.prop.weighted)
result$MuSiC$p[,1:10]
scData$Indata$P[,1:10]
result$MuSiC$result = getPearsonRMSE(result$MuSiC$p,scData$Indata$P)
result$MuSiC$result

#------------Bisque----------
require(Bisque)
result$Bisque$p <- BisqueRNA::ReferenceBasedDecomposition(scData$sc$T.eset, scData$sc$C.eset, 
                                                             markers=NULL, use.overlap=FALSE)$bulk.props 
#use.overlap is when there's both bulk and scRNA-seq for the same set of samples
result$Bisque$result = getPearsonRMSE(result$Bisque$p,scData$Indata$P)
result$Bisque$result

#--------SCDC-----
library(SCDC)
source("SCDC.R")
# result$SCDC$p <- t(SCDC::SCDC_prop(bulk.eset = scData$sc$T.eset,
result$SCDC$p <- t(SCDC_prop_my(bulk.eset = scData$sc$T.eset,
                                sc.eset = scData$sc$C.eset, 
                                ct.varname = "cellType", 
                                sample = "SubjectName", 
                                ct.sub = unique(as.character(scData$sc$phenoDataC$cellType)), 
                                iter.max = 1000)$prop.est.mvw)
result$SCDC$p[,1:5]
scData$Indata$P[,1:5]
result$SCDC$p <- result$SCDC$p[rownames(scData$Indata$P),]
cor(t(result$SCDC$p),t(scData$Indata$P))
colSums(scData$Indata$P)

result$SCDC$result = getPearsonRMSE(result$SCDC$p,scData$Indata$P)
result$SCDC$result

#---------DWLS--------
# require(DWLS)
source('DWLS.R')
getwd()
path=paste(getwd(),"/DWLS/results_",STRING_name,sep=""); path
if(! dir.exists(path)){ #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created
  
  dir.create(path)
  Signature <- buildSignatureMatrixMAST(scdata = scData$sc$C.eset@assayData$exprs,
                                        id = as.character(scData$sc$phenoDataC$cellType), 
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

result$DWLS$p <- apply(scData$sc$T,2, function(x){
  b = setNames(x, rownames(scData$sc$T))
  tr <- trimData(Signature, b)
  RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
})

rownames(result$DWLS$p) <- as.character(unique(scData$sc$phenoDataC$cellType))
result$DWLS$p = apply(result$DWLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DWLS$p = apply(result$DWLS$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DWLS$result = getPearsonRMSE(result$DWLS$p,scData$Indata$P)
result$DWLS$result

# saveRDS(result,paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))

#-----complete deconvolution methods----------
#--------DSA-------------------
require(CellMix)
md = scData$simulate1$markerslist
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

result$DSA = CellMix::ged(as.matrix(scData$Indata$T), ML, method = "DSA",verbose=TRUE,  
                          sscale = TRUE, exact = FALSE, maxIter=500, log=FALSE)
result$DSA$c = result$DSA@fit@W
result$DSA$p = result$DSA@fit@H
result$DSA$t = result$DSA$c%*%result$DSA$p

result$DSA$result$c = getPearsonRMSE(result$DSA$c,scData$Indata$C)
result$DSA$result$p = getPearsonRMSE(result$DSA$p,scData$Indata$P)
result$DSA$result$t = getPearsonRMSE(result$DSA$t,scData$Indata$T)

result$DSA$result$all <- cbind(result$DSA$result$p,result$DSA$result$c)
result$DSA$result$all <-  cbind(result$DSA$result$all, result$DSA$result$t)
colnames(result$DSA$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$DSA$result <- result$DSA$result$all
result$DSA$result

#---------ssKL------------------
require(CellMix)

result$ssKL <- CellMix::ged(as.matrix(scData$Indata$T), ML, method = "ssKL",verbose=TRUE, 
                            sscale = FALSE, maxIter=500, log = FALSE)
result$ssKL$c = result$ssKL@fit@W
result$ssKL$p = result$ssKL@fit@H
result$ssKL$t = result$ssKL$c%*%result$ssKL$p

result$ssKL$result$c = getPearsonRMSE(result$ssKL$c,scData$Indata$C)
result$ssKL$result$p = getPearsonRMSE(result$ssKL$p,scData$Indata$P)
result$ssKL$result$t = getPearsonRMSE(result$ssKL$t,scData$Indata$T)

result$ssKL$result$all <- cbind(result$ssKL$result$p,result$ssKL$result$c)
result$ssKL$result$all <-  cbind(result$ssKL$result$all,result$ssKL$result$t)
colnames(result$ssKL$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssKL$result <- result$ssKL$result$all
result$ssKL$result

#----------ssFrobenius-----------------
require(CellMix)

result$ssFrobenius <- CellMix::ged(as.matrix(scData$Indata$T), ML, method = "ssFrobenius",verbose=TRUE,  
                                   sscale = FALSE, maxIter = 500, log = FALSE) #equivalent to coef(CellMix::ged(T,...)
result$ssFrobenius$c = result$ssFrobenius@fit@W
result$ssFrobenius$p = result$ssFrobenius@fit@H
result$ssFrobenius$t = result$ssFrobenius$c%*%result$ssFrobenius$p

result$ssFrobenius$result$c = getPearsonRMSE(result$ssFrobenius$c,scData$Indata$C)
result$ssFrobenius$result$p = getPearsonRMSE(result$ssFrobenius$p,scData$Indata$P)
result$ssFrobenius$result$t = getPearsonRMSE(result$ssFrobenius$t,scData$Indata$T)

result$ssFrobenius$result$all <- cbind(result$ssFrobenius$result$p,result$ssFrobenius$result$c)
result$ssFrobenius$result$all <-  cbind(result$ssFrobenius$result$all,result$ssFrobenius$result$t)
colnames(result$ssFrobenius$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssFrobenius$result <- result$ssFrobenius$result$all
result$ssFrobenius$result

#-----------deconf------------
# library(deconf)
all(rownames(scData$Indata$T)== md$gene)
all(rownames(scData$Indata$T)== rownames(md))

# result$deconf <- CellMix::ged(as.matrix(scData$Indata$T), ML, method = "deconf")
# result$deconf$c = result$deconf@fit@W
# result$deconf$p = result$deconf@fit@H
# result$deconf$t = result$deconf$c%*%result$deconf$p
set.seed(1234)
result$deconf <- CellMix::ged(as.matrix(scData$Indata$T), ncol(scData$Indata$C_ref))

result$deconf$c = result$deconf@fit@W
result$deconf$p = result$deconf@fit@H
result$deconf$t = result$deconf$c%*%result$deconf$p
ctlabels <- Row_label(result$deconf$c,scData$Indata$C_ref,leastnum=3);ctlabels
colnames(result$deconf$c) <- ctlabels
rownames(result$deconf$p) <- ctlabels

result$deconf$result$c = getPearsonRMSE(result$deconf$c,scData$Indata$C)
result$deconf$result$p = getPearsonRMSE(result$deconf$p,scData$Indata$P)
result$deconf$result$t = getPearsonRMSE(result$deconf$t,scData$Indata$T)

result$deconf$result$all <- cbind(result$deconf$result$p,result$deconf$result$c)
result$deconf$result$all <-  cbind(result$deconf$result$all,result$deconf$result$t)
colnames(result$deconf$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$deconf$result <- result$deconf$result$all
result$deconf$result

#------TOAST + NMF-------
require(DeCompress)
source("DeCompress.R")
result$TOAST$toast.nmf <- csDeCompress_my(Y_raw = scData$Indata$T,
                                          K = dim(scData$Indata$C_ref)[2],
                                          nMarker = nrow(scData$Indata$T),
                                          # FUN = nmfOut,
                                          TotalIter = 30)
require(NMF)
result$TOAST$fin.nmf = nmf(x = scData$Indata$T,
                           rank = dim(scData$Indata$C_ref)[2])
result$TOAST$ppp = t(coef(result$TOAST$fin.nmf))
result$TOAST$ppp = t(apply(result$TOAST$ppp,1,function(c) c/sum(c)))
nmf.res = list(prop = t(result$TOAST$ppp),sig = basis(result$TOAST$fin.nmf))

result$TOAST$p <- t(result$TOAST$ppp)
result$TOAST$c <- basis(result$TOAST$fin.nmf)
result$TOAST$t <- result$TOAST$c%*%result$TOAST$p
nmf.res = list(prop = t(result$TOAST$ppp),
               sig = basis(result$TOAST$fin.nmf))

rownames(scData$Indata$P)
labels <- Row_label(t(result$TOAST$p),t(scData$Indata$P));labels
rownames(result$TOAST$p) <- labels
colnames(result$TOAST$c) <- labels
result$TOAST$p <- result$TOAST$p[rownames(scData$Indata$P),]
result$TOAST$c <- result$TOAST$c[,colnames(scData$Indata$C)]

result$TOAST$result$c = getPearsonRMSE(result$TOAST$c,scData$Indata$C)
result$TOAST$result$p = getPearsonRMSE(result$TOAST$p,scData$Indata$P)
result$TOAST$result$t = getPearsonRMSE(result$TOAST$t,scData$Indata$T)

result$TOAST$result$all <- NULL
result$TOAST$result$all <- cbind(result$TOAST$result$p,result$TOAST$result$c)
result$TOAST$result$all <-  cbind(result$TOAST$result$all,result$TOAST$result$t)
colnames(result$TOAST$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$TOAST$result <- result$TOAST$result$all
result$TOAST$result

#-------Linseed---------------
require(DeCompress)
result$Linseed$Linseed.rs = DeCompress::linCor(yref = scData$Indata$T,
                                               iters = 100,
                                               pval = 100,
                                               n.types = dim(scData$Indata$C_ref)[2],
                                               scree = 'drop',
                                               logTransform = T)
table(is.na(result$Linseed$Linseed.rs$sig))
result$Linseed$c <- result$Linseed$Linseed.rs$sig
result$Linseed$p <- result$Linseed$Linseed.rs$prop
result$Linseed$t <- result$Linseed$c%*%result$Linseed$p
rownames(scData$Indata$P)
labels <- Row_label(t(result$Linseed$p),t(scData$Indata$P));labels
rownames(result$Linseed$p) <- labels
colnames(result$Linseed$c) <- labels
result$Linseed$p <- result$Linseed$p[rownames(scData$Indata$P),]
result$Linseed$c <- result$Linseed$c[,colnames(scData$Indata$C)]

result$Linseed$result$c = getPearsonRMSE(result$Linseed$c,scData$Indata$C)
result$Linseed$result$p = getPearsonRMSE(result$Linseed$p,scData$Indata$P)
result$Linseed$result$t = getPearsonRMSE(result$Linseed$t,scData$Indata$T)

result$Linseed$result$all <- cbind(result$Linseed$result$p,result$Linseed$result$c)
result$Linseed$result$all <-  cbind(result$Linseed$result$all,result$Linseed$result$t)
colnames(result$Linseed$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$Linseed$result <- result$Linseed$result$all
result$Linseed$result

#-----CellDistinguisher-----
require(DeCompress)
result$CellDistinguisher <- CellDistinguisher::gecd_CellDistinguisher(
  scData$Indata$T,
  genesymb = rownames(scData$Indata$T),
  numCellClasses = dim(scData$Indata$C_ref)[2],
  minDistinguisherAlternatives=1,
  maxDistinguisherAlternatives=100,
  minAlternativesLengthsNormalized=0.5,
  expressionQuantileForScale = 0.75,
  expressionQuantileForFilter=0.999,
  expressionConcentrationRatio=0.333,
  probesWithGenesOnly = F,
  verbose=0)

# source("CellDistinguisher.R")
result$CellDist.deconv <-
  tryCatch(CellDistinguisher::gecd_DeconvolutionByDistinguishers(
    as.matrix(scData$Indata$T),
    result$CellDistinguisher$bestDistinguishers,
    nonNegativeOnly = T,
    convexSolution = T,
    verbose = 0),
    error = function(e) return(list(sampleCompositions =
                                      matrix(rep(1/dim(scData$Indata$C_ref)[2],
                                                 dim(scData$Indata$C_ref)[2]*ncol(scData$Indata$T)),
                                             ncol=dim(scData$Indata$C_ref)[2]))))
# library(CellMix)
# deconvolutionSSKL <- CellDistinguisher::gecd_DeconvolutionCellMix(
#   as.matrix(scData$Indata$T), 
#   result$CellDistinguisher$bestDistinguishers, 
#   method="ssKL", 
#   maxIter=dim(scData$Indata$C_ref)[2])

if (length(result$CellDist.deconv) > 1){
  result$CellDistinguisher$p <- result$CellDist.deconv$sampleComposition
  result$CellDistinguisher$c <- result$CellDist.deconv$cellSubclassSignatures
}else {
  result$CellDistinguisher$p <- nmf.res$prop
  result$CellDistinguisher$c <- nmf.res$sig
}
# result$CellDistinguisher$c[1:20,]
# result$CellDistinguisher$p[,1:20]
result$CellDistinguisher$t <- result$CellDistinguisher$c%*%result$CellDistinguisher$p
rownames(scData$Indata$P)
labels <- Row_label(t(result$CellDistinguisher$p),t(scData$Indata$P))
rownames(result$CellDistinguisher$p) <- labels
colnames(result$CellDistinguisher$c) <- labels
result$CellDistinguisher$p <- result$CellDistinguisher$p[rownames(scData$Indata$P),]
result$CellDistinguisher$c <- result$CellDistinguisher$c[,colnames(scData$Indata$C)]
result$CellDistinguisher$result$c = getPearsonRMSE(result$CellDistinguisher$c,scData$Indata$C)
result$CellDistinguisher$result$p = getPearsonRMSE(result$CellDistinguisher$p,scData$Indata$P)
result$CellDistinguisher$result$t = getPearsonRMSE(result$CellDistinguisher$t,scData$Indata$T)
result$CellDistinguisher$result$all <- cbind(result$CellDistinguisher$result$p,result$CellDistinguisher$result$c)
result$CellDistinguisher$result$all <-  cbind(result$CellDistinguisher$result$all,result$CellDistinguisher$result$t)
colnames(result$CellDistinguisher$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$CellDistinguisher$result <- result$CellDistinguisher$result$all
result$CellDistinguisher$result


#--------------all-----------------
MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,
                  "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",
                  "MuSiC","Bisque","SCDC", "DWLS",
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
result <- CountAllResults(result,MyMethodName)
result$all;dim(result$all)
paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep="")
saveRDS(result,paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))


#----------find Paramater ///----------
lambda1 <- c(0,10^-3,10^-2,10^-1,1,10)
lambda2 <- c(0,10^-3,10^-2,10^-1,1,10)
lambdaC <- c(0,10^4,1,10^1,10^2,10^3)
result$CDSC3$seedd = 44
scData$Indata$TerCondition = 10^-8

#------------cyclic to do deconvolution----------
pb <- txtProgressBar(style = 3)
star_time <- Sys.time()
result$CDSC3$Ss <- SM(t(scData$Indata$T))
result$CDSC3$Sg <- SM(scData$Indata$T)
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
      result_CDSC = CDSC_3(scData$Indata$T, scData$Indata$C_ref, dim(scData$Indata$C_ref)[2], 
                           lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k],
                           scData$Indata$TerCondition,result$CDSC3$seedd,
                           result$CDSC3$Ss,result$CDSC3$Sg,all_number = 3000)
      result_para_c[[num]] <- result_CDSC[[1]]
      result_para_p[[num]] <- result_CDSC[[2]]
      number_iter[num] <- result_CDSC[[3]]
      
      pearson_para_c[[num]] <- cor(result_para_c[[num]],scData$Indata$C_ref)
      result1 = NULL
      if(!all(is.na(pearson_para_c[[num]]) == FALSE) ){
        break
      }
      result1 <- calculate_result(result_para_c[[num]],result_para_p[[num]],
                                  scData$Indata$T,scData$Indata$C,scData$Indata$C_ref,scData$Indata$P,
                                  lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k],number_iter[num],
                                  result$CDSC3$seedd,scData$Indata$TerCondition)
      # result
      # rownames(result) = 1:length(rownames(result))
      para_lambda_44_8 =  rbind(para_lambda_44_8,result1)
      # cat("\nlambda1=",lambda1[dir_i],",lambda2=",lambda2[dir_j],
      #     ",lambdac=",lambdaC[dir_k],",杩唬娆℃暟",number_iter[num])
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



# para_lambda_44_8_12 <- para_lambda_44_8
saveRDS(para_lambda_44_8,
        paste(getwd(),"/paramater/para_lambda_",STRING_name,"_NEW.rds",sep=""))

# pheatMap
rm (list=ls ())

a = c("Nestorowa","Manno","Darmanis","Camp","Segerstolpe")

MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,     
                  "deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso" ,"EPIC",          
                  "MuSiC","Bisque","SCDC", "DWLS",                      
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
scData = NULL
result = NULL 
for(i in 1:length(a)){
  STRING_name = a[i]; STRING_name
  # scData[[i]] <- readRDS(paste(getwd(),"/dataSimulate/",STRING_name,".rds",sep=""))
  result[[i]] <- readRDS(paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))
}

names(result[[1]])
map = NULL
for(i in 1:length(a)){
  # result[[i]] <- CountAllResults(result[[i]],MyMethodName)
  result[[i]]$all;dim(result[[i]]$all)
  if(nrow(result[[i]]$all)!=length(MyMethodName)){
    print("ERROR")
    break;}
  map <- cbind(map,result[[i]]$all$RMSE_to_C)
}
rownames(map) <- MyMethodName
colnames(map) <- a
map;dim(map)
# "Pearson"  "RMSE"
library(pheatmap)#c("navy", "white", "firebrick3") color
plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = F
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
                         ,cellwidth = 30, cellheight = 12#,gaps_row = c(12, 17)
                         ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.3f"
                         ,border_color = "black",color = colorRampPalette(c("white", "firebrick3"))(30)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "Pearson of C")
plot_pheatmp <- pheatmap(map[c(1,17:24),],cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 30, cellheight = 12#,gaps_row = c(12, 17)
                         # ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.0f"
                         ,border_color = "black",color = colorRampPalette(c("white", "navy"))(30)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,main = "RMSE of C")

plot_pheatmp;

eoffice::topptx(plot_pheatmp,
                filename = "pictures/class1_NO_RMSE_C.pptx")

# boxplot--p----
mapBox <- NULL
mapBox =  data.frame(pearson = map[1,], method = rownames(map)[1],row.names = NULL)
nn1 <- nrow(mapBox)
for (i in 1:length(rownames(map))) {
  mapBox <- rbind(mapBox, data.frame(pearson = map[i,], method = rownames(map)[i],row.names = NULL))
}

mapBox <- mapBox[-(1:nn1),]
library(ggplot2)
plot_class1 =  ggplot(mapBox, aes(factor(method,levels=MyMethodName),pearson,fill=method)) +  # background
  stat_boxplot(geom = "errorbar",width=0.3)+
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.05),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson of P")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 搴曡壊
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #璁剧疆鏍囬灞呬腑
plot_class1;
eoffice::topptx(plot_class1,
                filename = "pictures/class1_NO_fanhua_p.pptx")

# boxplot----c-----
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
  ylim(0, 1)+
  labs(x="Methods",y = "Pearson of C")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 搴曡壊
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #璁剧疆鏍囬灞呬腑
plot_class1;
eoffice::topptx(plot_class1,
                filename = "pictures/class1_NO_fanhua_C.pptx")

# sample bar
rm (list=ls ())
setwd("1scSimulate")

source("CDSC.R")
source("CDSC_expand.R")

# myOneref <- intersect(Oneref,MethodName);myOneref=setdiff(myOneref,c("CDSC2"))
# myTworef <- intersect(NoCref,MethodName);myTworef=union(c("CDSC3"),myTworef);myTworef=setdiff(myTworef,c("CDSC2"))
a = gsub('.rds','',list.files("1scSimulate/dataSimulate"));a
a = c("Nestorowa","Manno","Darmanis","Camp","Segerstolpe")
scData = NULL
result = NULL 
for(i in 1:length(a)){
  STRING_name = a[i]; STRING_name
  scData[[i]] <- readRDS(paste(getwd(),"/dataSimulate/",STRING_name,".rds",sep=""))
  result[[i]] <- readRDS(paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))
}
map = NULL
for(i in 1:length(a)){
  map <- rbind(map,data.frame(group = a[i],
                              pearson = GetSampleCor(result[[i]]$CDSC3$dec$p,scData[[i]]$Indata$P,matrix = "p")))
}
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
  geom_bar(data=map_mean,mapping=aes(x=factor(group,levels = a),y=mean,fill=group), # fill濉厖
           position="dodge", # 鏌辩姸鍥炬牸寮�
           stat="identity", # 鏁版嵁鏍煎紡
           width = 0.7,
           show.legend = F)+  # 鏌辩姸鍥惧昂瀵�
  scale_fill_manual(values = c("#80AFBF","#80AFBF","#80AFBF","#80AFBF","#80AFBF"))+ # 鏌辩姸鍥鹃鑹�, "#DA635D","#B1938B"
  # geom_signif(data=plot_data2,mapping=aes(x=group,y=SOD), # 涓嶅悓缁勫埆鐨勬樉钁楁€�
  #             comparisons = list(c("C", "HT"), # 鍝簺缁勮繘琛屾瘮杈�
  #                                c("HI", "HT")),
  #             annotation=c("**"), # 鏄捐憲鎬у樊寮傚仛鏍囪
  #             map_signif_level=T, # T涓烘樉钁楁€э紝F涓簆 value
  #             tip_length=c(0.04,0.04,0.05,0.05), # 淇敼鏄捐憲鎬ч偅涓嚎鐨勯暱鐭�
  #             y_position = c(4100,3000), # 璁剧疆鏄捐憲鎬х嚎鐨勪綅缃珮搴�
  #             size=1, # 淇敼绾跨殑绮楃粏
  #             textsize = 10, # 淇敼*鏍囪鐨勫ぇ灏�
  #             test = "t.test")+ # 妫€楠岀殑绫诲瀷
  geom_errorbar(data=map_mean,mapping=aes(x = group,ymin = mean-sd, ymax = mean+sd), # 璇樊绾挎坊鍔�
                width = 0.1, #璇樊绾跨殑瀹藉害
                color = 'black', #棰滆壊
                size=0.8)+ #绮楃粏
  scale_y_continuous(limits =c(0, 1.1) ,expand = c(0,0))+ # y杞寸殑鑼冨洿
  theme_classic(  # 涓婚璁剧疆锛岃繖涓槸鏃犵嚎鏉′富棰�
    base_line_size = 1 # 鍧愭爣杞寸殑绮楃粏
  )+
  
  labs(title="",x="Dataset",y="Pearson of samples in P")+ # 娣诲姞鏍囬锛寈杞达紝y杞村唴瀹�
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
        legend.title = element_text(color="black", # 淇敼鍥句緥鐨勬爣棰�
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", # 璁剧疆鍥句緥鏍囩鏂囧瓧
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # 淇敼X杞翠笂瀛椾綋澶у皬锛�
                                   # family = "myFont", # 绫诲瀷
                                   color = "black", # 棰滆壊
                                   face = "bold", #  face鍙栧€硷細plain鏅€氾紝bold鍔犵矖锛宨talic鏂滀綋锛宐old.italic鏂滀綋鍔犵矖
                                   vjust = 1, # 浣嶇疆
                                   hjust = 1, 
                                   angle = 45), #瑙掑害
        axis.text.y = element_text(size = 13,  # 淇敼y杞翠笂瀛椾綋澶у皬锛�
                                   # family = "myFont", # 绫诲瀷
                                   color = "black", # 棰滆壊
                                   face = "bold", #  face鍙栧€硷細plain鏅€氾紝bold鍔犵矖锛宨talic鏂滀綋锛宐old.italic鏂滀綋鍔犵矖
                                   vjust = 0.5, # 浣嶇疆
                                   hjust = 0.5, 
                                   angle = 0) #瑙掑害
  ) 
# emf(file = "SOD.emf") # 鎵撳紑涓€涓煝閲忓浘鐢诲竷锛岃繖绉嶆牸寮忕殑鍥剧墖鏀惧湪word閲屼笉浼氬け鐪�
print(p1) # 鎵撳嵃鍥剧墖
eoffice::topptx(p1,
                filename = "pictures/class1_bar_NO_samples.pptx")

# ---- sample鎯呭喌鐨勬暎鐐瑰浘----
library(ggplot2)
plot_class2 =  ggplot(map, aes(factor(group,levels=a),pearson,fill=group)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 搴曡壊
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #璁剧疆鏍囬灞呬腑
plot_class2
eoffice::topptx(plot_class2,filename = "pictures/class1_boxplot_samples.pptx")

