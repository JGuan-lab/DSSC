source("CDSC.R")
source("CDSC_help.R")

bulkData <- readRDS("XXX.rds")

result <- list()
result$CDSC3$seedd = 44
result$CDSC3$TerCondition = 10^-8
result$CDSC3$Ss <- SM(t(bulkData$Indata$T))
result$CDSC3$Sg <- SM(bulkData$Indata$T)
# 
result$CDSC3$lambda1 <- 1e-03
result$CDSC3$lambda2 <- 0e+00
result$CDSC3$lambdaC <- 100

library(dplyr)

# result
result$CDSC3$dec = CDSC_3(bulkData$Indata$T, bulkData$Indata$C_ref, dim(bulkData$Indata$C_ref)[2], 
                          result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                          result$CDSC3$TerCondition,result$CDSC3$seedd,
                          result$CDSC3$Ss,result$CDSC3$Sg,all_number = 3000)

result$CDSC3$result <- calculate_result(result$CDSC3$dec$c,result$CDSC3$dec$p,
                                        bulkData$Indata$T,bulkData$Indata$C,bulkData$Indata$C_ref,bulkData$Indata$P,
                                        result$CDSC3$lambda1, result$CDSC3$lambda2, result$CDSC3$lambdaC,
                                        number_iter = result$CDSC3$dec$jump,
                                        seedd = result$CDSC3$seedd,
                                        TerCondition = result$CDSC3$TerCondition)
result$CDSC3$result

result$CDSC2$lambda1 <- 1e-02
result$CDSC2$lambda2 <- 1e+01
result$CDSC2$lambdaC <- 0
result$CDSC2$dec = CDSC_2(bulkData$Indata$T,  dim(bulkData$Indata$C_ref)[2], 
                          result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                          result$CDSC3$TerCondition,result$CDSC3$seedd,
                          result$CDSC3$Ss,result$CDSC3$Sg,all_number = 3000)

result$CDSC2$result <- calculate_result(result$CDSC2$dec$c,result$CDSC2$dec$p,
                                        bulkData$Indata$T, bulkData$Indata$C, bulkData$Indata$C_ref, bulkData$Indata$P,
                                        result$CDSC2$lambda1, result$CDSC2$lambda2, result$CDSC2$lambdaC,
                                        number_iter = result$CDSC2$dec$jump,
                                        seedd = result$CDSC3$seedd,
                                        TerCondition = result$CDSC3$TerCondition)
result$CDSC2$result

#----------NNLS----------
require(nnls)
result$NNLS$p = do.call(cbind.data.frame,lapply(apply(bulkData$Indata$T,2,function(x) nnls::nnls(as.matrix(bulkData$Indata$C_ref),x)),   function(y) y$x))
result$NNLS$p = apply(result$NNLS$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$NNLS$p) <- colnames(bulkData$Indata$C_ref)
result$NNLS$result = getPearsonRMSE(result$NNLS$p,bulkData$Indata$P)
result$NNLS$result
diag(cor(t(result$NNLS$p),t(bulkData$Indata$P)))

#--------OLS------------
result$OLS$p = apply(bulkData$Indata$T,2,function(x) lm(x ~ as.matrix(bulkData$Indata$C_ref))$coefficients[-1])
result$OLS$p = apply(result$OLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$OLS$p = apply(result$OLS$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$OLS$p) <- unlist(lapply(strsplit(rownames(result$OLS$p),")"),function(x) x[2]))
result$OLS$result = getPearsonRMSE(result$OLS$p,bulkData$Indata$P)
result$OLS$result


#---------FARDEEP-----------
library(FARDEEP)
#result_FARDEEP = t(FARDEEP(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = t(FARDEEP::fardeep(bulkData$Indata$C_ref, bulkData$Indata$T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
result$FARDEEP$p = apply(result$FARDEEP$p,2,function(x) x/sum(x)) #explicit STO constraint
result$FARDEEP$result = getPearsonRMSE(result$FARDEEP$p,bulkData$Indata$P)
result$FARDEEP$result


#----------CIBERSORT-----------------
source("CIBERSORT.R")
result$CIBERSORT$p = CIBERSORT(sig_matrix = bulkData$Indata$C_ref, mixture_file = bulkData$Indata$T, QN = FALSE)
result$CIBERSORT$p = t(result$CIBERSORT$p[,1:(ncol(result$CIBERSORT$p)-3)])
result$CIBERSORT$result = getPearsonRMSE(result$CIBERSORT$p,bulkData$Indata$P)
result$CIBERSORT$result
diag(cor(t(result$CIBERSORT$p),t(bulkData$Indata$P)))

#----------------DeconRNASeq------
#nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
#datasets and reference matrix: signatures, need to be non-negative. 
#"use.scale": whether the data should be centered or scaled, default = TRUE
unloadNamespace("Seurat") #needed for PCA step
library(pcaMethods) #needed for DeconRNASeq to work
result$deconRNASeq$p = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(bulkData$Indata$T), 
                                                  signatures = as.data.frame(bulkData$Indata$C_ref), 
                                                  proportions = NULL, checksig = FALSE, 
                                                  known.prop = FALSE, use.scale = TRUE,
                                                  fig = FALSE)$out.all)
colnames(result$deconRNASeq$p) = colnames(bulkData$Indata$T)
result$deconRNASeq$p[,1:10]
require(Seurat)
result$deconRNASeq$result = getPearsonRMSE(result$deconRNASeq$p,bulkData$Indata$P)
result$deconRNASeq$result


#-----------RLR----------------
require(MASS)
result$RLR$p = do.call(cbind.data.frame,lapply(apply(bulkData$Indata$T,2,function(x) MASS::rlm(x ~ as.matrix(bulkData$Indata$C_ref), maxit=100)), function(y) y$coefficients[-1]))
result$RLR$p = apply(result$RLR$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$RLR$p = apply(result$RLR$p,2,function(x) x/sum(x)) #explicit STO constraint
rownames(result$RLR$p) <- unlist(lapply(strsplit(rownames(result$RLR$p),")"),function(x) x[2]))
result$RLR$result = getPearsonRMSE(result$RLR$p,bulkData$Indata$P)
result$RLR$result


#-----------DCQ------------------
#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default
require(ComICS)
result$DCQ$p = t(ComICS::dcq(reference_data = bulkData$Indata$C_ref, 
                             mix_data = bulkData$Indata$T, 
                             marker_set = as.data.frame(row.names(bulkData$Indata$C_ref)) , 
                             alpha_used = 0.99, 
                             lambda_min = 0.1, 
                             number_of_repeats = 10)$average)
result$DCQ$p = apply(result$DCQ$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DCQ$p = apply(result$DCQ$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DCQ$result = getPearsonRMSE(result$DCQ$p,bulkData$Indata$P)
result$DCQ$result
# cor(t(result$DCQ$p),t(bulkData$Indata$P))

#-----------elastic_net-----------
#standardize = TRUE by default. lambda=NULL by default 
require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
result$elastic_net$p = apply(bulkData$Indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(bulkData$Indata$C_ref), y = z, alpha = 0.2, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(bulkData$Indata$C_ref), z)$lambda.1se))[1:ncol(bulkData$Indata$C_ref)+1,])
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$elastic_net$p = apply(result$elastic_net$p,2,function(x) x/sum(x)) #explicit STO constraint
result$elastic_net$result = getPearsonRMSE(result$elastic_net$p,bulkData$Indata$P)
result$elastic_net$result


#----------ridge----------------
# alpha=0
require(glmnet)
result$ridge$p = apply(bulkData$Indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(bulkData$Indata$C_ref), y = z, alpha = 0, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(bulkData$Indata$C_ref), z)$lambda.1se))[1:ncol(bulkData$Indata$C_ref)+1,])
result$ridge$p = apply(result$ridge$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$ridge$p = apply(result$ridge$p,2,function(x) x/sum(x)) #explicit STO constraint
result$ridge$result = getPearsonRMSE(result$ridge$p,bulkData$Indata$P)
result$ridge$result


#----------lasso-----------------
#alpha=1; shrinking some coefficients to 0. 
require(glmnet)
result$lasso$p = apply(bulkData$Indata$T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(bulkData$Indata$C_ref), y = z, alpha = 1, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(bulkData$Indata$C_ref), z)$lambda.1se))[1:ncol(bulkData$Indata$C_ref)+1,])
result$lasso$p = apply(result$lasso$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$lasso$p = apply(result$lasso$p,2,function(x) x/sum(x)) #explicit STO constraint
result$lasso$p[which(is.na(result$lasso$p) == TRUE)] <- 0 #Needed for models where glmnet drops all terms of a model and fit an intercept-only model (very unlikely but possible).
length(colSums(result$lasso$p) == 1)
result$lasso$result = getPearsonRMSE(result$lasso$p,bulkData$Indata$P)
result$lasso$result

# source("CDSC.R")
# source("CDSC_help.R")

#----------EPIC----------
require(EPIC)
C_EPIC <- list()
C_EPIC$marker = bulkData$simulate1$markerslist[bulkData$simulate1$markerslist$gene %in% rownames(bulkData$Indata$C_ref),]
C_EPIC$markers = as.character(C_EPIC$marker$gene)
C_EPIC$CTs <- intersect(colnames(bulkData$Indata$C_ref),colnames(bulkData$Indata$refProfiles.var))
C_EPIC[["sigGenes"]] <- rownames(bulkData$Indata$C_ref[C_EPIC$markers,C_EPIC$CTs])
C_EPIC[["refProfiles"]] <- as.matrix(bulkData$Indata$C_ref[C_EPIC$markers,C_EPIC$CTs])
C_EPIC[["refProfiles.var"]] <- bulkData$Indata$refProfiles.var[C_EPIC$markers,C_EPIC$CTs]

result$EPIC$p <- t(EPIC::EPIC(bulk=as.matrix(bulkData$Indata$T), 
                              reference=C_EPIC, withOtherCells=TRUE, 
                              scaleExprs=TRUE)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices
result$EPIC$p = result$EPIC$p[!rownames(result$EPIC$p) %in% "otherCells",]
result$EPIC$result = getPearsonRMSE(result$EPIC$p,bulkData$Indata$P)
result$EPIC$result

#--------sc methods_------------
bulkData$sc$keep <- which(rownames(bulkData$simulate1$split$train) %in% bulkData$simulate1$shif_marker$gene)
# bulkData$sc$keep <- which(rownames(bulkData$simulate1$split$train) %in% rownames(bulkData$simulate1$split$train))
bulkData$sc$T <- bulkData$simulate1$T[bulkData$sc$keep,]

bulkData$sc$C <- bulkData$simulate1$split$train[bulkData$sc$keep,]
colnames(bulkData$sc$C) <- bulkData$simulate1$original_cell_names[bulkData$simulate1$split$training]
bulkData$sc$C <- bulkData$sc$C[rownames(bulkData$sc$T),]
bulkData$sc$phenoDataC <- bulkData$simulate1$split$pData_train

#Bisque requires "SubjectName" in phenoDataC
if(length(grep("[N-n]ame",colnames(bulkData$sc$phenoDataC))) > 0){
  bulkData$sc$sample_column = grep("[N-n]ame",colnames(bulkData$sc$phenoDataC))
} else {
  bulkData$sc$sample_column = grep("[S-s]ample|[S-s]ubject",colnames(bulkData$sc$phenoDataC))
}

colnames(bulkData$sc$phenoDataC)[bulkData$sc$sample_column] = "SubjectName"
rownames(bulkData$sc$phenoDataC) = bulkData$sc$phenoDataC$cellID

require(xbioc)
colnames(bulkData$sc$C) <- bulkData$simulate1$original_cell_names[bulkData$simulate1$split$training]
bulkData$sc$C.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulkData$sc$C)
                                           ,phenoData = Biobase::AnnotatedDataFrame(bulkData$sc$phenoDataC))
bulkData$sc$T.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulkData$sc$T))

#---------MuSiC-----------
library(MuSiC)
result$MuSiC$p = t(MuSiC::music_prop(bulk.eset = bulkData$sc$T.eset, 
                                 sc.eset = bulkData$sc$C.eset,
                                 markers = NULL,
                                 clusters = "cellType",
                                 samples = 'SubjectName',
                                 select.ct = unique(as.character(bulkData$sc$phenoDataC$cellType)),
                                 verbose = F,
                                 normalize = F
                                 
)$Est.prop.weighted)
result$MuSiC$p[,1:10]
bulkData$Indata$P[,1:10]
result$MuSiC$result = getPearsonRMSE(result$MuSiC$p,bulkData$Indata$P)
result$MuSiC$result

#------------Bisque----------
require(Bisque)
result$Bisque$p <- BisqueRNA::ReferenceBasedDecomposition(bulkData$sc$T.eset, bulkData$sc$C.eset, 
                                                          markers=NULL, use.overlap=FALSE)$bulk.props 
#use.overlap is when there's both bulk and scRNA-seq for the same set of samples
result$Bisque$result = getPearsonRMSE(result$Bisque$p,bulkData$Indata$P)
result$Bisque$result

#--------SCDC-----
library(SCDC)
source("SCDC.R")
result$SCDC$p <- t(SCDC::SCDC_prop(bulk.eset = bulkData$sc$T.eset,
                                sc.eset = bulkData$sc$C.eset, 
                                ct.varname = "cellType", 
                                sample = "SubjectName", 
                                ct.sub = unique(as.character(bulkData$sc$phenoDataC$cellType)), 
                                iter.max = 1000)$prop.est.mvw)
result$SCDC$p[,1:5]
bulkData$Indata$P[,1:5]
result$SCDC$p <- result$SCDC$p[rownames(bulkData$Indata$P),]
cor(t(result$SCDC$p),t(bulkData$Indata$P))
colSums(bulkData$Indata$P)

result$SCDC$result = getPearsonRMSE(result$SCDC$p,bulkData$Indata$P)
result$SCDC$result

#---------DWLS--------
# require(DWLS)
source('DWLS.R')
getwd()
path=paste(getwd(),"/DWLS/results_",STRING_name,sep=""); path
if(! dir.exists(path)){ #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created
  
  dir.create(path)
  Signature <- buildSignatureMatrixMAST(bulkData = bulkData$sc$C.eset@assayData$exprs,
                                        id = as.character(bulkData$sc$phenoDataC$cellType), 
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

result$DWLS$p <- apply(bulkData$sc$T,2, function(x){
  b = setNames(x, rownames(bulkData$sc$T))
  tr <- trimData(Signature, b)
  RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
})

rownames(result$DWLS$p) <- as.character(unique(bulkData$sc$phenoDataC$cellType))
result$DWLS$p = apply(result$DWLS$p,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
result$DWLS$p = apply(result$DWLS$p,2,function(x) x/sum(x)) #explicit STO constraint
result$DWLS$result = getPearsonRMSE(result$DWLS$p,bulkData$Indata$P)
result$DWLS$result

# saveRDS(result,paste(getwd(),"/dataResult/result_",STRING_name,".rds",sep=""))

#-----complete deconvolution methods----------
#--------DSA-------------------
require(CellMix)
md = bulkData$simulate1$markerslist
ML = CellMix::MarkerList()
ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

result$DSA = CellMix::ged(as.matrix(bulkData$Indata$T), ML, method = "DSA",verbose=TRUE,  
                          sscale = TRUE, exact = FALSE, maxIter=500, log=FALSE)
result$DSA$c = result$DSA@fit@W
result$DSA$p = result$DSA@fit@H
result$DSA$t = result$DSA$c%*%result$DSA$p

result$DSA$result$c = getPearsonRMSE(result$DSA$c,bulkData$Indata$C)
result$DSA$result$p = getPearsonRMSE(result$DSA$p,bulkData$Indata$P)
result$DSA$result$t = getPearsonRMSE(result$DSA$t,bulkData$Indata$T)

result$DSA$result$all <- cbind(result$DSA$result$p,result$DSA$result$c)
result$DSA$result$all <-  cbind(result$DSA$result$all, result$DSA$result$t)
colnames(result$DSA$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$DSA$result <- result$DSA$result$all
result$DSA$result

#---------ssKL------------------
require(CellMix)

result$ssKL <- CellMix::ged(as.matrix(bulkData$Indata$T), ML, method = "ssKL",verbose=TRUE, 
                            sscale = FALSE, maxIter=500, log = FALSE)
result$ssKL$c = result$ssKL@fit@W
result$ssKL$p = result$ssKL@fit@H
result$ssKL$t = result$ssKL$c%*%result$ssKL$p

result$ssKL$result$c = getPearsonRMSE(result$ssKL$c,bulkData$Indata$C)
result$ssKL$result$p = getPearsonRMSE(result$ssKL$p,bulkData$Indata$P)
result$ssKL$result$t = getPearsonRMSE(result$ssKL$t,bulkData$Indata$T)

result$ssKL$result$all <- cbind(result$ssKL$result$p,result$ssKL$result$c)
result$ssKL$result$all <-  cbind(result$ssKL$result$all,result$ssKL$result$t)
colnames(result$ssKL$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssKL$result <- result$ssKL$result$all
result$ssKL$result

#----------ssFrobenius-----------------
require(CellMix)

result$ssFrobenius <- CellMix::ged(as.matrix(bulkData$Indata$T), ML, method = "ssFrobenius",verbose=TRUE,  
                                   sscale = FALSE, maxIter = 500, log = FALSE) #equivalent to coef(CellMix::ged(T,...)
result$ssFrobenius$c = result$ssFrobenius@fit@W
result$ssFrobenius$p = result$ssFrobenius@fit@H
result$ssFrobenius$t = result$ssFrobenius$c%*%result$ssFrobenius$p

result$ssFrobenius$result$c = getPearsonRMSE(result$ssFrobenius$c,bulkData$Indata$C)
result$ssFrobenius$result$p = getPearsonRMSE(result$ssFrobenius$p,bulkData$Indata$P)
result$ssFrobenius$result$t = getPearsonRMSE(result$ssFrobenius$t,bulkData$Indata$T)

result$ssFrobenius$result$all <- cbind(result$ssFrobenius$result$p,result$ssFrobenius$result$c)
result$ssFrobenius$result$all <-  cbind(result$ssFrobenius$result$all,result$ssFrobenius$result$t)
colnames(result$ssFrobenius$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$ssFrobenius$result <- result$ssFrobenius$result$all
result$ssFrobenius$result

#-----------deconf------------
# library(deconf)
all(rownames(bulkData$Indata$T)== md$gene)
all(rownames(bulkData$Indata$T)== rownames(md))

# result$deconf <- CellMix::ged(as.matrix(bulkData$Indata$T), ML, method = "deconf")
# result$deconf$c = result$deconf@fit@W
# result$deconf$p = result$deconf@fit@H
# result$deconf$t = result$deconf$c%*%result$deconf$p
set.seed(1234)
result$deconf <- CellMix::ged(as.matrix(bulkData$Indata$T), ncol(bulkData$Indata$C_ref))

result$deconf$c = result$deconf@fit@W
result$deconf$p = result$deconf@fit@H
result$deconf$t = result$deconf$c%*%result$deconf$p
ctlabels <- Row_label(result$deconf$c,bulkData$Indata$C_ref,leastnum=3);ctlabels
colnames(result$deconf$c) <- ctlabels
rownames(result$deconf$p) <- ctlabels

result$deconf$result$c = getPearsonRMSE(result$deconf$c,bulkData$Indata$C)
result$deconf$result$p = getPearsonRMSE(result$deconf$p,bulkData$Indata$P)
result$deconf$result$t = getPearsonRMSE(result$deconf$t,bulkData$Indata$T)

result$deconf$result$all <- cbind(result$deconf$result$p,result$deconf$result$c)
result$deconf$result$all <-  cbind(result$deconf$result$all,result$deconf$result$t)
colnames(result$deconf$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$deconf$result <- result$deconf$result$all
result$deconf$result

#------TOAST + NMF-------
require(DeCompress)
source("DeCompress.R")
result$TOAST$toast.nmf <- DeCompress::csDeCompress(Y_raw = bulkData$Indata$T,
                                          K = dim(bulkData$Indata$C_ref)[2],
                                          nMarker = nrow(bulkData$Indata$T),
                                          # FUN = nmfOut,
                                          TotalIter = 30)
require(NMF)
result$TOAST$fin.nmf = nmf(x = bulkData$Indata$T,
                           rank = dim(bulkData$Indata$C_ref)[2])
result$TOAST$ppp = t(coef(result$TOAST$fin.nmf))
result$TOAST$ppp = t(apply(result$TOAST$ppp,1,function(c) c/sum(c)))
nmf.res = list(prop = t(result$TOAST$ppp),sig = basis(result$TOAST$fin.nmf))

result$TOAST$p <- t(result$TOAST$ppp)
result$TOAST$c <- basis(result$TOAST$fin.nmf)
result$TOAST$t <- result$TOAST$c%*%result$TOAST$p
nmf.res = list(prop = t(result$TOAST$ppp),
               sig = basis(result$TOAST$fin.nmf))

rownames(bulkData$Indata$P)
labels <- Row_label(t(result$TOAST$p),t(bulkData$Indata$P));labels
rownames(result$TOAST$p) <- labels
colnames(result$TOAST$c) <- labels
result$TOAST$p <- result$TOAST$p[rownames(bulkData$Indata$P),]
result$TOAST$c <- result$TOAST$c[,colnames(bulkData$Indata$C)]

result$TOAST$result$c = getPearsonRMSE(result$TOAST$c,bulkData$Indata$C)
result$TOAST$result$p = getPearsonRMSE(result$TOAST$p,bulkData$Indata$P)
result$TOAST$result$t = getPearsonRMSE(result$TOAST$t,bulkData$Indata$T)

result$TOAST$result$all <- NULL
result$TOAST$result$all <- cbind(result$TOAST$result$p,result$TOAST$result$c)
result$TOAST$result$all <-  cbind(result$TOAST$result$all,result$TOAST$result$t)
colnames(result$TOAST$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$TOAST$result <- result$TOAST$result$all
result$TOAST$result

#-------Linseed---------------
require(DeCompress)
result$Linseed$Linseed.rs = DeCompress::linCor(yref = bulkData$Indata$T,
                                               iters = 100,
                                               pval = 100,
                                               n.types = dim(bulkData$Indata$C_ref)[2],
                                               scree = 'drop',
                                               logTransform = T)
table(is.na(result$Linseed$Linseed.rs$sig))
result$Linseed$c <- result$Linseed$Linseed.rs$sig
result$Linseed$p <- result$Linseed$Linseed.rs$prop
result$Linseed$t <- result$Linseed$c%*%result$Linseed$p
rownames(bulkData$Indata$P)
labels <- Row_label(t(result$Linseed$p),t(bulkData$Indata$P));labels
rownames(result$Linseed$p) <- labels
colnames(result$Linseed$c) <- labels
result$Linseed$p <- result$Linseed$p[rownames(bulkData$Indata$P),]
result$Linseed$c <- result$Linseed$c[,colnames(bulkData$Indata$C)]

result$Linseed$result$c = getPearsonRMSE(result$Linseed$c,bulkData$Indata$C)
result$Linseed$result$p = getPearsonRMSE(result$Linseed$p,bulkData$Indata$P)
result$Linseed$result$t = getPearsonRMSE(result$Linseed$t,bulkData$Indata$T)

result$Linseed$result$all <- cbind(result$Linseed$result$p,result$Linseed$result$c)
result$Linseed$result$all <-  cbind(result$Linseed$result$all,result$Linseed$result$t)
colnames(result$Linseed$result$all) <- c("RMSE_to_P", "Peason_to_P", "RMSE_to_C", "Peason_to_C", "RMSE_to_T", "Peason_to_T")
result$Linseed$result <- result$Linseed$result$all
result$Linseed$result

#-----CellDistinguisher-----
require(DeCompress)
result$CellDistinguisher <- CellDistinguisher::gecd_CellDistinguisher(
  bulkData$Indata$T,
  genesymb = rownames(bulkData$Indata$T),
  numCellClasses = dim(bulkData$Indata$C_ref)[2],
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
    as.matrix(bulkData$Indata$T),
    result$CellDistinguisher$bestDistinguishers,
    nonNegativeOnly = T,
    convexSolution = T,
    verbose = 0),
    error = function(e) return(list(sampleCompositions =
                                      matrix(rep(1/dim(bulkData$Indata$C_ref)[2],
                                                 dim(bulkData$Indata$C_ref)[2]*ncol(bulkData$Indata$T)),
                                             ncol=dim(bulkData$Indata$C_ref)[2]))))
# library(CellMix)
# deconvolutionSSKL <- CellDistinguisher::gecd_DeconvolutionCellMix(
#   as.matrix(bulkData$Indata$T), 
#   result$CellDistinguisher$bestDistinguishers, 
#   method="ssKL", 
#   maxIter=dim(bulkData$Indata$C_ref)[2])

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
rownames(bulkData$Indata$P)
labels <- Row_label(t(result$CellDistinguisher$p),t(bulkData$Indata$P))
rownames(result$CellDistinguisher$p) <- labels
colnames(result$CellDistinguisher$c) <- labels
result$CellDistinguisher$p <- result$CellDistinguisher$p[rownames(bulkData$Indata$P),]
result$CellDistinguisher$c <- result$CellDistinguisher$c[,colnames(bulkData$Indata$C)]
result$CellDistinguisher$result$c = getPearsonRMSE(result$CellDistinguisher$c,bulkData$Indata$C)
result$CellDistinguisher$result$p = getPearsonRMSE(result$CellDistinguisher$p,bulkData$Indata$P)
result$CellDistinguisher$result$t = getPearsonRMSE(result$CellDistinguisher$t,bulkData$Indata$T)
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

saveRDS(result,"XXX.result")