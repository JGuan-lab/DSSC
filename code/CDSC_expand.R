#
#--------the first step of QC--------------
QC_first <- function(data, full_phenoData){#
  require(dplyr); require(Matrix)
  # ???????���:???Ŀ???С????��?????????庬��????????MAD??ϸ??????
  filterCells <- function(filterParam){
    cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
    cellsToRemove
  }
  libSizes <- colSums(data)
  gene_names <- rownames(data)
  mtID <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
  rbID <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)
  mtPercent <- colSums(data[mtID, ])/libSizes
  rbPercent <- colSums(data[rbID, ])/libSizes
  lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>% 
    unlist() %>% 
    unique() -> cellsToRemove
  if(length(cellsToRemove) != 0){
    data <- data[,-cellsToRemove]
    full_phenoData <- full_phenoData[-cellsToRemove,]#
  }
  # Keep only "detectable" genes: at least 5% of cells (regardless of the group) have a read/UMI count different from 0
  #ֻ???����ɼ????ġ?????:????5%??ϸ??(??????????һ??)?Ķ?/UMI??????ͬ??0
  keep <- which(Matrix::rowSums(data > 0) >= round(0.05 * ncol(data)))#8415
  data = data[keep,]
  
 list(data = data, phenoData=full_phenoData)
  
  return(list(data= data, full_phenoData= full_phenoData))
}

#--------the second step of QC--------------
QC_second <- function(data, full_phenoData, leastNum = 50){
  set.seed(24)
  require(limma); require(dplyr); require(pheatmap)
  original_cell_names = colnames(data)
  colnames(data) <- as.character(full_phenoData$cellType[match(colnames(data),full_phenoData$cellID)])
  # Keep CTs with >= 50 cells after QC
  
  cell_counts = table(colnames(data))
  to_keep = names(cell_counts)[cell_counts >= leastNum]
  pData <- full_phenoData[full_phenoData$cellType %in% to_keep,]
  to_keep = which(colnames(data) %in% to_keep)   
  data <- data[,to_keep]
  original_cell_names <- original_cell_names[to_keep]
  return(list(data= data, pData= pData, original_cell_names= original_cell_names, cell_counts= cell_counts))
}

#--------Split data into train and test--------------
Split_data <- function(data, pData,cell_counts){
  training <- as.numeric(unlist(sapply(unique(colnames(data)), function(x) {
    sample(which(colnames(data) %in% x), cell_counts[x]/2) })))
  testing <- which(!1:ncol(data) %in% training)
  
  # Generate phenodata for reference matrix C
  pData_train = pData[training,]
  pData_test = pData[testing,]
  train <- data[,training]
  test <- data[,testing]
  return(list(train=train, test=test, pData_train=pData_train, pData_test=pData_test, training=training, testing=testing))
}
Normalization <- function(data){
  
  data <- edgeR::DGEList(data)
  CtrlGenes <- grep("ERCC-",rownames(data))
  
  if(length(CtrlGenes)>1){
    
    spikes <- data[CtrlGenes,]
    spikes <- edgeR::calcNormFactors(spikes, method = "TMM") 
    data$samples$norm.factors <- spikes$samples$norm.factors
    
  } else {
    
    data <- edgeR::calcNormFactors(data, method = "TMM")  
    
  }
  
  return(data)
  
}
#--------Find marker gene based on limma------------
Find_markerGene_limma <- function(train,plotmarker = F,norm1 = "TMM",log2.threshold = 1){
  
  library(limma)
  library(dplyr)
  library(scater)
  if(norm1 %in% c("cpm","CPM")){
    train2 <- calculateCPM(train)  #CPM  standardization
  }else if(norm1 %in% c("tmm","TMM")){
    train2 <- Normalization(train)
  }else if(norm1 %in% c("none","NONE")){
    train2 <- train
  }else{
    print("ERROR NORMIAZTION,please input right normiaztion names")
    break;
  }
  
  # INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account 
  #[compare one group with average expression of all other groups]
  annotation = factor(colnames(train2))
  design <- model.matrix(~0+annotation)
  colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
  dim(design)# [1] cellNum  * cellTypeNum
  # design is matrix contain 0 or 1,meaning which ct she gene is 
  
  cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
  colnames(cont.matrix) <- colnames(design)
  diag(cont.matrix) <- (ncol(design)-1)/ncol(design)
  # cont.matrix is Symmetrical matrix,which diag is 0.9,others is -0.1
  
  v <- limma::voom(train2, design=design, plot=plotmarker) 
  # voom 
  fit <- limma::lmFit(v, design)
  # lmFit
  # head(coef(fit))
  fit2 <- limma::contrasts.fit(fit, cont.matrix)# 
  # contrasts.fit?Ƚ?ÿ??????
  fit2 <- limma::eBayes(fit2, trend=TRUE)
  # ??׼?????ľ??鱴Ҷ˹ƽ?????????????????ı?׼????????С?ö??ı?׼??????С??ƽ????׼???
  # plotSA(fit2, main="Final model: Mean-variance trend")
  # ??ͼ
  
  # topTable ?г????????????? 
  # tT=topTable(fit2,adjust='BH',coef=1:ncol(fit2$coefficients),number=Inf,p.value=1) ##p.value?Լ?????
  tT=topTable(fit2,adjust='fdr',coef=1,number=Inf,p.value=1) ##p.value?Լ?????
  tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
  colnames(tT)=c("FDR","P.Value","logFC")
  tT$GENEs=rownames(tT)
  
  markers = marker.fc(fit2, log2.threshold)
  tT <- tT[markers$gene ,]
  markers$FDR <- tT$FDR
  markers$P.value <- tT$P.Value
  
  return(markers)
}

#--------Compose marker matrix------------
marker.fc <- function(fit2, log2.threshold = 1, output_name = "markers"){
  # log2FC > 1 && p.values < 0.1 && FDR < 0.1
  topTable_RESULTS = limma::topTable(fit2, coef = 1:ncol(fit2$coefficients), number = Inf, adjust.method = "BH", p.value = 0.1, lfc = log2.threshold)
  # FDR ???㲢??ȡ
  AveExpr_pval <- topTable_RESULTS[,(ncol(topTable_RESULTS)-3):ncol(topTable_RESULTS)]#??ȡ????????
  topTable_RESULTS <- topTable_RESULTS[,1:(ncol(topTable_RESULTS)-4)]#??ȡǰcellTypeNum??
  
  if(length(grep("ERCC-",topTable_RESULTS$gene)) > 0){ 
    topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-",topTable_RESULTS$gene),] }
  
  markers <- apply(topTable_RESULTS,1,function(x){
    temp = sort(x)
    ((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS)-1]) >= log2.threshold) | (abs(temp[1] - temp[2]) >= log2.threshold)
    # ????????????ǰ��???Ĳ??��???Ĳ????log2.threshold
  })
  topTable_RESULTS = topTable_RESULTS[markers,]
  
  markers <- cbind.data.frame(rownames(topTable_RESULTS),
                              t(apply(topTable_RESULTS, 1, function(x){
                                temp = max(x)
                                if(temp < log2.threshold){
                                  #????tempС??log2.threshold??ȡ??С??
                                  temp = c(min(x),colnames(topTable_RESULTS)[which.min(x)])
                                } else {
                                  #????temp????log2.threshold??ȡ??????
                                  temp = c(max(x),colnames(topTable_RESULTS)[which.max(x)])
                                } 
                                temp
                              })))
  
  colnames(markers) <- c("gene","log2FC","CT")
  markers$log2FC = as.numeric(as.character(markers$log2FC))
  markers <- markers %>% dplyr::arrange(CT,desc(log2FC)) 
  
  markers$AveExpr <- AveExpr_pval$AveExpr[match(markers$gene,rownames(AveExpr_pval))]
  markers$gene <- as.character(markers$gene)
  markers$CT <- as.character(markers$CT)
  
  #write.table(markers, file = output_name, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  return(markers)
}


#--------Compose C matrix genecelltype------------
Get_C <- function(train){
  # reference matrix (C) + refProfiles.var from TRAINING dataset
  cellType <- colnames(train)
  group = list()
  for(i in unique(cellType)){ 
    group[[i]] <- which(cellType %in% i)
  }
  # C_ref = lapply(group,function(x) Matrix::rowMeans(train[,x])) 
  #C should be made with the mean (not sum) to agree with the way markers were selected
  C_ref = lapply(group,function(x) Matrix::rowMeans(train[,x])) 
  C_ref = round(do.call(cbind.data.frame, C_ref),4)
  
  refProfiles.var = lapply(group,function(x) train[,x])
  refProfiles.var = lapply(refProfiles.var, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
  refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var))
  rownames(refProfiles.var) <- rownames(train)
  return(list(C_ref=C_ref, refProfiles.var=refProfiles.var))
}

#--------Compose T、P matrix ------------
Generator <- function(sce, phenoData, Num.mixtures = 1000, pool.size = 100, min.percentage = 1, max.percentage = 99, seed = 24){ 
  
  CT = unique(phenoData$cellType)
  ?stopifnot(length(CT) >= 2)
  
  set.seed(seed)
  require(dplyr)
  require(gtools)
  
  cell.distribution = data.frame(table(phenoData$cellType),stringsAsFactors = FALSE) 
  colnames(cell.distribution) = c("CT","max.n")
  
  Tissues = list()
  Proportions = list()
  
  for(y in 1:Num.mixtures){
    
    #Only allow feasible mixtures based on cell distribution
    while(!exists("P")){
      
      num.CT.mixture = sample(x = 2:length(CT),1)
      selected.CT = sample(CT, num.CT.mixture, replace = FALSE)
      
      P = runif(num.CT.mixture, min.percentage, max.percentage) 
      P = round(P/sum(P), digits = log10(pool.size))  #sum to 1
      P = data.frame(CT = selected.CT, expected = P, stringsAsFactors = FALSE)
      
      missing.CT = CT[!CT %in% selected.CT]
      missing.CT = data.frame(CT = missing.CT, expected = rep(0, length(missing.CT)), stringsAsFactors = FALSE)
      
      P = rbind.data.frame(P, missing.CT)
      potential.mix = merge(P, cell.distribution)
      potential.mix$size = potential.mix$expected * pool.size
      
      if( !all(potential.mix$max.n >= potential.mix$size) | sum(P$expected) != 1){
        rm(list="P") 
      }
      
    }
    
    # Using info in P to build T simultaneously
    chosen_cells <- sapply(which(P$expected != 0), function(x){
      
      n.cells = P$expected[x] * pool.size
      chosen = sample(phenoData$cellID[phenoData$cellType == P$CT[x]],
                      n.cells)
      
      chosen
    }) %>% unlist()
    
    
    T <- Matrix::rowSums(sce[,colnames(sce) %in% chosen_cells]) %>% as.data.frame()
    colnames(T) = paste("mix",y,sep="")
    
    P = P[,c("CT","expected")]
    P$mix = paste("mix",y,sep="")
    
    Tissues[[y]] <- T
    Proportions[[y]] <- P
    
    rm(list=c("T","P","chosen_cells","missing.CT"))
    
  }
  P = do.call(rbind.data.frame, Proportions)
  T = do.call(cbind.data.frame, Tissues)
  
  P = data.table::dcast(P, CT ~ mix, 
                        value.var = "expected",
                        fun.aggregate = sum) %>% data.frame(.,row.names = 1) 
  P = P[,gtools::mixedsort(colnames(P))]
  return(list(T = T, P = P))
  
} 

#-------shift markers based on requirement-----
marker_strategies <- function(marker_distrib, marker_strategy, C){
  set.seed(4)
  
  if(marker_strategy == "all"){
    
    #using all markers that were found
    markers = marker_distrib
    
  } else if (marker_strategy == "pos_fc"){
    
    # using only markers with positive FC (=over-expressed in cell type of interest)
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% as.data.frame()
    
  } else if (marker_strategy == "top_25p_logFC"){
    
    # top 50% of markers (per CT) based on logFC
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(ceiling(n()*0.25), wt=log2FC) %>% as.data.frame()
  
  } else if (marker_strategy == "top_50p_logFC"){
    
    # top 50% of markers (per CT) based on logFC
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(ceiling(n()*0.5), wt=log2FC) %>% as.data.frame()
  
  } else if (marker_strategy == "top_75p_logFC"){
    
    # top 50% of markers (per CT) based on logFC
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(ceiling(n()*0.75), wt=log2FC) %>% as.data.frame()
    
  } else if (marker_strategy == "bottom_50p_logFC"){
    
    # bottom 50% of markers based on logFC
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(floor(n()*-0.5), wt=log2FC) %>% as.data.frame()
    
  } else if (marker_strategy == "top_50p_AveExpr"){
    
    # top 50% of markers based on average gene expression (baseline expression)
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(AveExpr)) %>% group_by(CT) %>% dplyr::top_n(ceiling(n()*0.5), wt=log2FC) %>% as.data.frame()
    
  } else if (marker_strategy == "bottom_50p_AveExpr"){
    
    # low 50% based on average gene expression.
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(AveExpr)) %>% group_by(CT) %>% dplyr::top_n(floor(n()*-0.5), wt=log2FC) %>% as.data.frame()
    
  } else if (marker_strategy == "top_n2"){
    
    # using the top 2 genes/CT with highest log2FC
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(2, wt=log2FC) %>% as.data.frame()

  } else if (marker_strategy == "top_n1"){
    
    # using the top 3 genes/CT with highest log2FC
    markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(1, wt=log2FC) %>% as.data.frame()
    
  } else if (marker_strategy == "random5"){
    
    # using 5 random markers for each different cell types
    markers = marker_distrib[1:(ncol(C)*5),]
    markers$CT = rep(colnames(C),5) #labelling purposes: important for semi-supervised
    markers$gene = sample(rownames(C), nrow(markers), replace = FALSE)
    
  }
  
  return(markers) 
  
}

#-------combined two martix------------
combine_2 <- function(RESULTS,P){
  
  RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)),]
  RESULTS = reshape2::melt(RESULTS)
  colnames(RESULTS) <-c("CT","tissue","observed_values")
  
  if(!is.null(P)){
    P = as.matrix(P)
    P = P[gtools::mixedsort(rownames(P)),]
    #P$CT = rownames(P)
    P = reshape2::melt(P)#, id.vars="CT"
    colnames(P) <-c("CT","tissue","expected_values")
    
    RESULTS = merge(RESULTS,P)
    RESULTS$expected_values <-round(RESULTS$expected_values,3)
    RESULTS$observed_values <-round(RESULTS$observed_values,3)
  }
  return(RESULTS)
}

#-------calculate the RMSE and Pearson
getPearsonRMSE <- function(RESULTS,P){
  library(dplyr)
  RESULTS = combine_2(RESULTS,P)
  RESULTS = RESULTS %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4),Pearson=cor(observed_values,expected_values) %>% round(.,4))
  return(RESULTS)
}
Row_label_location <- function(num,X,k){
  if( num %% k == 0 ){
    location_y = num/(dim(X)[1])
    location_x = k
  }else{
    location_y = ceiling(num/(dim(X)[1])) 
    location_x = num - floor(num/dim(X)[1])*k
  }
  location = location_x
  location[2] = location_y
  return(location)
}

#-------find labels when the num of labels above 3-----------
Row_label <- function(result_c, C_ref, leastnum){
  X <- cor(result_c,C_ref)
  k = ncol(C_ref)
  if(ncol(C_ref) >= 4){
    num = k- 3
    Row_label = matrix(1:k,1,k)
    colnames_X = colnames(X)
    for (j in 1:num) {#先按大小赋值前k-3个最后3个再穷举
      # location <- Row_label_location(which(X==max(X),arr.ind=T),X,k)#????????
      location <- which(X==max(X),arr.ind=T)#????????
      Row_label[location[1]] = colnames_X[location[2]]
      X[location[1],] = -1
      X[,location[2]] = -1
    }
    Row_label_all = NULL
    vacant_y_1 <- colnames(X)[which(!( colnames(X) %in% Row_label))]
    vacant_x_1 <- which(!(Row_label %in% colnames(X)))
  }else if(ncol(C_ref) <= 3){
    #细胞类型数目大于等于1，小于等于三
    Row_label = matrix(1:k,1,k)
    Row_label_all = NULL
    vacant_y_1 <- colnames(X)
    vacant_x_1 <- 1:k
  }
  maybe = 1
  for (j in 1:length(vacant_x_1)){
    Row_label_ <- Row_label
    Row_label_[vacant_x_1[1]] = vacant_y_1[j]
    
    vacant_y_2 <- vacant_y_1[-j]
    vacant_x_2 <- which(!(Row_label_ %in% colnames(X)))
    for (kk in 1:length(vacant_x_2)) {
      Row_label_all[[  maybe ]] = Row_label_
      Row_label_all[[  maybe ]][ vacant_x_2[kk] ] = vacant_y_2[1]
      Row_label_all[[  maybe ]][ vacant_x_2[which(!(vacant_x_2 == vacant_x_2[kk]))] ] = vacant_y_2[2]
      maybe = maybe + 1
    }
  }

  ans = NULL
  for (i in 1:length(Row_label_all)) {
    XX= result_c
    colnames(XX) <- Row_label_all[[i]]
    XX <- XX[,colnames(C_ref) ]
    pear <- cor(XX,C_ref)
    pear[is.na(pear)] <- 0
    ans[i] = sum(diag(pear))
  }
  i = which.max(ans)
  Row_label_end = Row_label_all[[i]]
  #?Ƚ????Կ??ܽ????? ????ֵ???õ????????ŵ?label????????Row_label_end
  

  return(Row_label_end)
}
#------get cor matrix----

GetCorMatrix <- function(a,a_real,a_ref = NULL,leastnum = 3,matrix){
  if(matrix == "p"){
    ctlabels <- Row_label(t(a),t(a_real),leastnum)
    rownames(a) <- ctlabels
    ctlabels_ <- rownames(a_real)
    a <- a[ctlabels_, ]
    return(cor(t(a),t(a_real)))
  }else if(matrix == "c"){
    if(is.null(a_ref)){a_ref = a_real}
    ctlabels <- Row_label(a,a_ref,leastnum)
    colnames(a) <- ctlabels
    ctlabels_ <- intersect(colnames(a_real),as.character(ctlabels))
    a <- a[ ,ctlabels_]
    a_real = as.matrix(a_real[rownames(a),ctlabels_])
    
    return(cor(a,a_real))
  }else{
    print("please inpute matrix: p or c")
  }
}

GetCorMatrixLM22 <- function(a,a_real,a_ref = NULL,leastnum = 3,matrix){
  if(matrix == "p"){
    
    ctlabels <- Row_label(t(a),t(a_real),leastnum)
    rownames(a) <- ctlabels
    a_real = MergeCellType(a_real,'p')
    ctlabels_ <- rownames(a_real)
    a <- a[ctlabels_, ]
    return(cor(t(a),t(a_real)))
  }else if(matrix == "c"){
    if(is.null(a_ref)){a_ref = a_real}
    ctlabels <- Row_label(a,a_ref,leastnum)
    colnames(a) <- ctlabels
    a = MergeCellType(a,'c')
    ctlabels_ <- intersect(colnames(a_real), colnames(a))
    a <- a[ ,ctlabels_]
    a_real = as.matrix(a_real[rownames(a),ctlabels_])
    
    return(cor(a,a_real))
  }else{
    print("please inpute matrix: p or c")
  }
}
#------get sample cor----
GetSampleCor <- function(a,a_real,leastnum = 3,matrix){
  if(matrix == "p"){
    ctlabels <- Row_label(t(a),t(a_real),leastnum)
    rownames(a) <- ctlabels
    ctlabels_ <- rownames(a_real)
    a <- a[ctlabels_, ]
    return(diag(cor(a,a_real)))
  }else if(matrix == "c"){
    ctlabels <- Row_label(a,a_real,leastnum)
    colnames(a) <- ctlabels
    ctlabels_ <- colnames(a_real)
    a <- a[ ,ctlabels_]
    return(diag(cor(t(a),t(a_real))))
  }else{
    print("please inpute matrix: p or c")
  }
}
#--------calculate result-------------
calculate_result <- function(result_c,result_p,
                             T,C,C_ref,P, 
                             lambda1=NULL, lambda2=NULL, lambdaC=NULL,
                             number_iter=NULL,seedd=NULL,TerCondition=NULL,
                             leastnum = 3,ctlabels=NULL){
  result <- NULL
  if(is.null(ctlabels)){
      if(!is.null(C_ref)){
        ctlabels <- Row_label(result_c,C_ref,leastnum)
    }else{
        ctlabels <- Row_label(t(result_p),t(P),leastnum)
    }}
  
  rownames(result_p) <- ctlabels
  colnames(result_c) <- ctlabels
  ctlabels_ <- intersect(rownames(P),as.character(ctlabels))
  result_p <- result_p[ctlabels_, ]
  result_c <- result_c[,ctlabels_]
  
  C = as.matrix(C[rownames(result_c),colnames(result_c)])
  P = as.matrix(P[rownames(result_p),colnames(result_p)])
  library(dplyr)
  RESULT_c_all = combine_2( C,result_c)
  RESULT_c_all = RESULT_c_all %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4),Pearson=cor(observed_values,expected_values) %>% round(.,4))
  
  RESULT_p_all = combine_2( P,NORM(result_p))
  RESULT_p_all = RESULT_p_all %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4),Pearson=cor(observed_values,expected_values) %>% round(.,4))
  
  RESULT_t_all = combine_2( T,result_c%*% result_p)
  RESULT_t_all = RESULT_t_all %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4),Pearson=cor(observed_values,expected_values) %>% round(.,4))
  result = NULL
  result = cbind(RESULT_p_all,RESULT_c_all)
  result = cbind(result,RESULT_t_all)
  
  if(!is.null(lambda1)){
    result = cbind(lambdaC,result)
    result = cbind(lambda2,result)
    result = cbind(lambda1,result)
    result = cbind(result,number_iter)
    result = cbind(result,seedd)
    result = cbind(result,TerCondition)
    colnames(result) = c("lambda1", "lambda2","lambdaC","RMSE_to_P","Peason_to_P","RMSE_to_C","Peason_to_C","RMSE_to_T","Peason_to_T","number_iter","seedd","TerCondition")
  }else{
    colnames(result) = c("RMSE_to_P","Peason_to_P","RMSE_to_C","Peason_to_C","RMSE_to_T","Peason_to_T")
  }
  
  
  return(result)
}
#--------calculate result for PBMCs-------------
calculate_pbmcs_result <- function(result_c,result_p,
                             T,C,C_ref,P, 
                             lambda1=NULL, lambda2=NULL, lambdaC=NULL,
                             number_iter=NULL,seedd=NULL,TerCondition=NULL,
                             leastnum = 3,ctlabels=NULL){
  result <- NULL
  if(is.null(ctlabels)){
    if(!is.null(C_ref)){
      ctlabels <- Row_label(result_c,C_ref,leastnum)
    }else{
      ctlabels <- Row_label(t(result_p),t(P),leastnum)
    }}
  
  rownames(result_p) <- ctlabels
  colnames(result_c) <- ctlabels
  result_p = MergeCellType(result_p,'p')  
  result_c = MergeCellType(result_c,'c')

  ctlabels_ <- rownames(result_p)
  result_p <- result_p[ctlabels_, ]
  result_c <- result_c[,ctlabels_]
  
  C = as.matrix(C[rownames(result_c),colnames(result_c)])
  P = as.matrix(P[rownames(result_p),colnames(result_p)])
  library(dplyr)
  RESULT_c_all = combine_2( C,result_c)
  RESULT_c_all = RESULT_c_all %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4),Pearson=cor(observed_values,expected_values) %>% round(.,4))
  
  RESULT_p_all = combine_2( SumEqual_1(P), SumEqual_1(result_p))
  RESULT_p_all = RESULT_p_all %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4),Pearson=cor(observed_values,expected_values) %>% round(.,4))
  
  RESULT_t_all = combine_2( T,result_c%*% result_p)
  RESULT_t_all = RESULT_t_all %>% dplyr::summarise(RMSE = sqrt(mean((observed_values-expected_values)^2)) %>% round(.,4),Pearson=cor(observed_values,expected_values) %>% round(.,4))
  result = NULL
  result = cbind(RESULT_p_all,RESULT_c_all)
  result = cbind(result,RESULT_t_all)
  
  if(!is.null(lambda1)){
    result = cbind(lambdaC,result)
    result = cbind(lambda2,result)
    result = cbind(lambda1,result)
    result = cbind(result,number_iter)
    result = cbind(result,seedd)
    result = cbind(result,TerCondition)
    colnames(result) = c("lambda1", "lambda2","lambdaC","RMSE_to_P","Peason_to_P","RMSE_to_C","Peason_to_C","RMSE_to_T","Peason_to_T","number_iter","seedd","TerCondition")
  }else{
    colnames(result) = c("RMSE_to_P","Peason_to_P","RMSE_to_C","Peason_to_C","RMSE_to_T","Peason_to_T")
  }
  return(result)
}

MergeCellType <- function(result_p,matrix = 'p'){
  if(matrix %in% c('p', 'P')){
    names_ <- colnames(result_p)
    result_p = apply(result_p,2,function(x) x/sum(x))
    result_pp <- as.data.frame(t(result_p))
  }else if(matrix %in% c('c', 'C')){
    names_ <- rownames(result_p)
    result_pp <- as.data.frame(result_p)
    
  }else{
    print("please input right matrix")
    return()
  }
    result_p_new <- NULL
    result_p_new$B.cells <- (result_pp$B.cells.naive + result_pp$B.cells.memory)
    result_p_new$T.cells.CD8 <- (result_pp$T.cells.CD8 )
    result_p_new$T.cells.CD4 <- (result_pp$T.cells.CD4.naive + result_pp$T.cells.CD4.memory.resting +result_pp$T.cells.CD4.memory.activated)
    result_p_new$NK.cells <- (result_pp$NK.cells.resting + result_pp$NK.cells.activated)
    result_p_new$Monocytes <- (result_pp$Monocytes)
    if(matrix %in% c('p', 'P')){
      result_p <- as.matrix(do.call(rbind, result_p_new))
      colnames(result_p) <- names_
    }else if(matrix %in% c('c', 'C')){
      result_p <-  as.matrix(t(do.call(rbind, result_p_new)))
      rownames(result_p) <- names_
    }
    return(result_p)
}
#-------------simulate-----------
QCscDATA <- function(scData, leastNum=50){
  #this function donnot split data, but simulate T、P
  #------------QC_1------------
  table(scData$full_phenoData$cellType)
  scData_QC <- QC_first(scData$data, scData$full_phenoData)
  #------------QC_2------------
  table(scData_QC$full_phenoData$cellType)
  scData_QC <- QC_second(scData_QC$data, scData_QC$full_phenoData, leastNum)

  return(scData_QC)
}

#-------------simulate C-----------
scSimulateC <- function(scData, leastNum=50, plotmarker = F,
                        norm1 = "CPM",log2.threshold=log2(2)){
  #this function donnot split data, but simulate T、P
  # ------------QC_1------------
  table(scData$full_phenoData$cellType)
  scData_simualte <- QC_first(scData$data, scData$full_phenoData)
  # ------------QC_2------------
  table(scData_simualte$full_phenoData$cellType)
  scData_simualte <- QC_second(scData_simualte$data, scData_simualte$full_phenoData, leastNum)
  #------------find marker genes----------
  scData_simualte$markers <- Find_markerGene_limma(scData_simualte$data, 
                                                   plotmarker = F,
                                                   norm1,
                                                   log2.threshold)
  
  scData_simualte$df = which(scData_simualte$markers$FDR <= 0.1 & scData_simualte$markers$P.value<0.1)
  scData_simualte$markerslist = scData_simualte$markers[scData_simualte$df,]
  # scData_simualte$markerslist = scData_simualte$markers
  #------------GET c----------
  scData_simualte$C <- Get_C(scData_simualte$data)[[1]]
  scData_simualte$refProfiles.var <- Get_C(scData_simualte$data)[[2]]
  # # use test data to cpmposite C
  # require(dplyr);
  # colnames(scData_simualte$data) <- scData_simualte$original_cell_names
  # 
  # generator <- Generator(sce = scData_simualte$data, phenoData = scData_simualte$pData, Num.mixtures = 1000, pool.size = 100)
  # scData_simualte$T <- generator[["T"]]
  # scData_simualte$P <- generator[["P"]]
  return(scData_simualte)
}

#----------simulate_1sc data-----------
scSimulateSplit <- function(scData, leastNum=50, plotmarker = F, 
                            Num.mixtures = 1000, pool.size = 100,
                            norm1 = "CPM",log2.threshold){
  #------------QC_1------------
  table(scData$full_phenoData$cellType)
  scData_simualte <- QC_first(scData$data, scData$full_phenoData)
  #------------QC_2------------
  table(scData_simualte$full_phenoData$cellType)
  scData_simualte <- QC_second(scData_simualte$data, scData_simualte$full_phenoData, leastNum)
  #-----------split data----------
  # Data split into train & test  
  scData_simualte$split <- Split_data(scData_simualte$data, scData_simualte$pData, scData_simualte$cell_counts)
  #------------find marker genes----------
  scData_simualte$markers <- Find_markerGene_limma(scData_simualte$split$train, plotmarker,norm1,log2.threshold)
  scData_simualte$df = which(scData_simualte$markers$FDR <= 0.1 & scData_simualte$markers$P.value<0.1)
  scData_simualte$markerslist = scData_simualte$markers[scData_simualte$df,]
  # scData_simualte$markerslist = scData_simualte$markers
  #------------GET c----------
  scData_simualte$C_ref <- Get_C(scData_simualte$split$train)[[1]]
  scData_simualte$refProfiles.var <- Get_C(scData_simualte$split$train)[[2]]
  scData_simualte$C <- Get_C(scData_simualte$split$test)[[1]]
  # use test data to cpmposite C
  require(dplyr);
  colnames(scData_simualte$split$test) <- scData_simualte$original_cell_names[scData_simualte$split$testing]
  
  generator <- Generator(sce = scData_simualte$split$test, phenoData = scData_simualte$split$pData_test, 
                         Num.mixtures, pool.size)
  scData_simualte$T <- generator[["T"]]
  scData_simualte$P <- generator[["P"]]
  return(scData_simualte)
}

# scData = camp
scSimulateShift <- function(scData, marker_strategy = "all",standardization = TRUE ){
  # shift gene
  dim(scData$simulate1$T)
  dim(scData$simulate1$C)
  all(intersect(rownames(scData$simulate1$T), rownames(scData$simulate1$markerslist))==
        intersect(rownames(scData$simulate1$C), rownames(scData$simulate1$markerslist)))
  all(rownames(scData$simulate1$C_ref) == rownames(scData$simulate1$T))
  all(rownames(scData$simulate1$C_ref) == rownames(scData$simulate1$C))
  # save data with all or half gene
  cat('the number of marker genes: ' ,dim(scData$simulate1$markerslist)[1] )
  # all, top_50p_logFC
  scData$simulate1$shif_marker = marker_strategies(scData$simulate1$markerslist, marker_strategy, scData$simulate1$C_ref)
  cat('\nthe number of shifed marker genes: ' , dim(scData$simulate1$shif_marker)[1] )
  
  scData$Indata$gene_ <- scData$simulate1$shif_marker$gene
  length(scData$Indata$gene_)
  scData$Indata$T = scData$simulate1$T[scData$Indata$gene_,]
  scData$Indata$C = scData$simulate1$C[scData$Indata$gene_,]
  scData$Indata$C_ref = scData$simulate1$C_ref[scData$Indata$gene_,]
  scData$Indata$refProfiles.var = scData$simulate1$refProfiles.var[scData$Indata$gene_,]
  scData$Indata$P = scData$simulate1$P
  
  # shift celltype
  colnames(scData$Indata$C)
  colnames(scData$Indata$C_ref)
  rownames(scData$Indata$P)
  
  unique(scData$simulate1$markerslist$CT)
  scData$Indata$Cellltype_name <- unique(scData$simulate1$markerslist$CT)
  scData$Indata$C = scData$Indata$C[,scData$Indata$Cellltype_name]
  scData$Indata$C_ref = scData$Indata$C_ref[,scData$Indata$Cellltype_name]
  scData$Indata$P = scData$Indata$P[scData$Indata$Cellltype_name,]
  scData$Indata$refProfiles.var = scData$Indata$refProfiles.var[,scData$Indata$Cellltype_name]
  # scData$Indata$PData_train = scData$simulate1$split$pData_train
  
  scData$Indata$P = as.matrix(scData$Indata$P)
  scData$Indata$T = as.matrix(scData$Indata$T)
  scData$Indata$C = as.matrix(scData$Indata$C)
  scData$Indata$C_ref = as.matrix(scData$Indata$C_ref)
  
  #------------standardization----------
  library(scater)
  if(standardization){
    scData$Indata$T <- calculateCPM(scData$Indata$T)
    scData$Indata$C <- calculateCPM(scData$Indata$C)
    scData$Indata$C_ref <- calculateCPM(scData$Indata$C_ref)
  }

  length(which(colSums(scData$Indata$P) == 1))
  return(scData)
}

BulkMixturesShift <- function(bulkData, marker_strategy = "all"){
  # shift gene
  dim(bulkData$T)
  dim(bulkData$C)
  all(rownames(bulkData$C_ref) == rownames(bulkData$T))
  # save data with all or half gene
  cat('the number of marker genes: ' ,dim(bulkData$markerslist)[1] )
  bulkData$shif_marker = marker_strategies(bulkData$markerslist, marker_strategy="all", bulkData$C_ref)
  cat('\nthe number of shifed marker genes: ' , dim(bulkData$shif_marker)[1] )
  
  bulkData$Indata$gene_ <- bulkData$shif_marker$gene
  bulkData$Indata$T = bulkData$T[bulkData$Indata$gene_,]
  bulkData$Indata$C_ref = bulkData$C_ref[bulkData$Indata$gene_,]
  bulkData$Indata$refProfiles.var = bulkData$refProfiles.var[bulkData$Indata$gene_,]
  bulkData$Indata$P = bulkData$P
  
  # shift celltype
  colnames(bulkData$Indata$C_ref)
  rownames(bulkData$Indata$P)
  unique(bulkData$markerslist$CT)
  bulkData$Indata$Cellltype_name <- rownames(bulkData$Indata$P)
  bulkData$Indata$C_ref = bulkData$Indata$C_ref[,bulkData$Indata$Cellltype_name]
  bulkData$Indata$P = bulkData$Indata$P[bulkData$Indata$Cellltype_name,]
  bulkData$Indata$refProfiles.var = bulkData$Indata$refProfiles.var[,bulkData$Indata$Cellltype_name]
  
  bulkData$Indata$P = as.matrix(bulkData$Indata$P)
  bulkData$Indata$T = as.matrix(bulkData$Indata$T)
  bulkData$Indata$C_ref = as.matrix(bulkData$Indata$C_ref)
  
  #------------standardization----------
  library(scater)
  bulkData$Indata$T <- calculateCPM(bulkData$Indata$T)
  bulkData$Indata$C_ref <- calculateCPM(bulkData$Indata$C_ref)
  length(which(colSums(bulkData$Indata$P) == 1))
  bulkData$Indata$C =bulkData$Indata$C_ref
  return(bulkData)
}

Combine_2scdata_dream <- function(scData_1,scData_2,leastNum1=50,leastNum2=50, plotmarker = F,
                                  Num.mixtures = 1000, pool.size = 100,
                                  norm1 = "CPM", norm2 = "CPM",log2.threshold=1){
  #' @param scDATA_1是测试集，用于仿真生成bulk矩阵
  #' @param scDATA_2是训练集，用于产生参考C和求取markerGene
  #' @param leastNum1是测试集QC时保留的细胞类型的最小细胞数目
  #' @param leastNum2是训练集QC时保留的细胞类型的最小细胞数目
  #' @param plotmarker是训练集寻找markerGene时是否画拟合图像
  #' @param Num.mixtures是最终的bulkData的样本数目，也是仿真的抽样次数
  #' @param pool.size是每次仿真抽样的容积
  #' @param norm1对测试集结果标准化
  #' @param norm2对训练集结果标准化
  ComDate=NULL
  # 质量控制
  ComDate <- list(
    scData_1 = scData_1,
    Test = QCscDATA(scData_1,leastNum1),
    scData_2 = scData_2,
    Train = QCscDATA(scData_2,leastNum2)
    )
  
  # 理想情况下匹配gene和celltype
  cat( "The number of genes matching TrainSet and TestSet:",
         length(intersect(rownames(ComDate$Test$data), rownames(ComDate$Train$data))))
  cat( "\nCelltypes matching TrainSet and TestSet:",
       intersect(unique(ComDate$Test$pData$cellType), unique(ComDate$Train$pData$cellType)))
  keepTestGene <- which(rownames(ComDate$Test$data) %in% 
                        intersect(rownames(ComDate$Test$data), rownames(ComDate$Train$data)))
  keepTrainGene <- which(rownames(ComDate$Train$data) %in% 
                        intersect(rownames(ComDate$Test$data), rownames(ComDate$Train$data)))
  keepTestCell <- which(ComDate$Test$pData$cellType %in% 
                      intersect(unique(ComDate$Test$pData$cellType), unique(ComDate$Train$pData$cellType)))
  keepTrainCell <- which(ComDate$Train$pData$cellType %in% 
                      intersect(unique(ComDate$Test$pData$cellType), unique(ComDate$Train$pData$cellType)))
  ComDate$Test$data <- ComDate$Test$data[keepTestGene,keepTestCell]
  ComDate$Test$pData <- ComDate$Test$pData[keepTestCell,]
  ComDate$Test$original_cell_names <- ComDate$Test$original_cell_names[keepTestCell]

  ComDate$Train$data <- ComDate$Train$data[keepTrainGene,keepTrainCell]
  ComDate$Train$pData <- ComDate$Train$pData[keepTrainCell,]
  ComDate$Train$original_cell_names <- ComDate$Train$original_cell_names[keepTrainCell]
  
  # 寻找marker Gene
  ComDate$markers <- Find_markerGene_limma(ComDate$Train$data, plotmarker,norm2,log2.threshold)
  cat("\nDifferent celltypes's marker Gene number of RAW:");print(sort(table(ComDate$markers$CT)))
  ComDate$df = which(ComDate$markers$FDR <= 0.1 & ComDate$markers$P.value<0.1)
  ComDate$markerslist = ComDate$markers[ComDate$df,]
  # ComDate$markerslist = ComDate$markers
  cat("\nDifferent celltypes's marker Gene number after shift:");print(sort(table(ComDate$markerslist$CT)))
  
  # GET matrix c(reference from train and real from test)
  ComDate$C_ref <- Get_C(ComDate$Train$data)[['C_ref']]
  ComDate$refProfiles.var <- Get_C(ComDate$Train$data)[['refProfiles.var']]
  ComDate$C <- Get_C(ComDate$Test$data)[[1]]
  
  # simulate T,P from test data
  require(dplyr);
  colnames(ComDate$Test$data) <- ComDate$Test$original_cell_names
  generator <- Generator(sce = ComDate$Test$data, phenoData = ComDate$Test$pData, Num.mixtures, pool.size)
  ComDate$T <- generator[["T"]]
  ComDate$P <- generator[["P"]]
  
  # for T,C_ref,refProfiles.var: shift gene and ct(when markers list ct != markers ct)
  gene_ <- ComDate$markerslist$gene
  ct_ <- unique(ComDate$markerslist$CT)
  ComDate$indata$T <- ComDate$T[gene_,]
  ComDate$indata$C_ref <- ComDate$C_ref[gene_,ct_]
  ComDate$indata$refProfiles.var <- ComDate$refProfiles.var[gene_,ct_]
  
  # for C P :shift gene and ct(when markers list ct != markers ct)
  ComDate$indata$C <- ComDate$C[gene_,ct_]
  ComDate$indata$P <- ComDate$P[ct_,]
  
  # for T,P :shift samples,,,,,,没有涉及到样本情况的改变，事实上
  # ComDate$indata$P <- NORM(ComDate$indata$P)
  # cat("\nReal P's sample numbers after shift ct:",length(which(colSums(ComDate$indata$P) != 0)))
  # ComDate$indata$P <- ComDate$indata$P[, which(colSums(ComDate$indata$P) != 0)]
  # ComDate$indata$T <- ComDate$indata$T[, colnames(ComDate$P)]

  # CPM标准化
  library(scater)
  if(norm1 %in% c("cpm","CPM")){
      ComDate$indata$T <- as.matrix(calculateCPM(ComDate$indata$T))
      ComDate$indata$C <- as.matrix(calculateCPM(ComDate$indata$C))
  }else{
    ComDate$indata$T <- as.matrix(ComDate$indata$T)
    ComDate$indata$C <- as.matrix(ComDate$indata$C)
  }
  if(norm2 %in% c("cpm","CPM")){
    ComDate$indata$C_ref <- as.matrix(calculateCPM(ComDate$indata$C_ref))
  }else{
    ComDate$indata$C_ref <- as.matrix(ComDate$indata$C_ref)
  }
   ComDate$indata$P <- as.matrix(ComDate$indata$P)
 
  
  return(ComDate)
}


CutCTofscDATA <- function(data,cutCT=c("unknown")){
  library(dplyr)
  if(cutCT %in% unique(data$pData$cellType)){
    # 存在需要删除的细胞类型
    keepCell <- which(data$pData$cellType != cutCT);length(keepCell)
    data$pData <- data$pData[keepCell,]
    data$data <- data$data[,keepCell]
    data$original_cell_names <- data$original_cell_names[keepCell]
    data$cell_counts <- table(data$pData$cellType)
  }
  return(data)
}

# scData_1 = Macosko
# scData_2 = Shekhar
# leastNum1=50
# leastNum2=50
# plotmarker = F
# Num.mixtures = 1000
# pool.size = 100
Combine_2scdata_hard <- function(scData_1,scData_2,leastNum1=50,leastNum2=50, plotmarker = F,
                                 Num.mixtures = 1000, pool.size = 100,cutCT=c("unknown"),
                                 norm1 = 'CPM', norm2 = 'CPM',log2.threshold=1){
  #' @param scDATA_1是测试集，用于仿真生成bulk矩阵
  #' @param scDATA_2是训练集，用于产生参考C和求取markerGene
  #' @param leastNum1是测试集QC时保留的细胞类型的最小细胞数目
  #' @param leastNum2是训练集QC时保留的细胞类型的最小细胞数目
  #' @param plotmarker是训练集寻找markerGene时是否画拟合图像
  #' @param Num.mixtures是最终的bulkData的样本数目，也是仿真的抽样次数
  #' @param pool.size是每次仿真抽样的容积
  #' @param norm1对测试集结果标准化
  #' @param norm2对训练集结果标准化
  ComDate=NULL
  # 质量控制
  ComDate <- list(
    scData_1 = scData_1,
    Test = CutCTofscDATA(QCscDATA(scData_1,leastNum1),cutCT=c("unknown")),
    scData_2 = scData_2,
    Train = CutCTofscDATA(QCscDATA(scData_2,leastNum2),cutCT=c("unknown"))
  )
  ComDate$Train <- CutCTofscDATA(ComDate$Train)
  # 模拟真实情况下匹配gene，但是并不匹配celltype（事实上已知BULK数据的情况下，只有gene信息无celltype信息）
  cat( "The number of genes matching TrainSet and TestSet:",
       length(intersect(rownames(ComDate$Test$data), rownames(ComDate$Train$data))))
  cat( "\nCelltypes in TrainSet:"); print(table(ComDate$Train$pData$cellType))
  cat( "\nCelltypes in TestSet:"); print(table(ComDate$Test$pData$cellType))
  
  keepTestGene <- which(rownames(ComDate$Test$data) %in% 
                          intersect(rownames(ComDate$Test$data), rownames(ComDate$Train$data)))
  keepTrainGene <- which(rownames(ComDate$Train$data) %in% 
                           intersect(rownames(ComDate$Test$data), rownames(ComDate$Train$data)))
  
  ComDate$Test$data <- ComDate$Test$data[keepTestGene,]
  # ComDate$Test$pData <- ComDate$Test$pData[,]
  # ComDate$Test$original_cell_names <- ComDate$Test$original_cell_names[]
  ComDate$Train$data <- ComDate$Train$data[keepTrainGene,]
  # ComDate$Train$pData <- ComDate$Train$pData[,]
  # ComDate$Train$original_cell_names <- ComDate$Train$original_cell_names[]
  
  # 寻找marker Gene
  ComDate$markers <- Find_markerGene_limma(ComDate$Train$data, plotmarker,norm2,log2.threshold)
  cat("\nDifferent celltypes's marker Gene number of RAW:");print(sort(table(ComDate$markers$CT)))
  # ComDate$df = which(ComDate$markers$FDR <= 0.1 & ComDate$markers$P.value<0.1)
  # ComDate$markerslist = ComDate$markers[ComDate$df,]
  ComDate$markerslist = ComDate$markers
  cat("\nDifferent celltypes's marker Gene number after shift:");print(sort(table(ComDate$markerslist$CT)))
  
  # GET matrix c(reference from train and real from test)
  ComDate$C_ref <- Get_C(ComDate$Train$data)[['C_ref']]
  ComDate$refProfiles.var <- Get_C(ComDate$Train$data)[['refProfiles.var']]
  ComDate$C <- Get_C(ComDate$Test$data)[['C_ref']]
  
  # simulate T,P from test data
  require(dplyr);
  colnames(ComDate$Test$data) <- ComDate$Test$original_cell_names
  generator <- Generator(sce = ComDate$Test$data, phenoData = ComDate$Test$pData, Num.mixtures, pool.size)
  ComDate$T <- generator[["T"]]
  ComDate$P <- generator[["P"]]
  
  # for T,C_ref,refProfiles.var: shift gene and ct(when markers list ct != markers ct)
  ComDate$indata =NULL
  gene_ <- ComDate$markerslist$gene
  ct_ <- unique(ComDate$markerslist$CT)
  ComDate$indata$T <- ComDate$T[gene_,]
  ComDate$indata$C_ref <- ComDate$C_ref[gene_,ct_]
  ComDate$indata$refProfiles.var <- ComDate$refProfiles.var[gene_,ct_]
  
  ComDate$indata$C <- ComDate$C
  ComDate$indata$P <- ComDate$P
  
  
  # CPM标准化
  library(scater)
  if(norm1 %in% c("cpm","CPM")){
    ComDate$indata$T <- calculateCPM(as.matrix(ComDate$indata$T))
    ComDate$indata$C <- calculateCPM(as.matrix(ComDate$indata$C))
  }
  if(norm2 %in% c("cpm","CPM")){
    ComDate$indata$C_ref <- calculateCPM(as.matrix(ComDate$indata$C_ref))
  }
  ComDate$indata$P <- as.matrix(ComDate$indata$P)
  ComDate$indata$P <- as.matrix(ComDate$indata$P)
  
  return(ComDate)
}


CountAllResults <- function(result,MyMethodName){
  print('there are methods:')
  MethodName <- names(result);MethodName
  MatrixCref <- intersect(c("CDSC3","NNLS", "OLS" ,"FARDEEP", "CIBERSORT",        
                            "deconRNASeq","RLR" ,"DCQ" ,"elastic_net" ,"ridge","lasso"  ,"EPIC"),MyMethodName)
  ScCref <- intersect(c("MuSiC","Bisque","SCDC","DWLS"),MyMethodName)
  NoCref <- intersect(c("CDSC2","DSA" ,"ssKL" , "ssFrobenius","deconf" ,
                        "TOAST" ,"Linseed","CellDistinguisher"),MyMethodName)
  
  result$all <- NULL
  result$all <- result$CDSC3$result[4:9]
  myMatrixCref <- intersect(MatrixCref,MethodName)
  if(length( myMatrixCref ) != 0 ){
    for (i in 2:length(myMatrixCref)) {
    result$all[i,] <- cbind(result[[myMatrixCref[i]]]$result, matrix(0,1,6-length(result[[myMatrixCref[i]]]$result)))
    }}
  myScCref <- intersect(ScCref,MethodName)
  if(length( myScCref ) != 0 ){
    num = nrow(result$all)
    for (i in (num+1):(num+length(myScCref))) {
    result$all[i,] <- cbind(result[[myScCref[i-num]]]$result, matrix(0,1,6-length(result[[myScCref[i-num]]]$result)))
    }}
  myNoCref <- intersect(NoCref,MethodName)
  if(length( myNoCref ) != 0 ){
    num = nrow(result$all)
    for (i in (num+1):(num+length(myNoCref)) ) {
      if(myNoCref[i-num] == "CDSC2"){
        result$all[i,] <- result$CDSC2$result[4:9]
      }else{
        result$all[i,] <- cbind(result[[myNoCref[i-num]]]$result, matrix(0,1,6-length(result[[myNoCref[i-num]]]$result)))
      }
      
    }}
  
  rownames(result$all) <- c(myMatrixCref,myScCref,myNoCref)
  # result$all <- result$CDSC3$result[4:9]
  # result$all[2,] <- cbind(result$nnls$result, matrix(0,1,6-length(result$nnls$result)))
  # result$all[3,] <- cbind(result$ols$result,matrix(0,1,6-length(result$ols$result)))
  # result$all[4,] <- cbind(result$fardeep$result,matrix(0,1,6-length(result$fardeep$result)))
  # result$all[5,] <- cbind(result$cibersort$result,matrix(0,1,6-length(result$cibersort$result)))
  # result$all[6,] <- cbind(result$deconRNASeq$result,matrix(0,1,6-length(result$deconRNASeq$result)))
  # result$all[7,] <- cbind(result$rlr$result,matrix(0,1,6-length(result$rlr$result)))
  # result$all[8,] <- cbind(result$DCQ$result,matrix(0,1,6-length(result$DCQ$result)))
  # result$all[9,] <- cbind(result$elastic_net$result,matrix(0,1,6-length(result$elastic_net$result)))
  # result$all[10,] <- cbind(result$ridge$result,matrix(0,1,6-length(result$ridge$result)))
  # result$all[11,] <- cbind(result$lasso$result,matrix(0,1,6-length(result$lasso$result)))
  # result$all[12,] <- cbind(result$EPIC$result,matrix(0,1,6-length(result$EPIC$result)))
  # 
  # result$all[13,] <- cbind(result$MuSiC$result,matrix(0,1,6-length(result$MuSiC$result)))
  # result$all[14,] <- cbind(result$BisqueRNA$result,matrix(0,1,6-length(result$BisqueRNA$result)))
  # result$all[15,] <- cbind(result$deconvSeq$result,matrix(0,1,6-length(result$deconvSeq$result)))
  # result$all[16,] <- cbind(result$SCDC$result,matrix(0,1,6-length(result$SCDC$result)))
  # result$all[17,] <- cbind(result$DWLS$result,matrix(0,1,6-length(result$DWLS$result)))
  #   
  # result$all[18,] <- result$CDSC2$result[4:9]
  # result$all[19,] <- cbind(result$dsa$result,matrix(0,1,6-length(result$dsa$result)))
  # result$all[20,] <- cbind(result$sskl$result,matrix(0,1,6-length(result$sskl$result)))
  # result$all[21,] <- cbind(result$ssFrobenius$result,matrix(0,1,6-length(result$ssFrobenius$result)))
  # result$all[22,] <- cbind(result$deconf$result,matrix(0,1,6-length(result$deconf$result)))
  # 
  # result$all[23,] <- cbind(result$CDSeq$result,matrix(0,1,6-length(result$CDSeq$result)))
  # result$all[24,] <- cbind(result$TOAST$result,matrix(0,1,6-length(result$TOAST$result)))
  # result$all[25,] <- cbind(result$Linseed$result,matrix(0,1,6-length(result$Linseed$result)))
  # result$all[26,] <- cbind(result$cd$result,matrix(0,1,6-length(result$cd$result)))
  # 
  # 
  # rownames(result$all) <- c("CDSC3",
  #                           "nnls",
  #                           "ols",
  #                           "fardeep",
  #                           "cibersort",
  #                           "deconRNASeq",
  #                           "rlr",
  #                           "DCQ",
  #                           "elastic_net",
  #                           "ridge",
  #                           "lasso",
  #                           "EPIC",
  #                           
  #                           "MuSiC",
  #                           "BisqueRNA",
  #                           "deconvSeq",
  #                           "SCDC",
  #                           "DWLS",
  #                           
  #                           "CDSC2",
  #                           "dsa",
  #                           "sskl",
  #                           "ssFrobenius",
  #                           "deconf",
  #                           "CDSeq",
  #                           "TOAST",
  #                           "Linseed",
  #                           "CellDistinguisher")
  return(result)
}

Deal7DataSet <- function(bulkData){
  bulkData$full_phenoData
  bulkData$full_phenoData$sample <- 1
  colnames(bulkData$full_phenoData) <- c("cellID",	"cellType",	"sampleID")
  bulkData$full_phenoData <- as.data.frame(bulkData$full_phenoData)
  bulkData$data <- as.matrix(bulkData$data)
  bulkData$full_phenoData
  
  bulkData$simulate <- scSimulateC(bulkData, leastNum=0, plotmarker = F)
  bulkData$shifGene <- bulkData$simulate$markerslist$gene
  length(bulkData$shifGene )
  length(which(rownames(bulkData$T) %in% bulkData$shifGene ))
  bulkData$T <- bulkData$T[which(rownames(bulkData$T) %in% bulkData$shifGene ),]
  bulkData$C <- bulkData$C[which(rownames(bulkData$C) %in% bulkData$shifGene ),]
  getPearsonRMSE(as.matrix(bulkData$C), as.matrix(bulkData$simulate$C))
  
  bulkData$T <- as.matrix(bulkData$T)
  bulkData$C <- as.matrix(bulkData$C)
  bulkData$P <- as.matrix(bulkData$P)
  bulkData$C_ref <- bulkData$C
  return(bulkData)
}

squash_axis <- function(from, to, factor) { 
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  library(scales)
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}

