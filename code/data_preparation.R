rm (list=ls ())
set.seed(44)
source("CDSC.R")
source("function_help.R")

##------------Intra-data----
Segerstolpe <- list(data = readRDS("Segerstolpe.rds"), full_phenoData = readRDS("Segerstolpe_phenoData.rds"))
Camp <- list(data = readRDS("Camp.rds"),full_phenoData = readRDS("Camp_phenoData.rds"))
Darmanis <- list(data = readRDS("Darmanis.rds"),full_phenoData = readRDS("Darmanis_phenoData.rds"))
Manno <- list(data = readRDS("Manno.rds"),full_phenoData = readRDS("Manno_phenoData.rds"))
Nestorowa <- list(data = readRDS("Nestorowa.rds"),full_phenoData = readRDS("Nestorowa_phenoData.rds"))
Manno <- list(data = readRDS("Manno.rds"),full_phenoData = readRDS("Manno_phenoData.rds"))

#---Segerstolpe--
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

#---Camp--------
table(Camp$full_phenoData$cellType)
Camp$keep = which(Camp$full_phenoData$cellType != "Unknown")
Camp$full_phenoData <- Camp$full_phenoData[Camp$keep, ]
Camp$data <- Camp$data[,Camp$keep]
table(Camp$full_phenoData$cellType)
Camp$full_phenoData$cellType <- gsub(" ","_",Camp$full_phenoData$cellType)
table(Camp$full_phenoData$cellType)
# # Camp$simualte <- scSimulate(Camp, leastNum=20, plotmarker = F)
type(Camp$data)
Camp$simulate1 <- scSimulateSplit(Camp,
                                  leastNum=20,plotmarker = F,
                                  norm1 = "none",log2.threshold = log2(2))
# Camp$simulate1 <- scSimulateSplit(Camp, leastNum=20, plotmarker = F,norm1 = "cpm")
nrow(Camp$simulate1$markerslist)
table(Camp$simulate1$markerslist$CT);nrow(Camp$simulate1$markerslist)
Camp <- scSimulateShift(Camp,"all",FALSE)
Camp$simulate1$T["RTN1",1:5]
Camp$Indata$T["RTN1",1:5]
saveRDS(Camp,"XXX.rds")

#---Darmanis--------
table(Darmanis$full_phenoData$cellType)

table(Darmanis$full_phenoData$cellType)
# Darmanis$full_phenoData$cellType <- gsub(" ","_",Darmanis$full_phenoData$cellType)
table(Darmanis$full_phenoData$cellType)
# Darmanis$simualte <- scSimulate(Darmanis, leastNum=20, plotmarker = F)
Darmanis$simulate1 <- scSimulateSplit(Darmanis,
                                      leastNum=20, plotmarker = F,
                                      norm1 = "CPM",log2.threshold = 1)
# table(Darmanis$simulate1$markers$CT)
nrow(Darmanis$simulate1$markers)
table(Darmanis$simulate1$markers$CT);nrow(Darmanis$simulate1$markers)
Darmanis <- scSimulateShift(Darmanis,"all",TRUE)
Darmanis$simulate1$T["NEUROD6",1:5]
Darmanis$Indata$T["NEUROD6",1:5]
saveRDS(Darmanis,"XXX.rds")

#---Nestorowa--------
table(Nestorowa$full_phenoData$cellType)
Nestorowa$keep = which(Nestorowa$full_phenoData$cellType != "Unknown")
Nestorowa$full_phenoData <- Nestorowa$full_phenoData[Nestorowa$keep, ]
Nestorowa$data <- Nestorowa$data[,Nestorowa$keep]
table(Nestorowa$full_phenoData$cellType)
# Nestorowa$full_phenoData$cellType <- gsub(" ","_",Nestorowa$full_phenoData$cellType)
table(Nestorowa$full_phenoData$cellType)
# Nestorowa$simualte <- scSimulate(Nestorowa, leastNum=20, plotmarker = F)
Nestorowa$data[1:10,1:10]
Nestorowa$simulate1 <- scSimulateSplit(Nestorowa,
                                       leastNum=20, plotmarker = F,
                                       norm1 = "none",log2.threshold = log2(1),
)
nrow(Nestorowa$simulate1$markerslist)
table(Nestorowa$simulate1$markerslist$CT);nrow(Nestorowa$simulate1$markerslist)
Nestorowa <- scSimulateShift(Nestorowa,"all",FALSE)# FALSE \C\C_REF
Nestorowa$simulate1$T["ENSMUSG00000060131",1:5]
Nestorowa$Indata$T["ENSMUSG00000060131",1:5]
# 
saveRDS(Nestorowa,"XXX.rds")

#---Manno--------
which(rownames(Manno$data) == "SEPT3")#[1] 15670
which(rownames(Manno$data) == "'SEPT3'")#[1] 21
all(Manno$data[which(rownames(Manno$data) == "'SEPT3'"),] ==
      Manno$data[which(rownames(Manno$data) == "SEPT3"),])#[1] FALSE
names <- rownames(Manno$data);names[21] = c("SEPT3_2")
rownames(Manno$data) = names

table(Manno$full_phenoData$cellType)
Manno$keep = which(Manno$full_phenoData$cellType != "Unk")
Manno$full_phenoData <- Manno$full_phenoData[Manno$keep, ]
Manno$data <- Manno$data[,Manno$keep]
sort(table(Manno$full_phenoData$cellType))

Manno$simulate1 <- scSimulateSplit(Manno,
                                   leastNum=60, plotmarker = F,
                                   norm1 = "CPM",log2.threshold = log2(2))
nrow(Manno$simulate1$markers)
table(Manno$simulate1$markers$CT);nrow(Manno$simulate1$markers);
length(unique(Manno$simulate1$markers$CT))
Manno <- scSimulateShift(Manno,"all",TRUE)
Manno$simulate1$T["ONECUT2",1:5]
Manno$Indata$T["ONECUT2",1:5]

saveRDS(Manno,"XXX.rds")

##-----------Inter-data----
#--human pancreas--
Baron <- list(data = readRDS("Baron.rds"),full_phenoData = readRDS("Baron_phenoData.rds"))
Muraro <- list(data = readRDS("Muraro.rds"),full_phenoData = readRDS("Muraro_phenoData.rds"))
Segerstolpe <- list(data = readRDS("Segerstolpe.rds"),full_phenoData = readRDS("Segerstolpe_phenoData.rds"))
#--Mouse Retina--
Macosko <- list(data = readRDS("Macosko.rds"),full_phenoData = readRDS("Macosko_phenoData.rds"))
Shekhar <- list(data = readRDS("shekhar.rds"),full_phenoData = readRDS("shekhar_phenoData.rds"))

a = c("Baron_Muraro","Muraro_Baron","Baron_Segerstolpe","Segerstolpe_Baron","Segerstolpe_Muraro",
      "Muraro_Segerstolpe","Macosko_Shekhar","Shekhar_Macosko" );
#easy
ComDate1 <-NULL
ComDate1 <- Combine_2scdata_dream(Muraro,Segerstolpe,norm1 =  'none',
                                  norm2 ='CPM',log2.threshold = log(2))
saveRDS(ComDate1,"Muraro_Segerstolpe.rds")

#hard
ComDate2 <-NULL
ComDate2 <- Combine_2scdata_hard( Baron,Muraro,norm1 = 'cpm',
                                  norm2 =  'none',log2.threshold = log(2))
saveRDS(ComDate2,"Baron_Muraro.rds")

##------------mixture data----
BreastBlood <- list(
  T = read.table("BreastBlood_GSE29832/mix.txt",row.names = 1,header = T),
  C = read.table("BreastBlood_GSE29832/sig.txt",row.names = 1,header = T,sep = "\t"),
  P = read.table("BreastBlood_GSE29832/coef.txt",row.names = 1,header = T,sep = "\t"),
  data = read.table("BreastBlood_GSE29832/pure.txt",row.names = 1,header = T,sep = "\t"),
  full_phenoData = read.table("BreastBlood_GSE29832/pure_annotations.txt",header = T,sep = "\t"),
  ShiftGene = read.table("BreastBlood_GSE29832/ShiftGene.txt",sep = "\t"),
  CT = read.table("BreastBlood_GSE29832/CT.txt",sep = "\t")
  )
# GSE29832 <- getGEO("GSE29832")
#
CellLines <- list(
  T = read.table("CellLines_GSE11058/mix.txt"),
  C = read.table("CellLines_GSE11058/sig.txt",row.names = 1,header = T,sep = "\t"),
  P = read.table("CellLines_GSE11058/coef.txt",row.names = 1,header = T,sep = "\t"),
  data = read.table("CellLines_GSE11058/pure.txt",row.names = 1,header = T,sep = "\t"),
  full_phenoData = read.table("CellLines_GSE11058/pure_annotations.txt",header = T,sep = "\t"),
  ShiftGene = read.table("CellLines_GSE11058/ShiftGene.txt",sep = "\t"),
  CT = read.table("CellLines_GSE11058/CT.txt",sep = "\t"),

)
# GSE11058 <- getGEO("GSE11058")

LiverBrainLung <- list(
  T = read.table("LiverBrainLung_GSE19830/mix.txt",row.names = 1,header = T),
  C = read.table("LiverBrainLung_GSE19830/sig.txt",row.names = 1,header = T,sep = "\t"),
  P = read.table("LiverBrainLung_GSE19830/coef.txt",row.names = 1,header = T,sep = "\t"),
  data = read.table("LiverBrainLung_GSE19830/pure.txt",row.names = 1,header = T,sep = "\t"),
  full_phenoData = read.table("LiverBrainLung_GSE19830/pure_annotations.txt",header = T,sep = "\t"),
  ShiftGene = read.table("LiverBrainLung_GSE19830/ShiftGene.txt",sep = "\t"),
  CT = read.table("LiverBrainLung_GSE19830/CT.txt",sep = "\t")

)
# GSE19830<- getGEO("GSE19830")

RatBrain<- list(
  T = read.table("RatBrain_GSE19380/mix.txt",row.names = 1,header = T),
  C = read.table("RatBrain_GSE19380/sig.txt",row.names = 1,header = T,sep = "\t"),
  P = read.table("RatBrain_GSE19380/coef.txt",row.names = 1,header = T,sep = "\t"),
  data = read.table("RatBrain_GSE19380/pure.txt",row.names = 1,header = T,sep = "\t"),
  full_phenoData = read.table("RatBrain_GSE19380/pure_annotations.txt",header = T,sep = "\t"),
  ShiftGene = read.table("RatBrain_GSE19380/ShiftGene.txt",sep = "\t"),
  CT = read.table("RatBrain_GSE19380/CT.txt",sep = "\t")

)
# GSE19380 <- getGEO("GSE19380")

Retina <- list(
  T = read.table("Retina_GSE33076/mix.txt",row.names = 1,header = T),
  C = read.table("Retina_GSE33076/sig.txt",row.names = 1,header = T,sep = "\t"),
  P = read.table("Retina_GSE33076/coef.txt",row.names = 1,header = T,sep = "\t"),
  data = read.table("Retina_GSE33076/pure.txt",row.names = 1,header = T,sep = "\t"),
  full_phenoData = read.table("Retina_GSE33076/pure_annotations.txt",header = T,sep = "\t"),
  ShiftGene = read.table("Retina_GSE33076/ShiftGene.txt",sep = "\t"),
  CT = read.table("Retina_GSE33076/CT.txt",sep = "\t")
)
# GSE33076 <- getGEO("GSE33076")
# ----------------BreastBlood-----------------
# 
colnames(BreastBlood$data) <- BreastBlood$full_phenoData$Class
BreastBlood$markers <- Find_markerGene_limma(BreastBlood$data)
BreastBlood$ShiftGene <- BreastBlood$markers$gene

all(rownames(BreastBlood$T) == rownames(BreastBlood$C) )
all(rownames(BreastBlood$T) == rownames(BreastBlood$data))

BreastBlood$ShiftGene <- as.matrix(BreastBlood$ShiftGene)
BreastBlood$CT <- as.matrix(BreastBlood$CT)

length(which(rownames(BreastBlood$T) %in% BreastBlood$ShiftGene))

BreastBlood$indata$T = as.matrix(BreastBlood$T[which(rownames(BreastBlood$T) %in% BreastBlood$ShiftGene),])
BreastBlood$indata$C = as.matrix(BreastBlood$C[which(rownames(BreastBlood$C) %in% BreastBlood$ShiftGene),])
BreastBlood$indata$C_ref = BreastBlood$indata$C

BreastBlood$indata$P = as.matrix(BreastBlood$P)
BreastBlood$indata$gene = BreastBlood$ShiftGene

# BreastBlood <- Deal7DataSet(BreastBlood)
saveRDS(BreastBlood,"BreastBlood.rds")

#----------------CellLines-----------------
#
colnames(CellLines$data) <- CellLines$full_phenoData$Class
CellLines$markers <- Find_markerGene_limma(CellLines$data)
CellLines$ShiftGene <- CellLines$markers$gene

all(rownames(CellLines$T) == rownames(CellLines$C) )
all(rownames(CellLines$T) == rownames(CellLines$data))

CellLines$ShiftGene <- as.matrix(CellLines$ShiftGene)
length(which(rownames(CellLines$T) %in% CellLines$ShiftGene))
CellLines$indata$T = as.matrix(CellLines$T[which(rownames(CellLines$T) %in% CellLines$ShiftGene),])
CellLines$indata$C = as.matrix(CellLines$C[which(rownames(CellLines$C) %in% CellLines$ShiftGene),])
CellLines$indata$C_ref = CellLines$indata$C
CellLines$indata$P = as.matrix(CellLines$P)
CellLines$indata$gene = CellLines$ShiftGene

CellLines <- Deal7DataSet(CellLines)
saveRDS(bulkData,"CellLines.rds")

#----------------LiverBrainLung-----------------
#
all(rownames(LiverBrainLung$T) == rownames(LiverBrainLung$C) )
all(rownames(LiverBrainLung$T) == rownames(LiverBrainLung$data))
LiverBrainLung$ShiftGene <- as.matrix(LiverBrainLung$ShiftGene)
length(which(rownames(LiverBrainLung$T) %in% LiverBrainLung$ShiftGene))
LiverBrainLung$indata$T = as.matrix(LiverBrainLung$T[which(rownames(LiverBrainLung$T) %in% LiverBrainLung$ShiftGene),])
LiverBrainLung$indata$C = as.matrix(LiverBrainLung$C[which(rownames(LiverBrainLung$C) %in% LiverBrainLung$ShiftGene),])
LiverBrainLung$indata$C_ref = LiverBrainLung$indata$C
LiverBrainLung$indata$P = as.matrix(LiverBrainLung$P)
LiverBrainLung$indata$gene = LiverBrainLung$ShiftGene

# LiverBrainLung <- Deal7DataSet(LiverBrainLung)
saveRDS(LiverBrainLung,"LiverBrainLung.rds")

#----------------RatBrain-----------------

all(rownames(RatBrain$T) == rownames(RatBrain$C) )
all(rownames(RatBrain$T) == rownames(RatBrain$data))
RatBrain$ShiftGene <- as.matrix(RatBrain$ShiftGene)
length(which(rownames(RatBrain$T) %in% RatBrain$ShiftGene))
RatBrain$indata$T = as.matrix(RatBrain$T[which(rownames(RatBrain$T) %in% RatBrain$ShiftGene),])
RatBrain$indata$C = as.matrix(RatBrain$C[which(rownames(RatBrain$C) %in% RatBrain$ShiftGene),])
RatBrain$indata$C_ref = RatBrain$indata$C
RatBrain$indata$P = as.matrix(RatBrain$P)
RatBrain$indata$gene = RatBrain$ShiftGene

# RatBrain <- Deal7DataSet(RatBrain)
saveRDS(RatBrain,"RatBrain.rds")

#----------------Retina-----------------
#
all(rownames(Retina$T) == rownames(Retina$C) )
all(rownames(Retina$T) == rownames(Retina$data))
Retina$ShiftGene <- as.matrix(Retina$ShiftGene)
length(which(rownames(Retina$T) %in% Retina$ShiftGene))
Retina$indata$T = as.matrix(Retina$T[which(rownames(Retina$T) %in% Retina$ShiftGene),])
Retina$indata$C = as.matrix(Retina$C[which(rownames(Retina$C) %in% Retina$ShiftGene),])
Retina$indata$C_ref = Retina$indata$C
Retina$indata$P = as.matrix(Retina$P)
Retina$indata$gene = Retina$ShiftGene

# Retina <- Deal7DataSet(Retina)
saveRDS(Retina,"Retina.rds")

##------------whole blood----
WholeBlood <- list(
  T = read.table("WholeBlood_bulk.txt",row.names = 1,header = T),
  P = read.table("WholeBlood_groundtruth.txt",row.names = 1,header = T,sep = "\t")
)
# --LM22-----------
WholeBlood$C$LM22 <- read.table("/LM22.txt",row.names = 1,header = T,sep = "\t")
library(readxl)
WholeBlood$C$LM22_gene_CT <- read_excel("/LM22_geneToCT.xlsx",col_names = T)

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
saveRDS(WholeBlood,"WholeBlood_lm22.rds")

# --3'PBMCs-----------
WholeBlood$C$PBMCs_3 <- read.table("/3PBMCs_scRNAseq_matrix.txt",row.names = 1,header = T,sep = "\t")

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
saveRDS(WholeBlood,"WholeBlood_3pbmcs.rds")

# --5'PBMCs-------------
WholeBlood$C$PBMCs_5 <- read.table("/5PBMCs_scRNAseq_matrix.txt",row.names = 1,header = T,sep = "\t")

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

saveRDS(WholeBlood,"WholeBlood_5pbmcs.rds")
