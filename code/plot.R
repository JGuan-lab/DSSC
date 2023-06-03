
#pheatMap of C's ct
library(pheatmap)
STRING_name = ''
map = GetCorMatrix(result$CDSC3$dec$c,scData$Indata$C,matrix = "c");map
plot_pheatmp <- pheatmap(map,cluster_row = FALSE,cluster_cols = FALSE
                         ,cellwidth = 15, cellheight =15#,gaps_row = c(12, 17)
                         ,fontsize = 6
                         # ,breaks = c(seq(0.4,1,by = 0.02))
                         ,display_numbers = TRUE, number_format = "%.2f" # 数字
                         # ,border_color = "black",color = colorRampPalette(MyColor)(60)
                         ,fontsize_row = 9,fontsize_col = 9,angle_col = 45
                         ,annotation_legend = T
                         ,main = STRING_name)

eoffice::topptx(plot_pheatmp,filename = 
                  paste("XXX/pictures/class1_ct_NO_",
                        STRING_name,".pptx",sep=""))




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
                filename = "XXX/pictures/class1_NO_RMSE_C.pptx")

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
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_class1;
eoffice::topptx(plot_class1,
                filename = "XXX/pictures/class1_NO_fanhua_p.pptx")

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
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_class1;
eoffice::topptx(plot_class1,
                filename = "XXX/pictures/class1_NO_fanhua_C.pptx")

# sample bar
rm (list=ls ())
setwd("XXX/1scSimulate")

source("XXX/CDSC.R")
source("XXX/CDSC_expand.R")

# myOneref <- intersect(Oneref,MethodName);myOneref=setdiff(myOneref,c("CDSC2"))
# myTworef <- intersect(NoCref,MethodName);myTworef=union(c("CDSC3"),myTworef);myTworef=setdiff(myTworef,c("CDSC2"))
a = gsub('.rds','',list.files("XXX/1scSimulate/dataSimulate"));a
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
  geom_bar(data=map_mean,mapping=aes(x=factor(group,levels = a),y=mean,fill=group),
           position="dodge", 
           stat="identity", 
           width = 0.7,
           show.legend = F)+  
  scale_fill_manual(values = c("#80AFBF","#80AFBF","#80AFBF","#80AFBF","#80AFBF"))+ #  "#DA635D","#B1938B"

  geom_errorbar(data=map_mean,mapping=aes(x = group,ymin = mean-sd, ymax = mean+sd),
                width = 0.1, 
                color = 'black',
                size=0.8)+ 
  scale_y_continuous(limits =c(0, 1.1) ,expand = c(0,0))+ 
  theme_classic(  
    base_line_size = 1 
  )+
  
  labs(title="",x="Dataset",y="Pearson of samples in P")+ 
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
        legend.title = element_text(color="black", 
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", 
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  
                                   # family = "myFont",
                                   color = "black", 
                                   face = "bold", 
                                   vjust = 1, 
                                   hjust = 1, 
                                   angle = 45), 
        axis.text.y = element_text(size = 13, 
                                   # family = "myFont", 
                                   color = "black", 
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  ) 
# emf(file = "SOD.emf") 真
print(p1) 
eoffice::topptx(p1,
                filename = "XXX/pictures/class1_bar_NO_samples.pptx")

# ---- sample scatter plot----
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
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5))
plot_class2
eoffice::topptx(plot_class2,filename = "XXX/pictures/class1_boxplot_samples.pptx")



#------------WholeBlood-----
rm (list=ls ())
setwd("3realBulk")
source("CDSC.R")
source("CDSC_help.R")
library(ggplot2)
STRING_name <- c("WholeBlood")
WH <- list()

WH$result_lm22 <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_lm22.rds",sep=""))
WH$result_3pbmcs <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_3pbmcs.rds",sep=""))
WH$result_5pbmcs <- readRDS(paste(getwd(),"/bulkResult/result_",STRING_name,"_5pbmcs.rds",sep=""))

RS <- list()
RS$result_lm22 <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_lm22.rds",sep=""))
RS$result_3pbmcs <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_3pbmcs.rds",sep=""))
RS$result_5pbmcs <- readRDS(paste(getwd(),"/bulkData/",STRING_name,"_5pbmcs.rds",sep=""))

#-----------PLOT-----
Refnames <- names(WH);Refnames
# [1] "result_nsclc_sc" "result_lm22"     "result_3pbmcs"   "result_5pbmcs"   "result_FL" 
Refnames <- c("result_lm22","result_3pbmcs","result_5pbmcs")
names(RS)
# names(Aresult)
#-----first----------
OverAllResult <- NULL
OverAllResult <- WH[[Refnames[1]]]$all$Peason_to_P
for (i in 1:length(names(RS))) {
  OverAllResult <- cbind(OverAllResult, WH[[Refnames[i]]]$all$Peason_to_P)
}
OverAllResult <- OverAllResult[,-1]
rownames(OverAllResult) <- rownames(WH[[Refnames[1]]]$all); colnames(OverAllResult) <- Refnames
OverAllResult

OverAllResultC <- NULL
OverAllResultC <- WH[[Refnames[1]]]$all$Peason_to_C
for (i in 1:length(names(RS))) {
  OverAllResultC <- cbind(OverAllResultC, WH[[Refnames[i]]]$all$Peason_to_C)
}
OverAllResultC <- OverAllResultC[,-1]
rownames(OverAllResultC) <- rownames(WH[[Refnames[1]]]$all);colnames(OverAllResultC) <- Refnames
OverAllResultC
### boxplot
MyMethodName <- c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT" ,     
                  "deconRNASeq","RLR","DCQ","elastic_net","ridge", "lasso" ,       
                  "CDSC2","DSA","ssKL","ssFrobenius","deconf","TOAST","Linseed","CellDistinguisher")
MethodName <- names(WH[[Refnames[1]]]);MethodName=setdiff(MethodName, c("all","CellDist.deconv"));MethodName
MatrixCref <- intersect(c("CDSC3","NNLS","OLS","FARDEEP","CIBERSORT","deconRNASeq","RLR","DCQ","elastic_net","ridge","lasso","EPIC"),MyMethodName)
ScCref <- intersect(c("MuSiC","Bisque","deconvSeq","SCDC","DWLS"),MyMethodName)
NoCref <- intersect(c("CDSC2","DSA","ssKL","ssFrobenius","deconf","CDSeq","TOAST","Linseed","CellDistinguisher"),MyMethodName)

Oneref <- union(MatrixCref,NoCref)

myOneref <- intersect(Oneref,MethodName);#myOneref=setdiff(myOneref,c("CDSC2"))
myTworef <- intersect(NoCref,MethodName);myTworef=union(c("CDSC3"),myTworef);#myTworef=setdiff(myTworef,c("CDSC2"))
myMatrixCref <- intersect(MatrixCref,MethodName)
myScCref <- intersect(ScCref,MethodName)
myNoCref <- intersect(NoCref,MethodName)
# Assign the correct labels to the CDSC3 and CDSC2 methods first, but the order of labels varies for different data.

Refnames
for (i in Refnames) {
  ctlabels <- Row_label(WH[[i]]$CDSC3$dec$c,RS[[i]]$indata$C_ref,leastnum=3)
  rownames(WH[[i]]$CDSC3$dec$p) <- ctlabels
  colnames(WH[[i]]$CDSC3$dec$c) <- ctlabels
  
  if(i == "result_lm22"){
    WH[[i]]$CDSC3$p <- MergeCellType(WH[[i]]$CDSC3$dec$p,'p')
    WH[[i]]$CDSC3$c <- MergeCellType(WH[[i]]$CDSC3$dec$c,'c')
  }else{
    WH[[i]]$CDSC3$p <- WH[[i]]$CDSC3$dec$p
    WH[[i]]$CDSC3$c <- WH[[i]]$CDSC3$dec$c
  }
  
  ctlabels_ <- c("B.cells", "T.cells.CD8", "T.cells.CD4", "NK.cells", "Monocytes")
  WH[[i]]$CDSC3$p <- WH[[i]]$CDSC3$p[ctlabels_, ]
  WH[[i]]$CDSC3$c <- WH[[i]]$CDSC3$c[ ,ctlabels_]
}
for (i in Refnames) {
  ctlabels <- Row_label(WH[[i]]$CDSC2$dec$c,RS[[i]]$indata$C_ref,leastnum=3)
  rownames(WH[[i]]$CDSC2$dec$p) <- ctlabels
  colnames(WH[[i]]$CDSC2$dec$c) <- ctlabels
  
  if(i == "result_lm22"){
    WH[[i]]$CDSC2$p <- MergeCellType(WH[[i]]$CDSC2$dec$p,'p')
    WH[[i]]$CDSC2$c <- MergeCellType(WH[[i]]$CDSC2$dec$c,'c')
  }else{
    WH[[i]]$CDSC2$p <- WH[[i]]$CDSC2$dec$p
    WH[[i]]$CDSC2$c <- WH[[i]]$CDSC2$dec$c
  }
  
  ctlabels_ <- c("B.cells", "T.cells.CD8", "T.cells.CD4", "NK.cells", "Monocytes")
  WH[[i]]$CDSC2$p <- WH[[i]]$CDSC2$p[ctlabels_, ]
  WH[[i]]$CDSC2$c <- WH[[i]]$CDSC2$c[ ,ctlabels_]
}
# Assign consistent cell type order
Refnames
RS_ok<-list()
for(i_method in MethodName){
  for (i_data in Refnames){
    WH[[i_data]][[i_method]]$p <- WH[[i_data]][[i_method]]$p[ctlabels_, ]
    WH[[i_data]][[i_method]]$c <- WH[[i_data]][[i_method]]$c[ ,ctlabels_]
    RS[[i_data]]$indata$P <- RS[[i_data]]$indata$P[ctlabels_, ]
    RS[[i_data]]$indata$C <- RS[[i_data]]$indata$C[,ctlabels_]
  }
}
# RS_ok$indata$P <- RS[[Refnames[1]]]$indata$P[ctlabels_, ]

Refnames 
MethodName
library(RColorBrewer)
## partial-P-sample

PartSampleResult <- NULL
plot_1 <- NULL
for (i in Refnames) {
  PartSampleResult[[i]] <- NULL
  PartSampleResult[[i]] <- data.frame(pearson = diag(cor(WH[[i]]$CDSC3$p,RS[[i]]$indata$P)),
                                      method = c("CDSC3"),row.names = NULL)
  nn1 <- nrow(PartSampleResult[[i]])
  for (i_data in i) {
    for (i_method in myOneref) {
      PartSampleResult[[i]] <- rbind(PartSampleResult[[i]], 
                                     data.frame(pearson = diag(cor(WH[[i_data]][[i_method]]$p,RS[[i_data]]$indata$P)),
                                                method = i_method,row.names = NULL))
    }
  }
  PartSampleResult[[i]] <- PartSampleResult[[i]][-(1:nn1),]
  names_1 <-myOneref 
  names_1 [which(names_1 == c("CDSC3"))] <- c("CDSC3")
  names_1 [which(names_1 == c("CDSC2"))] <- c("CDSC2")
  PartSampleResult[[i]][which(PartSampleResult[[i]]$method == c("CDSC3")),2] <- c("CDSC3")
  PartSampleResult[[i]][which(PartSampleResult[[i]]$method == c("CDSC2")),2] <- c("CDSC2")
  plot_1[[i]] = ggplot(PartSampleResult[[i]], aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
    # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
    geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#BEBADA") + # Boxplot 
    geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
    # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
    # sacle_color_brewer(palette='set1')+
    labs(x="Methods",y = "Pearson")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white")) +   #panel.background = element_rect(fill = '#d8dfea')) # 底色
    ggtitle(strsplit(i,"_")[[1]][2]) +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
}
plot_1[["result_lm22"]]
eoffice::topptx(plot_1[["result_lm22"]],filename =
                  "pictures/class4_NO_no1-1.pptx")
plot_1[["result_3pbmcs"]]
eoffice::topptx(plot_1[["result_3pbmcs"]],filename =
                  "pictures/class4_NO_no1-2.pptx")
plot_1[["result_5pbmcs"]]
eoffice::topptx(plot_1[["result_5pbmcs"]],filename =
                  "pictures/class4_NO_no1-3.pptx")

#------plot all----------
PartSampleResultALL <- PartSampleResult[["result_lm22"]]
PartSampleResultALL <- rbind(PartSampleResultALL,PartSampleResult[["result_3pbmcs"]])
PartSampleResultALL <- rbind(PartSampleResultALL,PartSampleResult[["result_5pbmcs"]])

plot_1_all = ggplot(PartSampleResultALL, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white")) +   #panel.background = element_rect(fill = '#d8dfea')) # 底色
  ggtitle("All") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_1_all
eoffice::topptx(plot_1_all,filename = 
                  "pictures/class4_no1-all_samples.pptx")


colorRampPalette(brewer.pal(9,'Reds')[c(1,2,4,6,9)])(5)
facebook = c("#3b5998","#6d84b4", "#afbdd4", "#d8dfea")
google = c("#5380E4", "#E12A3C", "#FFBF03","#00B723")
etsy = c("#F14000", "#67B6C3", "#F0DA47", "#EBEBE6", "#D0D0CB")
twitter = c("#55ACEE", "#292f33", "#8899a6", "#e1e8ed")

## complere-C-celltype
CompCTResultC <- NULL
plot_2 <- NULL
for (i in Refnames) {
  CompCTResultC[[i]] <- NULL
  CompCTResultC[[i]] <- data.frame(pearson = diag(cor(WH[[i]]$CDSC3$c,RS[[i]]$indata$C)),
                                   method = c("CDSC3"),row.names = NULL)
  nn1 <- nrow(CompCTResultC[[i]])
  for (i_data in i) {
    for (i_method in myTworef) {
      CompCTResultC[[i]] <- rbind(CompCTResultC[[i]], 
                                  data.frame(pearson = diag(cor(WH[[i_data]][[i_method]]$c,RS[[i_data]]$indata$C)),
                                             method = i_method,row.names = NULL))
    }
  }
  CompCTResultC[[i]] <- CompCTResultC[[i]][-(1:nn1),]
  names_2 <- myTworef 
  names_2 [which(names_2 == c("CDSC3"))] <- c("CDSC3")
  names_2 [which(names_2 == c("CDSC2"))] <- c("CDSC2")
  CompCTResultC[[i]][which(CompCTResultC[[i]]$method == c("CDSC3")),2] <- c("CDSC3")
  CompCTResultC[[i]][which(CompCTResultC[[i]]$method == c("CDSC2")),2] <- c("CDSC2")
  
  plot_2[[i]] = ggplot(CompCTResultC[[i]], aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
    # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
    geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#91D1C27F") + # Boxplot 
    geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
    # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
    # scale_colour_manual(values = c("black"))
    labs(x="Methods",y = "Pearson")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white")) +   #panel.background = element_rect(fill = '#d8dfea')) # 底色
    ggtitle(strsplit(i,"_")[[1]][2]) +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
}
plot_2[["result_lm22"]]
eoffice::topptx(plot_2[["result_lm22"]],
                filename = "pictures/class4_NO_no1-4.pptx")
plot_2[["result_3pbmcs"]]
eoffice::topptx(plot_2[["result_3pbmcs"]],
                filename = "pictures/class4_NO_no1-5.pptx")
plot_2[["result_5pbmcs"]]
eoffice::topptx(plot_2[["result_5pbmcs"]],
                filename = "pictures/class4_NO_no1-6.pptx")

CompCTResultCALL <- CompCTResultC[["result_lm22"]]
CompCTResultCALL <- rbind(CompCTResultCALL,CompCTResultC[["result_3pbmcs"]])
CompCTResultCALL <- rbind(CompCTResultCALL,CompCTResultC[["result_5pbmcs"]])

plot_2_all = ggplot(CompCTResultCALL, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#80AFBF") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white")) +   #panel.background = element_rect(fill = '#d8dfea')) # 底色
  ggtitle("All") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_2_all
eoffice::topptx(plot_2_all,filename = "pictures/class4_no1-all_ct.pptx")

#-----The second method - first finding the average, then finding the evaluation index-----
Refnames 
MethodName 

WH_OK <- list()
RS_OK <- list()
# MethodName <- names(WH[[Refnames[1]]]);MethodName=setdiff(MethodName,c("all","lasso","CellDist.deconv"));MethodName
MethodName <- names(WH[[Refnames[1]]]);MethodName=setdiff(MethodName, c("all","CellDist.deconv"));MethodName

for (i_method in MethodName) {
  WH_OK[[i_method]]$p <- SumEqual_1(WH[[Refnames[1]]][[i_method]]$p)
  for (i in 2:length(Refnames)){
    WH_OK[[i_method]]$p <- WH_OK[[i_method]]$p + SumEqual_1(WH[[Refnames[i]]][[i_method]]$p)
  }
  WH_OK[[i_method]]$p = WH_OK[[i_method]]$p/5
  WH_OK[[i_method]]$p <- WH_OK[[i_method]]$p[ctlabels_, ]
}

RS_OK$indata$P <- RS[[Refnames[1]]]$indata$P[ctlabels_, ]
##  Overall indicators after averaging
for(i_method in MethodName){
  WH_OK[[i_method]]$result_p = getPearsonRMSE(WH_OK[[i_method]]$p,RS_OK$indata$P)
}
WH_OK$all_result <- NULL
WH_OK$all_result <- WH_OK$CDSC3$result_p
for (i_method in MethodName) {
  WH_OK$all_result <- rbind(WH_OK$all_result,WH_OK[[i_method]]$result_p)
};WH_OK$all_result <- WH_OK$all_result[-1,]
rownames(WH_OK$all_result) <- MethodName
WH_OK$all_result

##  partial-P-sample
EqualPartSampleResult <- NULL
EqualPartSampleResult <- data.frame(pearson = diag(cor(WH_OK[[MethodName[1]]]$p,RS_OK$indata$P)),
                                    method = MethodName[1],row.names = NULL)
nn1 <- nrow(EqualPartSampleResult)
for (i_method in myOneref) {
  EqualPartSampleResult <- rbind(EqualPartSampleResult, data.frame(pearson = diag(cor(WH_OK[[i_method]]$p,RS_OK$indata$P)),
                                                                   method = i_method,row.names = NULL))
}
EqualPartSampleResult <-EqualPartSampleResult[-(1:nn1),]
names_1 <-myOneref 
names_1 [which(names_1 == c("CDSC3"))] <- c("CDSC3")
names_1 [which(names_1 == c("CDSC2"))] <- c("CDSC2")
EqualPartSampleResult[which(EqualPartSampleResult$method == c("CDSC3")),2] <- c("CDSC3")
EqualPartSampleResult[which(EqualPartSampleResult$method == c("CDSC2")),2] <- c("CDSC2")
plot_3 =  ggplot(EqualPartSampleResult, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#4DBBD57F") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_3
eoffice::topptx(plot_3,filename = "pictures/class4_NOL_no2-1.pptx")


#-----The third type - calculate the indicators first; Reaverage-----
# After calculating the overall similarity indicators of the P matrix for each reference, average the indicators
ThreeOverAllResult <- as.data.frame(round(as.numeric(rowSums(OverAllResult)/ncol(OverAllResult)),8))
rownames(ThreeOverAllResult) <- rownames(OverAllResult)
ThreeOverAllResult
ThreeOverAllResultC <- as.data.frame(round(as.numeric(rowSums(OverAllResultC)/ncol(OverAllResult)),8))
rownames(ThreeOverAllResultC) <- rownames(OverAllResultC)
ThreeOverAllResultC

## 部分反卷积-P-samples
TwoPartSampleResult <- NULL
TwoPartSampleResult$all <- NULL
for (i_method in myOneref) {
  TwoPartSampleResult[[i_method]] <- data.frame(pearson = diag(cor(WH[[Refnames[1]]]$CDSC3$p,RS[[Refnames[1]]]$indata$P)),
                                                row.names = NULL)
  for (i_data in Refnames) {
    TwoPartSampleResult[[i_method]] <- cbind(TwoPartSampleResult[[i_method]],
                                             data.frame(Pearson = diag(cor(WH[[i_data]][[i_method]]$p,RS[[i_data]]$indata$P)),
                                                        row.names = NULL))
  }
  TwoPartSampleResult[[i_method]] <- TwoPartSampleResult[[i_method]][,-1]
  colnames(TwoPartSampleResult[[i_method]]) <- Refnames
  TwoPartSampleResult[[i_method]] <- as.data.frame(round(as.numeric(rowSums(TwoPartSampleResult[[i_method]])/length(Refnames)), 8))
  TwoPartSampleResult$all <- rbind(TwoPartSampleResult$all,
                                   data.frame(pearson = TwoPartSampleResult[[i_method]],
                                              method = i_method,row.names = NULL))
}
colnames(TwoPartSampleResult$all) <- c("pearson","method")
TwoPartSampleResult$all[which(TwoPartSampleResult$all$method == c("CDSC3")),2] <- c("CDSC3")
TwoPartSampleResult$all[which(TwoPartSampleResult$all$method == c("CDSC2")),2] <- c("CDSC2")

plot_4 = ggplot(TwoPartSampleResult$all, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#00A0877F") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_4;
eoffice::topptx(plot_4,filename = "pictures/class4_NO_no3-1.pptx")

## complete-C-celltype
TwoCompCTResultC <- NULL
TwoCompCTResultC$all <- NULL
for (i_method in myTworef) {
  TwoCompCTResultC[[i_method]] <- data.frame(pearson = diag(cor(WH[[Refnames[1]]]$CDSC2$c,RS[[Refnames[1]]]$indata$C)),
                                             row.names = NULL)
  for (i_data in Refnames) {
    TwoCompCTResultC[[i_method]] <- cbind(TwoCompCTResultC[[i_method]],
                                          data.frame(Pearson = diag(cor(WH[[i_data]][[i_method]]$c,RS[[i_data]]$indata$C)),
                                                     row.names = NULL))
  }
  TwoCompCTResultC[[i_method]] <- TwoCompCTResultC[[i_method]][,-1]
  colnames(TwoCompCTResultC[[i_method]]) <- Refnames
  TwoCompCTResultC[[i_method]] <- as.data.frame(round(as.numeric(rowSums(TwoCompCTResultC[[i_method]])/length(Refnames)), 8))
  TwoCompCTResultC$all <- rbind(TwoCompCTResultC$all,
                                data.frame(pearson = TwoCompCTResultC[[i_method]],
                                           method = i_method,row.names = NULL))
}
colnames(TwoCompCTResultC$all) <- c("pearson","method")
TwoCompCTResultC$all[which(TwoCompCTResultC$all$method == c("CDSC3")),2] <- c("CDSC3")
TwoCompCTResultC$all[which(TwoCompCTResultC$all$method == c("CDSC2")),2] <- c("CDSC2")

plot_5 = ggplot(TwoCompCTResultC$all, aes(factor(method,levels=names_1),pearson,fill=method)) +  # background
  # geom_boxplot(notch = T,width=0.5,outlier.color = "red",outlier.shape = 2,outlier.size = 3) + # Boxplot 
  geom_boxplot(notch = F,width=0.5,outlier.shape = NA,fill="#00A0877F") + # Boxplot 
  geom_jitter(shape=16, position=position_jitter(0.1),show.legend = F) + #plot Scatter diagram
  # stat_summary(fun = mean, geom = "point", shape = 23, size=4, aes(color=paste("mean","black")),show.legend = F)+ # add Mean
  # scale_colour_manual(values = c("black"))
  labs(x="Methods",y = "Pearson")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),plot.background = element_rect(fill = "white"))    #panel.background = element_rect(fill = '#d8dfea')) # 底色
# ggtitle("I'm a titile") +theme(plot.title = element_text(hjust = 0.5)) #设置标题居中
plot_5
eoffice::topptx(plot_5,filename = "pictures/class4_NO_no3-2.pptx")

