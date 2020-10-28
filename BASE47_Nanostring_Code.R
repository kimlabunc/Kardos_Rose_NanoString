setwd(".")

library(pamr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(gridExtra)
library(xlsx)
library(DESeq2)
library(edgeR)
library(plotrix)
library(RColorBrewer)
library(vegan)
library(VennDiagram)
library(survival)
library(survminer)
library(gplots)

rm(list=ls())
graphics.off()

##########
#Figure 1#
##########

#Was generated without any code or analysis

##########
#Figure 2#
##########

#Classify UNC Training Samples with RNA and Nanostring Data

#Read in the RNASeq BASE47 Prediction Matrix
x <- read.table("BASE47_RNASEQ_PREDICTOR_DO_NOT_EDIT.txt",sep="\t",row.names=1,header=T)
classes <- as.vector(t(x[1,]))
xn <- apply(x[-1,],2,function(x){as.numeric(as.vector(x))})
rownames(xn) <- rownames(x)[-1]


#Read in and classify the RNASeq BASE47 File for the UNC-Training Dataset
test <- read.table("UNC_TMA_RNASeq_BASE47_INPUT.txt", sep="\t",row.names=1,header=T) 
trainSet <- list(x=scale(xn), y=classes, geneid=rownames(xn), genenames=rownames(xn))
mytrain <-pamr.train(trainSet)
pred.class <- pamr.predict(mytrain,scale(test),0)
pred.prob <- pamr.predict(mytrain,scale(test),0,type="posterior")
unc.tma.rna.calls <- as.data.frame(cbind(pred.prob,pred.class))
colnames(unc.tma.rna.calls) <- c("Basal.Correlation","Luminal.Correlation","BASE47.Subtype")
unc.tma.rna.calls$BASE47.Subtype[unc.tma.rna.calls$BASE47.Subtype == 1] <- "Basal"
unc.tma.rna.calls$BASE47.Subtype[unc.tma.rna.calls$BASE47.Subtype == 2] <- "Luminal"

#Read in and classify the Nanostring BASE47 File for the UNC-Training Dataset
test <- read.table("UNC_TMA_Nanostring_BASE47_Null_INPUT.txt", sep="\t",row.names=1,header=T) 
trainSet <- list(x=scale(xn), y=classes, geneid=rownames(xn), genenames=rownames(xn))
mytrain <-pamr.train(trainSet)
pred.class <- pamr.predict(mytrain,scale(test),0)
pred.prob <- pamr.predict(mytrain,scale(test),0,type="posterior")
unc.tma.nano.calls <- as.data.frame(cbind(pred.prob,pred.class))
colnames(unc.tma.nano.calls) <- c("Basal.Correlation","Luminal.Correlation","BASE47.Subtype")
unc.tma.nano.calls$BASE47.Subtype[unc.tma.nano.calls$BASE47.Subtype == 1] <- "Basal"
unc.tma.nano.calls$BASE47.Subtype[unc.tma.nano.calls$BASE47.Subtype == 2] <- "Luminal"

rna.null.model <- merge(x = unc.tma.rna.calls,y = unc.tma.nano.calls,by="row.names")

#This is the table from Figure 2A
table(rna.null.model$BASE47.Subtype.x,rna.null.model$BASE47.Subtype.y)



#Read in the Nanostring BASE47 Prediction Matrix
x <- read.table("BASE47_Nanostring_Predictor_DO_NOT_EDIT.txt",sep="\t",row.names=1,header=T)
classes <- as.vector(t(x[1,]))
xn <- apply(x[-1,],2,function(x){as.numeric(as.vector(x))})
rownames(xn) <- rownames(x)[-1]

#Read in and classify the Nanostring BASE47 File for the UNC-Training Dataset
test <- read.table("UNC_TMA_Nanostring_BASE47_Nano_INPUT.txt", sep="\t",row.names=1,header=T) 
trainSet <- list(x=scale(xn), y=classes, geneid=rownames(xn), genenames=rownames(xn))
mytrain <-pamr.train(trainSet)
pred.class <- pamr.predict(mytrain,scale(test),0)
pred.prob <- pamr.predict(mytrain,scale(test),0,type="posterior")
unc.tma.nano.calls <- as.data.frame(cbind(pred.prob,pred.class))
colnames(unc.tma.nano.calls) <- c("Basal.Correlation","Luminal.Correlation","BASE47.Subtype")
unc.tma.nano.calls$BASE47.Subtype[unc.tma.nano.calls$BASE47.Subtype == 1] <- "Basal"
unc.tma.nano.calls$BASE47.Subtype[unc.tma.nano.calls$BASE47.Subtype == 2] <- "Luminal"

rna.nano.model <- merge(x = unc.tma.rna.calls,y = unc.tma.nano.calls,by="row.names")

#This is the table from Figure 2B, 3 samples were filtered for lack of clinical info
table(rna.nano.model$BASE47.Subtype.x,rna.nano.model$BASE47.Subtype.y)

##########
#Figure 3#
##########

rm(list=ls())

tma.nano <- as.matrix(read.table(file = "TMA_Heatmap_Nanostring_BASE47_Ordered_Med.txt",header = T,sep = "\t",row.names = 1))
tma.pred <- read.delim(file = "TMA_Heatmap_BASE47_Predictions.txt",header = T,sep = "\t")
tma.info <- read.delim(file = "TMA_Heatmap_Clinical_Info.txt",header = T,sep = "\t")

p1 <- ggplot(data=tma.pred, aes(x=reorder(Sample_ID, -Centroid_Correlation), y=Centroid_Correlation, fill=BASE47_Call)) +
  geom_bar(stat="identity") + 
  theme_classic() + 
  xlab("") + 
  ylab("Correlation to Basal Centroid") +
  scale_fill_manual(values=c("red", "blue"))
#Barplot of Figure 3A
p1

#Heatmap of Figure 3A
heatmap.2(tma.nano,
          Colv=FALSE, 
          dendrogram="row",
          col="greenred",
          trace="none",
          labCol=FALSE)


rm(list=ls())

validation.nano <- as.matrix(read.table(file = "Validation_Metadataset_Heatmap_Input.txt",header = T,sep = "\t",row.names = 1))
validation.pred <- read.delim(file = "Validation_Metadataset_Clinical.txt",header = T,sep = "\t")
validation.info <- read.delim(file = "Validation_Metadataset_Heatmap_Colorbar.txt",header = T,sep = "\t")

p1 <- ggplot(data=validation.pred, aes(x=reorder(Sample_ID, -Centroid_Correlation), y=Centroid_Correlation, fill=BASE47_Call)) +
  geom_bar(stat="identity") + 
  theme_classic() + 
  xlab("") + 
  ylab("Correlation to Basal Centroid") +
  scale_fill_manual(values=c("red", "blue"))
#Barplot of Figure 3B
p1

#Heatmap of Figure 3B
heatmap.2(validation.nano,
          Colv=FALSE, 
          dendrogram="row",
          col="greenred",
          trace="none",
          labCol=FALSE)


##########
#Figure 4#
##########

rm(list=ls())

#Read in the RNASeq BASE47 Prediction Matrix
x <- read.table("BASE47_RNASEQ_PREDICTOR_DO_NOT_EDIT.txt",sep="\t",row.names=1,header=T)
classes <- as.vector(t(x[1,]))
xn <- apply(x[-1,],2,function(x){as.numeric(as.vector(x))})
rownames(xn) <- rownames(x)[-1]


#Read in and classify the RNASeq BASE47 File for the UNC-Training Dataset
test <- read.table("UNC_GC_RNASEQ_BASE47_INPUT.txt", sep="\t",row.names=1,header=T) 
trainSet <- list(x=scale(xn), y=classes, geneid=rownames(xn), genenames=rownames(xn))
mytrain <-pamr.train(trainSet)
pred.class <- pamr.predict(mytrain,scale(test),0)
pred.prob <- pamr.predict(mytrain,scale(test),0,type="posterior")
unc.gc.rna.calls <- as.data.frame(cbind(pred.prob,pred.class))
colnames(unc.gc.rna.calls) <- c("Basal.Correlation","Luminal.Correlation","BASE47.Subtype")
unc.gc.rna.calls$BASE47.Subtype[unc.gc.rna.calls$BASE47.Subtype == 1] <- "Basal"
unc.gc.rna.calls$BASE47.Subtype[unc.gc.rna.calls$BASE47.Subtype == 2] <- "Luminal"

#Read in the Nanostring BASE47 Prediction Matrix
x <- read.table("BASE47_Nanostring_Predictor_DO_NOT_EDIT.txt",sep="\t",row.names=1,header=T)
classes <- as.vector(t(x[1,]))
xn <- apply(x[-1,],2,function(x){as.numeric(as.vector(x))})
rownames(xn) <- rownames(x)[-1]

#Read in and classify the Nanostring BASE47 File for the UNC-Training Dataset
test <- read.table("UNC_GC_Nanostring_BASE47_INPUT.txt", sep="\t",row.names=1,header=T) 
trainSet <- list(x=scale(xn), y=classes, geneid=rownames(xn), genenames=rownames(xn))
mytrain <-pamr.train(trainSet)
pred.class <- pamr.predict(mytrain,scale(test),0)
pred.prob <- pamr.predict(mytrain,scale(test),0,type="posterior")
unc.gc.nano.calls <- as.data.frame(cbind(pred.prob,pred.class))
colnames(unc.gc.nano.calls) <- c("Basal.Correlation","Luminal.Correlation","BASE47.Subtype")
unc.gc.nano.calls$BASE47.Subtype[unc.gc.nano.calls$BASE47.Subtype == 1] <- "Basal"
unc.gc.nano.calls$BASE47.Subtype[unc.gc.nano.calls$BASE47.Subtype == 2] <- "Luminal"

rna.nano.model <- merge(x = unc.gc.rna.calls,y = unc.gc.nano.calls,by="row.names")

#This is the table from Figure 4A
table(rna.nano.model$BASE47.Subtype.x,rna.nano.model$BASE47.Subtype.y)
#Figures 4B and 4C are just standard correlation plots

gene.cor <- read.delim(file = "UNC_GC_Gene_Correlations.txt",header = T,sep = "\t")

barplot(gene.cor$Pearson_R)

p1 <- ggplot(gene.cor, aes(x = Pearson_R)) + 
  geom_histogram(bins = 12, color = "black", fill = "gray") +
  theme_classic()
#This is Figure 4D
p1

##########
#SupFig 1#
##########

rm(list=ls())

tech.var <- read.delim(file = "SuppFig1A_Input.txt",header = T,sep = "\t")

p1 <- ggplot(data=tech.var, aes(x=reorder(Name, -WYK3426_CoV), y=WYK3426_CoV, fill=Code.Class)) +
  geom_bar(stat="identity") + 
  theme_classic() + 
  xlab("") + 
  ylab("Coefficient of Variance") +
  scale_fill_manual(values=c("blue", "red"))
#Barplot of Sup Fig 1A
p1

bio.var <- read.delim(file = "SuppFig1B_Input.txt",header = T,sep = "\t")

p2 <- ggplot(data=bio.var, aes(x=reorder(Name, -CoV), y=CoV, fill=CodeClass)) +
  geom_bar(stat="identity") + 
  theme_classic() + 
  xlab("") + 
  ylab("Coefficient of Variance") +
  scale_fill_manual(values=c("blue", "red"))
#Barplot of Sup Fig 1B
p2

##########
#SupFig 2#
##########
rm(list=ls())

temp <- read.delim(file = "SuppFig2A_Input.txt",header = T,sep = "\t")

cor.test(temp$True_Basal_Correlation,temp$DamrauerLog_Basal_Correlation)
p1 <- ggplot(data=temp, aes(x=True_Basal_Correlation, y=DamrauerLog_Basal_Correlation)) +
  geom_point() + 
  theme_classic() + 
  xlab("RNASeq Basal Correlation") + 
  ylab("Damrauer Basal Correlation") + 
  geom_hline(yintercept=0.5) + 
  annotate(geom="text", x=0.9, y=0.1, label="Pearson R=0.254 \n p=0.07",color="black")
p1  


cor.test(temp$True_Basal_Correlation,temp$Nanostring_Basal_Correlation)
p2 <- ggplot(data=temp, aes(x=True_Basal_Correlation, y=Nanostring_Basal_Correlation)) +
  geom_point() + 
  theme_classic() + 
  xlab("RNASeq Basal Correlation") + 
  ylab("Nanostring Basal Correlation") + 
  geom_hline(yintercept=0.5) + 
  annotate(geom="text", x=0.9, y=0.1, label="Pearson R=0.877 \n p < 2.2e-16",color="black")
p2

#Correlation Plots from Sup Fig 2A
ggarrange(p1,p2,ncol = 2,nrow = 1)


temp <- read.delim(file = "SuppFig2B_Input.txt",header = T,sep = "\t")

cor.test(temp$True_Basal_Correlation,temp$Basal)
p3 <- ggplot(data=temp, aes(x=True_Basal_Correlation, y=Basal)) +
  geom_point() + 
  theme_classic() + 
  xlab("True Basal Correlation") + 
  ylab("Simulated Basal Correlation") + 
  geom_hline(yintercept=0.5) + 
  annotate(geom="text", x=0.9, y=0.1, label="Pearson R=0.847 \n p < 2.2e-16",color="black")
#Correlation Plot from Sup Fig 2B
p3

##########
#SupFig 3#
##########
rm(list=ls())

temp <- read.delim(file = "SuppFig3_Input.txt",header = T,sep = "\t")

p1 <- ggplot(data=temp, aes(x=Count, y=temp$Damrauer_Basal)) +
  geom_bar(stat="identity", fill = "red") + 
  theme_classic() + 
  xlab("") + 
  ylab("Basal")
p1

p2 <- ggplot(data=temp, aes(x=Count, y=temp$Damrauer_Luminal)) +
  geom_bar(stat="identity", fill = "blue") + 
  theme_classic() + 
  xlab("") + 
  ylab("Luminal")
p2

p3 <- ggplot(data=temp, aes(x=Count, y=temp$Nanostring_Basal)) +
  geom_bar(stat="identity", fill = "red") + 
  theme_classic() + 
  xlab("") + 
  ylab("Basal")
p3

p4 <- ggplot(data=temp, aes(x=Count, y=temp$Nanostring_Luminal)) +
  geom_bar(stat="identity", fill = "blue") + 
  theme_classic() + 
  xlab("") + 
  ylab("Luminal")
p4

#Centroid Plots from SupFig3A
ggarrange(p1,p2,ncol = 1,nrow = 2)

#Centroid Plots from SupFig3B
ggarrange(p3,p4,ncol = 1,nrow = 2)







