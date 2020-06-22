rm(list = ls())
library(RColorBrewer)
library(scales)
library(rgl)
library(NbClust)
library(pheatmap)
#H-cluster
#load data
data <- as.data.frame(read.csv("gene_count_matrix.csv",row.names="gene_id"))
#转置数据并标准化
dataN <- scale(t(data))
#计算样本间的距离
datad <- dist(dataN)

#进行聚类分析并画图
d <- hclust(datad,method = "complete",members = NULL)
plot(d,labels = NULL,hang = -1)
datadm <- as.matrix(datad)
colData <- read.csv("coldata.txt", sep="\t", row.names=1)
pheatmap(datadm,clustering_distance_rows=datad,
         clustering_distance_cols=datad,
         show_rownames=FALSE,
         cluster_cols=TRUE, 
         border_color = "white",
         annotation_col=colData,scale ="row")
#尝试不同的聚类方法
opar <- par(mfrow = c(2,2))
d1 <- hclust(datad,method = "ward.D",members = NULL)
plot(d1,labels = NULL,hang = -1)
d2 <- hclust(datad,method = "mcquitty",members = NULL)
plot(d2,labels = NULL,hang = -1)
d3 <- hclust(datad,method = "single",members = NULL)
plot(d3,labels = NULL,hang = -1)
d4 <- hclust(datad,method = "average",members = NULL)
plot(d4,labels = NULL,hang= -1)
par(opar)

#判断类数目
rect.hclust(d,k=3)

#k-means聚类
#过滤数据
data$mean <- apply(data,1,mean)
jud <- data[,13] > 10
data <- data[jud,]
dataf <- data[,-13]

dim <- dim(dataf)
sample.row <- sample(seq(1,dim[1]),2000,replace = F)
dataf <- as.data.frame(dataf[sample.row,])
dataf$WT <- apply(dataf[1:4],1,mean)
dataf$C1 <- apply(dataf[5:8],1,mean)
dataf$C2 <- apply(dataf[9:12],1,mean)
dataA <- dataf[,13:15]

#确定最优的K

dataA <- as.data.frame(scale((dataA)))
bk <- NbClust(dataA,min.nc = 2,max.nc=10,method = "ward.D")

#no nstart
wss <- (nrow(dataA)-1)*sum(apply(dataA,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(dataA,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
#have nstart
wss <- (nrow(dataA)-1)*sum(apply(dataA,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(dataA,
                                     centers=i,nstart = 20)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#k-means 聚类

p <- kmeans(dataA,3,nstart = 25,iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(dataA, col=p$clust, pch=16)
plot3d(dataA, col=p$clust)
dataC <- as.data.frame(p$centers)

#data for mev
library(tidyverse)
data$mean <- apply(data,1,mean)
jud <- data[,13] > 10
data <- data[jud,]
dataf <- data[,-13]
dim <- dim(dataf)
sample.row <- sample(seq(1,dim[1]),2000,replace = F)
dataf <- as.data.frame(dataf[sample.row,])

dataf$WT <- apply(dataf[1:4],1,mean)
dataf$C1 <- apply(dataf[5:8],1,mean)
dataf$C2 <- apply(dataf[9:12],1,mean)
dataA <- dataf[,13:15]


dataA <- scale(t(dataA))
dataA <- as.data.frame(scale(t(dataA)))
colnames(dataA)[0] <- c("gene_id")
write.table(dataA,"Mev.test.txt")
