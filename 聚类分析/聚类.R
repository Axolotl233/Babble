rm(list = ls())

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
jud <- data[,13] > 2500
data <- data[jud,]
dataf <- data[,-13]
dataf$A <- apply(dataf[1:4],1,mean)
dataf$B <- apply(dataf[5:8],1,mean)
dataf$C <- apply(dataf[9:12],1,mean)
dataA <- dataf[,13:15]
dataA <- dataA[1:200,]

#确定最优的K
library(NbClust)
dataA <- scale(t(dataA))
bk <- NbClust(dataA,min.nc = 2,max.nc=10,method = "complete")

#k-means 聚类
p <- kmeans(dataA,2)
plot(dataA,pch=p$cluster)
