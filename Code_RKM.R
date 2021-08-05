
library("clustrd")
library("FactoMineR")

set.seed("1234")
### Reads student satisfaction file
data = read.csv("C:\\Users\\siets\\Documents\\E&OR\\Thesis\\Rcode\\StudentSF.csv")
matrix = data[,9:22]

res.pca <- PCA(matrix, scale = TRUE)
res.pca$eig[,2]

var= res.pca$sdev^2/sum(res.pca$sdev^2)
### Plots amount of variance explained by each dimension
plot(c(1:14), res.pca$eig[,2], type = "l", xlab = "Added value of dimensions", ylab = "Percentage of variance explained", main="Percentage of variance explained by each dimension")

crit = data.frame()
for (i in c(2:8)){
  outRKM = cluspca(matrix, i, 2, method = "RKM", nstart = 100)
  crit = rbind(crit, outRKM$criterion)
}
### Plots amount of error with different amount of clusters
plot(x=c(2:8),y=t(crit), type="l", xlab = "Amount of clusters", ylab = "Error from loss function", main="Amount of loss for RKM using different clusters")

outRKM = cluspca(matrix, 5, 2, method = "RKM", nstart = 100)
### Plots the cluster allocation using RKM
plot(outRKM$obscoord, pch = outRKM$cluster,  xlab = "Scores for dimension 1", ylab = "Scores for dimension 2",  main = "Cluster allocation for the respondents")

path = "C:\\Users\\siets\\Documents\\E&OR\\Thesis\\Results\\contrib_dim.csv"
write.csv(outRKM$attcoord, file= path)

