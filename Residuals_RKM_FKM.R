##### Check RKM vs FKM
library("fastDummies")
library("philentropy")
library(stats)
library("factoextra")
library("gplots")
library("FactoMineR")
library("ggfortify")
library("pdfCluster")
library("clustrd")
library("tictoc")
library("cluster")
library(doParallel)
library(plyr)



Tandem_res = matrix( nrow =0 , ncol =2)
# Run FKM and RKm on all the simulated data
for (c in c(1,2)){
  for (i in c(16,32,64)){
    for (j in c(3,5,9)){
      for (std in c("low", "high")){
      am_sim=50
      path = paste("C:\\Users\\siets\\Documents\\E&OR\\Thesis\\Proposal\\Simulatie_case",c,i,j,std,sep="_")
      dim_DR=1
      res = RKM_FKM(am_sim, path, dim_DR)}
      Tandem_res = rbind(Tandem_res, res)
}
}
}
Tandem_res


RKM_FKM = function(am_sim, path, dim_DR){
  ##Function used to obtain the complement residuals and the subspace residuals
  RKMcol = data.frame()
  FKMcol = data.frame()
  for (i in c(1:am_sim)){
    num = i
    name = paste(num,"sim.csv", sep="_")
    pa_na = paste(path, name, sep="\\")
    matrix = read.csv(pa_na)
    dt <- as.table(as.matrix(matrix[,1:(dim(matrix)[2]-1)]))
    outRKM = cluspca(dt, 4, dim_DR, method = "RKM", nstart = 25)
    ARIRKM = adj.rand.index(c(outRKM$cluster), c(matrix[,dim(matrix)[2]]))
    
    outFKM = cluspca(dt, 4, dim_DR, method = "FKM",rotation = "varimax",  nstart = 25)
    ARIFKM = adj.rand.index(c(outFKM$cluster), c(matrix[,dim(matrix)[2]]))
    
    #### resid XAA-UFA following FKM
    #compl res = X - XAA
    compl_res = dt - outRKM$obscoord %*%  t(outRKM$attcoord)
    #subsp_resid = XAA' - UFA'
    subsp_resid = outRKM$obscoord %*%  t(outRKM$attcoord) -  to.indicators(outRKM$cluster, exclude.base = FALSE)%*% outRKM$centroid %*%t(outRKM$attcoord)
    obscoors = t(outFKM$attcoord) %*% outFKM$attcoord
    varcomp = var(c(compl_res))
    varsubsp = var(c(subsp_resid))
    RKMcol = rbind(RKMcol, varcomp)
    FKMcol = rbind(FKMcol, varsubsp)
  }
  ARICIL = list("compl_res"=mean(RKMcol), "subsp_res"=mean(FKMcol))
  return(ARICIL)}



