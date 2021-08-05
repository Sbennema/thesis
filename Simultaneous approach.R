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



ARI_RKM = function(am_sim, path, dim_DR){
  # Returns average ARI and ASW for RKM in 50 simulations
  ARIcol = data.frame()
  SILcol = data.frame()
  for (i in c(1:am_sim)){
    num = i
    name = paste(num,"sim.csv", sep="_")
    pa_na = paste(path, name, sep="\\")
    matrix = read.csv(pa_na)
    dt <- as.table(as.matrix(matrix[,1:(dim(matrix)[2]-1)]))
    outRKM = cluspca(dt, 4, dim_DR, method = "RKM", nstart = 25)
    sil = silhouette(c(outRKM$cluster), dist(dt))
    meansil = mean(c(sil[,3]))
    ARI = adj.rand.index(c(outRKM$cluster), c(matrix[,dim(matrix)[2]]))
    ARIcol = rbind(ARIcol, ARI)
    SILcol = rbind(SILcol, meansil)
  }
  ARICIL = cbind(ARIcol, SILcol)
  return(ARICIL)
}

ARI_FKM = function(am_sim, path, dim_DR){
  # Returns average ARI and ASW for FKM in 50 simulations
  ARIcol = data.frame()
  SILcol = data.frame()
  for (i in c(1:am_sim)){
    num = i
    name = paste(num,"sim.csv", sep="_")
    pa_na = paste(path, name, sep="\\")
    matrix = read.csv(pa_na)
    dt <- as.table(as.matrix(matrix[,1:(dim(matrix)[2]-1)]))
    outFKM = cluspca(dt, 4, dim_DR, method = "FKM",rotation = "varimax",  nstart = 25)
    sil = silhouette(c(outFKM$cluster), dist(dt))
    meansil = mean(c(sil[,3]))
    ARI = adj.rand.index(c(outFKM$cluster), c(matrix[,dim(matrix)[2]]))
    ARIcol = rbind(ARIcol, ARI)
    SILcol = rbind(SILcol, meansil)
  }
  ARICIL = cbind(ARIcol, SILcol)
  return(ARICIL)
}

ARI_MCA_kM = function(am_sim, path, dim_DR){
  # Returns average ARI and ASW for MCA_KM in 50 simulations
  ARIcol = data.frame()
  SILcol = data.frame()
  
  for (i in c(1:am_sim)){
    num = i
    name = paste(num,"sim.csv", sep="_")
    pa_na = paste(path, name, sep="\\")
    matrix = read.csv(pa_na)
    dt <- as.table(as.matrix(matrix[,1:(dim(matrix)[2]-1)]))
    mat = c(1:(dim(matrix)[1]))
    # Makes an indicator matrix of the numerical matrix
    for (i in c(1:(dim(matrix)[2]-1))){
      column = data.frame(factor(dt[,i]))
      mat = cbind(mat, column)
    }
    out_MCA_kM = clusmca(mat[,2:dim(matrix)[2]], 4, dim_DR, method = "MCAk", nstart = 25, seed = 1234)
    sil = silhouette(c(out_MCA_kM$cluster), dist(dt))
    meansil = mean(c(sil[,3]))
    ARI = adj.rand.index(c(out_MCA_kM$cluster), c(matrix[,dim(matrix)[2]]))
    ARIcol = rbind(ARIcol, ARI)
    SILcol = rbind(SILcol, meansil)
  }
  ARICIL = cbind(ARIcol, SILcol)
  return(ARICIL)
}

ARI_CCA = function(am_sim, path, dim_DR){
  # Returns average ARI and ASW for MCA_KM in 50 simulations
  ARIcol = data.frame()
  SILcol = data.frame()
  for (i in c(1:am_sim)){
    num = i
    name = paste(num,"sim.csv", sep="_")
    pa_na = paste(path, name, sep="\\")
    matrix = read.csv(pa_na)
    dt <- as.table(as.matrix(matrix[,1:(dim(matrix)[2]-1)]))
    doubled = cbind(dt, max(dt)+1-dt)
    CCAZK = CCA(doubled, dim_DR,25)
    sil = silhouette(c(CCAZK$alloc%*%c(1:4)), dist(dt))
    meansil = mean(c(sil[,3]))
    ARI = adj.rand.index(c(CCAZK$alloc%*%c(1:4)), c(matrix[,dim(matrix)[2]]))
    ARIcol = rbind(ARIcol, ARI)
    SILcol = rbind(SILcol, meansil)
  }
  ARICIL = cbind(ARIcol, SILcol)
  return(ARICIL)
}

CCA = function(doubled, dim_DR, ran_starts){
  # Code for own version of CCA
  minloss= Inf
  # Ran starts indicate the amount random initialisations
  for (i in c(1:ran_starts)){
    ret = ran_start(doubled, dim_DR)
    if (ret$loss < minloss){
      minloss = ret$loss
      optimAlloc = ret$alloc
      coord = ret$coord
    }
  }
  return(list("alloc"=optimAlloc, "coord"=coord))
}

ran_start = function(doubled, dim_DR){
  # In this function, we try to obtain the best possible solution given the random start
  random= floor(runif(500,min=1, max=5))
  ZK = as.matrix(dummy_cols(random, remove_selected_columns = TRUE))
  ZKprev = as.matrix(matrix(1:2000, nrow = 500, ncol = 4))
  am_var = dim(doubled)[2]/2
  iter=0
  # In this while loop, we iterate between K-means and CA until convergence is reached
  while((!identical(ZK,ZKprev)) && iter<300 ){ 
    ZKprev = ZK
    iter= iter+1
    ### obtain G and Z using CA
    N = as.matrix(doubled)
    ZKN = t(ZK) %*% N
    P = ZKN
    res.ca  <- CA(P, ncp=dim_DR, graph = FALSE)
    if (am_var==16){
      colnames(N) = c(1:32)
    }
    else if(am_var==32){
      colnames(N) = c(1:64)
    }else{
      colnames(N) = c(1:128)
    }
    G = res.ca$row$coord
    M = diag(500) - 1/500
    Dc = sqrt(solve(diag(colSums(N))))
    B = res.ca$col$coord
    ### Y = M N B
    Y = M %*% N %*% B
    # execute Km-means using Y and G
    Km = kmeans(Y, centers=G, algorithm="Lloyd", iter.max=20)
    ZK = Km$cluster
    ZK = as.matrix(dummy_cols(ZK, remove_selected_columns = TRUE))
    loss = Km$tot.withinss
    if(dim(ZK)[2]<4)
    {loss=Inf}
    my_list= list("alloc" = ZK, "loss" = loss, "coord"=Y)
  }
  return(my_list)
}




Simul = function(am_sim, path, dim_DR){
  #Run all simultaneous function and obtain average for all functions
  ARIcol_RKM = ARI_RKM(am_sim, path, dim_DR)
  ARIcol_FKM = ARI_FKM(am_sim, path, dim_DR)
  ARIcol_CCA=ARI_CCA(am_sim, path, dim_DR)
  ARIcol_MCA_kM = ARI_MCA_kM(am_sim, path, dim_DR)
  return(c( colMeans(ARIcol_RKM), colMeans(ARIcol_FKM), colMeans(ARIcol_CCA), colMeans(ARIcol_MCA_kM)))
}



Simul_Res = function(am_sim){
  ## Execute all simultaneous methods on all the simulated data
  Simul_Res = matrix( nrow =0 , ncol =12)
  dim_DR=3
  colnames(Simul_Res)=c("case", "am_var", "am_rat", "Std_dev", "RKM_ARI", "RKM_SIL", "FKM_ARI", "FKM_SIL", "CCA_ARI", "CCA_SIL", "MCAKM_ARI", "MCAKM_SIL")
   for (c in c(1, 2)){
    for (i in c(16, 32, 64)){
      for (j in c(3,5,9)){
        for (std in c("low", "high")){
          set.seed("1234")
          path = paste("C:\\Users\\siets\\Documents\\E&OR\\Thesis\\Proposal\\Simulatie_case",c,i,j,std,sep="_")
          dim_DR=3
          if (std=="low"){
            std_dev=1
          }else{
            std_dev=2
          }
          res = cbind(c, i, j, std_dev, t(as.numeric(Simul(am_sim, path, dim_DR))))
          Simul_Res = rbind(Simul_Res, res)
          ### Give an update of all the results
          print(Simul_Res)
        }
      }
    }
  }
  return(Simul_Res)
}



