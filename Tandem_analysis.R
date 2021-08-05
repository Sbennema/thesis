Tan_Res = function(am_sim){
  ## Tandem clustering
  library("factoextra")
  library("gplots")
  library("FactoMineR")
  library("ggfortify")
  library("pdfCluster")
  library("cluster")
  library(ncpen)
  
  dim_DR = 3
  ARI_PCA = function(am_sim, path, dim_DR){
    # Returns average ARI and ASW using tandem PCA in 50 simulations
    ARIcol = data.frame()
    SILcol = data.frame()
    for (i in c(1:am_sim)){
      num = i
      name = paste(num,"sim.csv", sep="_")
      pa_na = paste(path, name, sep="\\")
      matrix = read.csv(pa_na)
      dt <- as.table(as.matrix(matrix[,1:(dim(matrix)[2]-1)]))
      res.pca <- prcomp(dt, scale = FALSE)
      df = res.pca$x[,1:dim_DR]
      k2 <- kmeans(df, centers = 4, nstart = 25)
      ARI = adj.rand.index(c(k2$cluster), c(matrix[,dim(matrix)[2]]))
      sil = silhouette(c(k2$cluster), dist(dt))
      meansil = mean(sil[,3])
      ARIcol = rbind(ARIcol, ARI)
      SILcol = rbind(SILcol, meansil)
    }
    ARICIL = cbind(ARIcol, SILcol)
    return(ARICIL)
  }
  
  ARI_CA = function(am_sim, path, dim_DR){
    # Returns average ARI and ASW for Tandem CA in 50 simulations
    ARIcol = data.frame()
    SILcol = data.frame()
    for (i in c(1:am_sim)){
      num = i
      name = paste(num,"sim.csv", sep="_")
      pa_na = paste(path, name, sep="\\")
      matrix = read.csv(pa_na)
      dt <- as.table(as.matrix(matrix[,1:(dim(matrix)[2]-1)]))
      doubled = cbind(dt, max(dt)+1-dt)
      res.ca  <- CA(doubled, ncp=dim_DR, graph = FALSE)
      df = res.ca$row$coord
      k2 <- kmeans(df, centers = 4, nstart = 25)
      ARI = adj.rand.index(c(k2$cluster), c(matrix[,dim(matrix)[2]]))
      sil = silhouette(c(k2$cluster), dist(dt))
      meansil = mean(sil[,3])
      ARIcol = rbind(ARIcol, ARI)
      SILcol = rbind(SILcol, meansil)
    }
    ARICIL = cbind(ARIcol, SILcol)
    return(ARICIL)
  }
  
  
  ARI_MCA = function(am_sim, path, dim_DR){
    # Returns average ARI and ASW for Tandem MCA in 50 simulations
    ARIcol = data.frame()
    SILcol = data.frame()
    for (i in c(1:am_sim)){
      num = i
      name = paste(num,"sim.csv", sep="_")
      pa_na = paste(path, name, sep="\\")
      matrix = read.csv(pa_na)
      dt <- as.table(as.matrix(matrix[,1:(dim(matrix)[2]-1)]))
      # dt = data.frame(factor(dt[,1]), factor(dt[,2]), factor(dt[,3]), factor(dt[,4]))
      mat = c(1:(dim(matrix)[1]))
      for (i in c(1:(dim(matrix)[2]-1))){
        column = data.frame(factor(dt[,i]))
        mat = cbind(mat, column)
      }
      res.mca  <- MCA(mat[,2:dim(matrix)[2]], ncp=dim_DR, graph = FALSE)
      df = res.mca$ind$coord
      k2 <- kmeans(df, centers = 4, nstart = 25)
      ARI = adj.rand.index(c(k2$cluster), c(matrix[,dim(matrix)[2]]))
      sil = silhouette(c(k2$cluster), dist(dt))
      meansil = mean(sil[,3])
      ARIcol = rbind(ARIcol, ARI)
      SILcol = rbind(SILcol, meansil)
    }
    ARICIL = cbind(ARIcol, SILcol)
    return(ARICIL)
  }
  
  ARI_FDC = function(am_sim, path, dim_DR){
    # Returns average ARI and ASW for Full-dimensional clustering in 50 simulations
    ARIcol = data.frame()
    SILcol = data.frame()
    for (i in c(1:am_sim)){
      num = i
      name = paste(num,"sim.csv", sep="_")
      pa_na = paste(path, name, sep="\\")
      matrix = read.csv(pa_na)
      dt <- as.matrix(matrix[,1:(dim(matrix)[2]-1)])
      k2 <- kmeans(dt, centers = 4, nstart = 25)
      sil = silhouette(c(k2$cluster), dist(dt))
      meansil = mean(sil[,3])
      ARI = adj.rand.index(c(k2$cluster), c(matrix[,dim(matrix)[2]]))
      ARIcol = rbind(ARIcol, ARI)
      SILcol = rbind(SILcol, meansil)
    }
    ARICIL = cbind(ARIcol, SILcol)
    return(ARICIL)
  }
  
  tandem =function(am_sim, path, dim_DR){
    #Run all tandem functions and full dimensional clustering and obtain average for all functions
    ARIcol_PCA = ARI_PCA(am_sim, path, dim_DR)
    ARIcol_CA = ARI_CA(am_sim, path, dim_DR)
    ARIcol_MCA = ARI_MCA(am_sim, path, dim_DR)
    ARIcol_FDC = ARI_FDC(am_sim, path, dim_DR)
    return(c(colMeans(ARIcol_PCA),  colMeans(ARIcol_CA),  colMeans(ARIcol_MCA),  colMeans(ARIcol_FDC)))
  }  
  ## Execute all tandem methods and full dimensional clustering on all the simulated data
  Tandem_res = matrix( nrow =0 , ncol =12)
  colnames(Tandem_res)=c("case", "am_var", "am_rat", "Std_dev", "PCA_ARI", "PCA_SIL", "CA_ARI", "CA_SIL", "MCA_ARI", "MCA_SIL", "FDC_ARI", "FDC_SIL")
  for (c in c(1,2)){
    for (i in c(16, 32, 64)){
      for (j in c(3,5,9)){
        for (std in c("low", "high")){
        set.seed("1234")
        path = paste("C:\\Users\\siets\\Documents\\E&OR\\Thesis\\Proposal\\Simulatie_case",c,i,j,std,sep="_")
        if (std=="low"){
          std_dev=1
        }else{
          std_dev=2
        }
        res = cbind(c, i, j, std_dev, t(as.numeric(tandem(am_sim, path, dim_DR))))
        Tandem_res = rbind(Tandem_res, res)
        }
      }
    }
  }
  name =paste("Tandem_Res")
  write.csv(Tandem_res, file = name, row.names = FALSE)
  return(Tandem_res)
}




