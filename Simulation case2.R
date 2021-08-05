### Data simulation case 2

library(matrixStats)
library(csv)


dist_emo = function(emo, am_cat, am_rat, std_dev){
  # Function that returns value given the emotion, the rating scale, and whether the std_dev is high/low
  if (am_rat==3){
    mult_rat= 0.5
  }
  else if(am_rat==5){
    mult_rat=1
  }else{
    mult_rat=2
  }
    
  if(std_dev=="low"){
    #USed to bed std_dev = am_rat/4
    std_dev = am_rat/4
  } 
  else{
    std_dev = am_rat/3
  }
  # given the emotion, the data is simulated from the normal distribution
  if(emo=="neg"){
    a= rnorm(am_cat, (1+mult_rat), std_dev)
  }
  if(emo=="neu"){
    a= rnorm(am_cat, (1+am_rat)/2, std_dev*(2))
  }
  if (emo=="pos"){
    a= rnorm(am_cat, am_rat-mult_rat, std_dev)
  }
  # Data is rounded such we get rating data
  a = round(a)
  a[a<1] =1
  a[a>am_rat] = am_rat
  return(a)
}


sim_mat = function(am_var, am_rat, std_dev){
  # Funtion that simulates 500 obsevrations, where the observations belong to one of the clusters
  # Amount of variables, rating scale and high or low standard deviation is dependent on the input
  n_obs=500
  mat = data.frame()
  for (i in c(1:4)){
    for (j in c(1:(n_obs/4))){
      ob_i = cl_i(i, am_var, am_rat, std_dev)
      mat = rbind(mat, ob_i)
    }}
  return(mat)
}

cl_i = function(clus_i, am_var, am_rat, std_dev){
  # Given the cluster of an observation, the variables are simulated
  am_relev = (6+am_var/8)/2
  am_neutral = am_var-2*am_relev
  if (clus_i==1){
  row_v= c(dist_emo("pos", am_relev, am_rat, std_dev), dist_emo("neg", am_relev, am_rat, std_dev), dist_emo("neu", am_neutral, am_rat, std_dev))
  }
  if (clus_i==2){
    row_v= c(dist_emo("neg", am_relev, am_rat, std_dev), dist_emo("pos", am_relev, am_rat, std_dev), dist_emo("neu", am_neutral, am_rat, std_dev))
  }
  if (clus_i==3){
    row_v= c(dist_emo("pos", am_relev, am_rat, std_dev), dist_emo("pos", am_relev, am_rat, std_dev), dist_emo("neu", am_neutral, am_rat, std_dev))
  }
  if (clus_i==4){
    row_v= c(dist_emo("neg", am_relev, am_rat, std_dev), dist_emo("neg", am_relev, am_rat, std_dev), dist_emo("neu", am_neutral, am_rat, std_dev))
  }
  
  val = cbind(t(c(row_v)), clus_i)
  return(val)
}

er <- function(x, n = 1) {
  ## makes sure the simulations are different for the different functions
  if (n == 0) x else c(tail(x, -n), head(x, n))
}



write_file = function(am_sim, am_var, am_rat, std_dev){
  ### given the settings, the files are generated to this path
  path = paste("C:\\Users\\siets\\Documents\\E&OR\\Thesis\\Proposal\\Simulatie_case_2",am_var,am_rat,std_dev,sep="_")
  dir.create(path)
  for (i in c(1:am_sim))
  {
    sample = sim_mat(am_var, am_rat, std_dev)
    # print(dim(sample))
    # colMeans(sample)
    # colSds(data.matrix(sample))
    file = paste(i, "sim.csv", sep="_")
    pa_file= paste(path, file, sep="\\")
    write.csv(sample, file = pa_file, row.names = FALSE)
  }
  return(print("done"))
}


Sim_case_2 = function(am_sim){
  ### simulates all datasets with the different settings, am_sim times
  for (i in c(3,5,9)){
    for (j in c(16, 32, 64)){
      for (std in c("low", "high")){
        set.seed("1234")
        sample = write_file(am_sim, am_var=j, am_rat=i, std_dev=std)
      }
    }
  }
  return(print("case_2_done"))
}
