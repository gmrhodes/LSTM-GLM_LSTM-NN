# This file creates the FPCA vectors \zeta^{(F)}(\tau) for the simulations.
library(dplyr)
library(fdapace)
library(foreach)
library(doParallel)
set.seed(5678)


############################### USER SET-UP ##############################
#Survival model ("aft" or "cox")
surv_model = "aft"

#Prediction time
t = 5

#Number of iterations
iter = 500

#Number of cores to parallelize across
numCores = 15

#Directory path of folder containing data generated by 01_generate_data.R
##Must contain sub-directories "aft" and "cox"
input_path = "D:/Simulations/data/"

#Directory path of output folder to write results to
##Must contain sub-directories "aft" and "cox"
output_path = "D:/Simulations/fpca/"


##########################  FUNCTION TO MIN-MAX SCALE DATA ########################## 
#Vector of continuous variables: x
minMaxScale = function(x) {
  #Min-max scale data
  return( (x- min(x)) / (max(x)-min(x)) )
}


##########################  FUNCTION TO CREATE FPCA DATA SET FOR GIVEN ITERATION ########################## 
#Survival time data set: timeDF
#Covariate data set: covDF
iterFN = function(timeDF, covDF, i) {
  library(fdapace)
  set.seed(5678)
  print(paste("Start:",i))
  
  #Merge survival/censoring time data & covariate data
  simDF = merge(timeDF, covDF, by="id")
  
  #Min-max scale longitudinal covariate
  simDF$long_cov = minMaxScale(simDF$long_cov)
  
  #Keep only measurements taken prior to t on patients at risk at t
  simDF_t = simDF[which((simDF$meas_time<t) & (simDF$Y>t)),]
  
  #Create restricted residual life variable
  simDF_t$restricted_residual_life = simDF_t$Y - t
  
  #Create dataframe of IDs, response variables, and baseline covariates
  fpcaDF = simDF_t[which(simDF_t$block==1),c("id","restricted_residual_life","delta_star","base_cov")]
  
  #Obtain FPC scores via PACE algorithm
  covList = split(simDF_t$long_cov, simDF_t$id)
  tList = split(simDF_t$meas_time, simDF_t$id)
  fpca_fit = FPCA(covList, tList, list(dataType='Sparse'))
  scores = fpca_fit$xiEst
  
  #Create FPC scores data set
  colnames(scores) = paste0("long_cov_", c(1:ncol(scores)))
  fpcaDF = cbind(fpcaDF, scores)
  
  print(paste("End:",i))
  return(fpcaDF)
}


############################### MAIN ##############################
#Input covariate data
covDF = read.table(paste0(input_path, "covariates.csv"), sep=",", header=T)

#Input survival time data for each simulation
time_list = vector(mode="list", length=iter)
for(i in 1:iter) {
  time_list[[i]] = read.table(paste0(input_path, surv_model, "/survTimes", i, ".csv"), sep=",", header=T)
}

#Set up parallelization
myCluster = makeCluster(numCores, type="PSOCK", outfile="") 
registerDoParallel(myCluster)

#Create FPC data set for each simulation
results = foreach(i=1:iter) %dopar% {iterFN(time_list[[i]], covDF, i)}
stopCluster(myCluster)

#Write data sets of FPC scores for each simulation
for(i in 1:iter) {
  write.table(results[[i]], paste0(output_path, surv_model, "/fpca",i,".txt"), sep=",", row.names=F)
}

