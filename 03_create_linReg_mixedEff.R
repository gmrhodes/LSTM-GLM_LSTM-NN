# This file creates the linear regression vectors \zeta^{(S)}(\tau) 
# and the mixed effects vectors \zeta^{(M)}(\tau) for the simulations.
library(dplyr)
library(nlme)
library(lme4)
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
numCores = 25

#Directory path of folder containing data generated by 01_generate_data.R
##Must contain sub-directories "aft" and "cox"
input_path = "D:/Simulations/data/"

#Directory path of output folder to write results to
##Must contain sub-directories "aft" and "cox"
output_path = "D:/Simulations/linReg_mixedEff/"


##########################  FUNCTION TO MIN-MAX SCALE DATA ########################## 
#Vector of continuous variables: x
minMaxScale = function(x) {
  #Min-max scale data
  return( (x- min(x)) / (max(x)-min(x)) )
}


##########################  FUNCTION TO CREATE LINEAR REGRESSION & MIXED EFFECTS DATA SETS FOR GIVEN ITERATION ########################## 
#Survival time data set: timeDF
#Covariate data set: covDF
iterFN = function(timeDF, covDF) {
  library(nlme)
  library(lme4)
  set.seed(5678)
  
  #Merge survival/censoring time data & covariate data
  simDF = merge(timeDF, covDF, by="id")
  
  #Min-max scale longitudinal covariate
  simDF$long_cov = minMaxScale(simDF$long_cov)
  
  #Keep only measurements taken prior to t on patients at risk at t
  simDF_t = simDF[which((simDF$meas_time<t) & (simDF$Y>t)),]
  
  #Create restricted residual life variable
  simDF_t$restricted_residual_life = simDF_t$Y - t
  
  #Create linear regression data set shell
  lrDF = simDF_t[which(simDF_t$block==1),c("id","restricted_residual_life","delta_star","base_cov")]
  
  #Fit linear regression model on each patient
  lrFit = lmList(long_cov ~ meas_time | id, data=simDF_t[,c("id","meas_time","long_cov")])
  
  #Extract slopes & intercepts
  lrMat = matrix(NA, ncol=2, nrow=length(lrFit))
  for(j in 1:length(lrFit)) {
    lrMat[j,1] = lrFit[[j]]$coefficients[1]
    lrMat[j,2] = lrFit[[j]]$coefficients[2]
  }
  lrDF[,"long_interc"] = lrMat[,1]
  lrDF[,"long_slope"] = lrMat[,2]
  
  #Create mixed effects data set shell
  mixedDF = simDF_t[which(simDF_t$block==1),c("id","restricted_residual_life","delta_star","base_cov")]
  
  #Fit mixed effects model 
  lmeFit = lmer(long_cov ~ meas_time + (meas_time|id), data=simDF_t[,c("id","meas_time","long_cov")], REML=F, control=lmerControl(optimizer="bobyqa"))

  #Extract slopes & intercepts
  mixedDF[,c("long_interc","long_slope")] = coef(lmeFit)[[1]]
  
  return(list(lrDF, mixedDF))
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
myCluster = makeCluster(numCores, type="PSOCK") 
registerDoParallel(myCluster)

#Create linear regression & mixed effects data sets for each simulation
results = foreach(i=1:iter) %dopar% {iterFN(time_list[[i]], covDF)}
stopCluster(myCluster)

#Write data sets for each simulation
for(i in 1:iter) {
  write.table(results[[i]][[1]], paste0(output_path, surv_model, "/linReg",i,".txt"), sep=",", row.names=F)
  write.table(results[[i]][[2]], paste0(output_path, surv_model, "/mixedEff",i,".txt"), sep=",", row.names=F)
}


