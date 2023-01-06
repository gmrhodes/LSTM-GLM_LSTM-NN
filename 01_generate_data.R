#This file generates the data sets for the simulations.
library(MASS)
set.seed(5678)

################################################### FUNCTIONS ################################################### 
#Function to compute longitudinal biomarker 
#Vector of observed measurement times: t
#Fixed intercept of piecewise constant function: a
#Vector of 10 fixed slopes of piecewise constant function: cVec
#Matrix of 11 random effects for m patients: randeff
Bt = function(t, a, cVec, randeff) {
  a + randeff[,1] + (cVec[1]+randeff[,2])*(t) +
    as.integer(t-1>0)*(cVec[2]+randeff[,3])*(t-1) +
    as.integer(t-2>0)*(cVec[3]+randeff[,4])*(t-2) +
    as.integer(t-3>0)*(cVec[4]+randeff[,5])*(t-3) +
    as.integer(t-4>0)*(cVec[5]+randeff[,6])*(t-4) +
    as.integer(t-5>0)*(cVec[6]+randeff[,7])*(t-5) +
    as.integer(t-6>0)*(cVec[7]+randeff[,8])*(t-6) +
    as.integer(t-7>0)*(cVec[8]+randeff[,9])*(t-7) +
    as.integer(t-8>0)*(cVec[9]+randeff[,10])*(t-8)  +
    as.integer(t-9>0)*(cVec[10]+randeff[,11])*(t-9) 
}


#Function to generate baseline covariates & longitudinal biomarkers 
#Sequence of fixed times: tSeq 
#Sample size: m
#Fixed intercept of piecewise constant function: a
#Vector of 10 fixed slopes of piecewise constant function: cVec
gen_covs = function(tSeq, m, a, cVec) {
  #Generate baseline covariates 
  xVec = runif(m, 0, 1) 
  
  #Generate random effects from N(0,D)
  ##Define the variance matrix
  sigma = matrix(0, 11, 11)
  diag(sigma) = c(1, rep(0.05,5), rep(0.01,5))
  ##Define the mean vector
  mu = rep(0,11)
  ##Draw m random effects
  b_i = mvrnorm(m, mu, sigma)
  
  #Generate longitudinal biomarkers
  ##Matrix to store longitudinal biomarkers with measurement errors
  bioMat = matrix(NA, nrow=m, ncol=length(tSeq))
  colnames(bioMat) = paste0("meas_",c(1:length(tSeq)))
  ##Matrix to store longitudinal biomarkers without measurement errors
  BtMat = matrix(NA, nrow=m, ncol=length(tSeq))
  colnames(BtMat) = paste0("val_",c(1:length(tSeq)))
  ##Matrix to store measurement times  
  tMat = matrix(NA, nrow=m, ncol=length(tSeq))
  colnames(tMat) = paste0("t_",c(1:length(tSeq)))
  
  ##For each fixed time
  for(i in 1:length(tSeq)) {
    ##Generate measurement errors
    err = rnorm(m, 0, 0.5)
    
    ##Generate jittered measurement times
    time_jitter = rnorm(m, 0, 0.05)
    tMat[,i] = time_jitter + tSeq[i]
    tMat[which(tMat[,i]<0),i] = 0
    
    ##Generate values from requested model
    BtMat[,i] = Bt(tMat[,i], a, cVec, b_i)
    bioMat[,i] = BtMat[,i] + err
  }
  
  #Return results
  res = list(xVec, bioMat, tMat, b_i, BtMat)
  names(res) = c("base_cov", "long_cov", "meas_time", "rand_effs", "err_free_long")
  return( res )
}


#Function to compute nu_K(t)
#Scalar interval: K
#Scalar time: t
#Fixed intercept of piecewise constant function: a
#Vector of 10 fixed slopes of piecewise constant function: cVec
#Survival time model parameters: beta1, beta2, lambda (lambda=1 for AFT)
#Matrix of 11 random effects for m patients: b_i
#Vector of m baseline covariates: xVec
nu_K_t = function(K, t, a, cVec, beta1, beta2, lambda=1, b_i, xVec) {
  #Compute summations
  sum1 = 0
  for(j in 1:K) {
    sum1 = sum1 + (cVec[j]+b_i[,j+1])*(j-1)
  }
  sum2 = 0
  for(j in 1:K) {
    sum2 = sum2 + (cVec[j]+b_i[,j+1])
  }
  
  #Compute multiplier
  mult = lambda*exp(beta2*xVec + beta1*(a+b_i[,1]) - beta1*sum1) / (beta1*sum2)
  
  #Compute exp terms
  term1 = exp(t*beta1*sum2)
  term2 = exp((K-1)*beta1*sum2)
  
  #Compute nu_K(t)
  return( mult*(term1-term2) )
}


#Function to compute nu(t)
#Scalar time: t
#Fixed intercept of piecewise constant function: a
#Vector of 10 fixed slopes of piecewise constant function: cVec
#Survival time model parameters: beta1, beta2, lambda (lambda=1 for AFT)
#Matrix of 11 random effects for patients: b_i
#Vector of m baseline covariates: xVec
nu_t = function(t, a, cVec, beta1, beta2, lambda=1, b_i, xVec) {
  #Identify interval
  if(t<=1) {
    K=1
  } else if (t>1 & t<=2) {
    K=2
  } else if (t>2 & t<=3) {
    K=3
  } else if (t>3 & t<=4) {
    K=4
  } else if (t>4 & t<=5) {
    K=5
  } else if (t>5 & t<=6) {
    K=6
  } else if (t>6 & t<=7) {
    K=7
  } else if (t>7 & t<=8) {
    K=8
  } else if (t>8 & t<=9) {
    K=9
  } else if (t>9) {
    K=10
  }
  
  #Compute nu(t)
  sum = 0
  if(K>1) {
    for(j in 1:(K-1)) {
      sum = sum + nu_K_t(j, j, a, cVec, beta1, beta2, lambda, b_i, xVec)
    }
  }
  nu_K_t_curr = nu_K_t(K, t, a, cVec, beta1, beta2, lambda, b_i, xVec)
  return( sum + nu_K_t_curr )
}


#Function to compute inverse of nu_i(t)
#Scalar nu for patient i: nu
#Fixed intercept of piecewise constant function: a
#Vector of 10 fixed slopes of piecewise constant function: cVec
#Survival time model parameters: beta1, beta2, lambda (lambda=1 for AFT)
#Vector of 11 random effects for patient i: b_i_vec
#Scalar baseline covariate for patient i: xVec_i
#Vector of 9 knots for patient i: knot_i
nu_t_inv = function(nu, a, cVec, beta1, beta2, lambda=1, b_i_vec, xVec_i, knot_i) {
  #Identify interval 
  if(nu<=knot_i[1]) {
    K=1
  } else if (nu>knot_i[1] & nu<=knot_i[2]) {
    K=2
  } else if (nu>knot_i[2] & nu<=knot_i[3]) {
    K=3
  } else if (nu>knot_i[3] & nu<=knot_i[4]) {
    K=4
  } else if (nu>knot_i[4] & nu<=knot_i[5]) {
    K=5
  } else if (nu>knot_i[5] & nu<=knot_i[6]) {
    K=6
  } else if (nu>knot_i[6] & nu<=knot_i[7]) {
    K=7
  } else if (nu>knot_i[7] & nu<=knot_i[8]) {
    K=8
  } else if (nu>knot_i[8] & nu<=knot_i[9]) {
    K=9
  } else if (nu>knot_i[9]) {
    K=10
  }
  
  #Compute inverse of nu_i(t)
  ##Compute summations
  b_i_mat = matrix(b_i_vec, nrow=1, ncol=11)
  sum0 = 0
  if(K>1) {
    for(j in 1:(K-1)) {
      sum0 = sum0 + nu_K_t(j, j, a, cVec, beta1, beta2, lambda, b_i_mat, xVec_i)
    }
  }
  sum1 = 0
  for(j in 1:K) {
    sum1 = sum1 + (cVec[j]+b_i_vec[j+1])*(j-1)
  }
  sum2 = 0
  for(j in 1:K) {
    sum2 = sum2 + (cVec[j]+b_i_vec[j+1])
  }
  
  ##Compute terms
  term1 = ((nu-sum0)*beta1*sum2) / (lambda*exp(beta2*xVec_i + beta1*(a+b_i_vec[1]) - beta1*sum1))
  term2 = exp((K-1)*beta1*sum2)
  
  ##Compute inverse
  nu_inv = log(term1+term2) / (beta1*sum2)
  
  return( nu_inv )
}


#Function to generate survival times from a
#Cox proportional hazards model with an exponential baseline hazard
#Sample size: m
#Fixed intercept of piecewise constant function: a
#Vector of 10 fixed slopes of piecewise constant function: cVec
#AFT model parameters: beta1, beta2, lambda
#Matrix of 11 random effects for m patients: b_i
#Vector of m baseline covariates: xVec
gen_cox = function(m, a, cVec, beta1, beta2, lambda, b_i, xVec) {
  #Compute knots of linear piecewise function for each patient
  knotMat = matrix(NA, nrow=m, ncol=9)
  for(j in 1:9) {
    knotMat[,j] = nu_t(j, a, cVec, beta1, beta2, lambda, b_i, xVec)
  }
  
  #Generate survival times
  ##Generate random uniform variables
  uVec = runif(m, 0, 1) 
  hVec = -log(1-uVec)
  
  ##Compute survival times
  surv_time = rep(NA, m)
  for(i in 1:m) {
    surv_time[i] = nu_t_inv(hVec[i], a, cVec, beta1, beta2, lambda, b_i[i,], xVec[i], knotMat[i,])
  }
  surv_time[which(is.na(surv_time))] = 9999
  
  #Generate censoring times 
  cens_time = runif(m, 0, 100)
  evt_times = cbind(surv_time, cens_time)
  
  #Return results
  return( evt_times )
}


#Function to generate survival times from an accelerated failure time model
#Sample size: m
#Fixed intercept of piecewise constant function: a
#Vector of 10 fixed slopes of piecewise constant function: cVec
#AFT model parameters: beta1, beta2
#Matrix of 11 random effects for m patients: b_i
#Vector of m baseline covariates: xVec
gen_aft = function(m, a, cVec, beta1, beta2, b_i, xVec) {
  #Compute knots of linear piecewise function for each patient
  knotMat = matrix(NA, nrow=m, ncol=9)
  for(j in 1:9) {
    knotMat[,j] = nu_t(j, a, cVec, beta1, beta2, 1, b_i, xVec)
  }
  
  #Generate survival times
  ##Generate random normal variables
  epsVec = rnorm(m, 3, 1)
  nuVec = exp(epsVec)
  
  ##Compute survival times
  surv_time = rep(NA, m)
  for(i in 1:m) {
    surv_time[i] = nu_t_inv(nuVec[i], a, cVec, beta1, beta2, 1, b_i[i,], xVec[i], knotMat[i,])
  }
  surv_time[which(is.na(surv_time))] = 9999
  
  #Generate censoring times 
  cens_time = runif(m, 0, 100)
  evt_times = cbind(surv_time, cens_time)
  
  #Return results
  return( evt_times )
}


#Function to re-shape dataframe from wide-format to long-format
#List returned from gen_covs: covs
reshape_w2l = function(covs) {
  #Bind measurement times, longitudinal covariates, & baseline covariates
  df = as.data.frame(cbind(covs$long_cov, covs$meas_time))
  df$base_cov = covs$base_cov
  
  #Add id column
  df$id = c(1:nrow(df))
  
  #Put longitudinal covariates in long format
  longDF_meas = reshape(df, varying=names(df)[startsWith(names(df),"meas_")], v.names="long_cov",
                        timevar="block", times=c(1:length(tSeq)), direction="long")
  
  #Put measurement times in long format
  longDF_t = reshape(df, varying=names(df)[startsWith(names(df),"t_")], v.names="meas_time",
                     timevar="block", times=c(1:length(tSeq)), direction="long")
  
  #Merge & sort
  longDF = merge(longDF_meas[,c("id","block","base_cov","long_cov")],
                 longDF_t[,c("id","block","meas_time")], by=c("id","block"))
  longDF = longDF[order(longDF$id, longDF$block),]
  return(longDF)
}


#Function to create analysis variables
#Matrix of survival/censoring times from gen_aft or gen_cox: timeMat
#Restricted lifetime: L
analysis_vars = function(timeMat, L) {
  #Add id 
  timeDF = as.data.frame(timeMat)
  timeDF$id = c(1:nrow(timeDF))
  
  #Create analysis variables
  timeDF$U = pmin(timeDF$surv_time, timeDF$cens_time) #U=min(T,C)
  timeDF$delta = ifelse(timeDF$surv_time<=timeDF$cens_time, 1, 0) #delta=I(T<=C)
  timeDF$Y = pmin(timeDF$U, L) #Y=min(T,C,L)
  timeDF$delta_L = ifelse(L<=timeDF$U, 1, 0)
  timeDF$delta_star = timeDF$delta + timeDF$delta_L*(1-timeDF$delta) #delta*=I(min(T,L)<=C)
  
  return(timeDF)
}


################################################### MAIN ################################################### 
#Define output folder where data will be written to
output_path = "D:/Simulations/data/"

#Define number of simulation iterations
iter = 500

#Define total sample size (training + testing)
m = 5000

#Define sequence of fixed times for measurements
tSeq = seq(0, 9, by=0.5)

#Define restricted lifetime
L = 50

#Define survival time model parameters
beta1 = 1
beta2 = 1
lambda = 0.05

#Generate covariates
a = -2
cVec = c(4,-7,5,-2.5,3.5,-5,1.5,2,-2,1)
covs = gen_covs(tSeq, m, a, cVec)

#Reshape covariates to long-format and write
covs_long = reshape_w2l(covs) 
write.table(covs_long, paste0(output_path, "covariates.csv"), sep=",", row.names=F)

#Create lists to store survival/censoring times generated by each model
aft = vector(mode="list", length=iter)
cox = vector(mode="list", length=iter)

#Generate m survival/censoring times for iter iterations from covs
for(i in 1:iter) {
  #Generate and write survival times from accelerated failure time model
  aft_time = gen_aft(m, a, cVec, beta1, beta2, covs$rand_effs, covs$base_cov)
  aft[[i]] = analysis_vars(aft_time, L)
  write.table(aft[[i]], paste0(output_path, "aft/survTimes", i, ".csv"), sep=",", row.names=F)
  
  #Generate and write survival times from Cox model
  cox_time = gen_cox(m, a, cVec, beta1, beta2, lambda, covs$rand_effs, covs$base_cov)
  cox[[i]] = analysis_vars(cox_time, L)
  write.table(cox[[i]], paste0(output_path, "cox/survTimes", i, ".csv"), sep=",", row.names=F)
}

#Save survival/censoring time data as RDS
saveRDS(aft, paste0(output_path,"aft/survTimes.RDS"))
saveRDS(cox, paste0(output_path,"cox/survTimes.RDS"))
