library(sbgcop)
library(mixedgcImp)
library(clusterGeneration) # for generate random correlation matrices
library(mvtnorm)
# change the working directory to "Missing-Value-Imputation-for-Mixed-Data-via-Gaussian-Copula"
source('fun_auxi.R')
source('fun_eval.R')

misRatio = 0.3
nrep = 100 
n = 2000
p = 15
nIterSbgcop = c(100,500,1000,2000,4000,8000)
nIterEM = c(3,5,15,25,35,45)

res = array(0, dim = c(nrep,5,6,2), 
                 dimnames = list(NULL, 
                                 c('continuous', 'binary', 'ordinal', 'corr', 'time'),
                                 NULL, c('Copula-EM', 'sbgcop')))

for (i in 1:nrep){
  # generate data 
  gen = generateDataSim(n=2000,p=15,seed=10*i-1)
  X_true = gen$X
  Sigma = gen$Sigma
  
  # mask
  X_obs = mask(X_true, ratio = 0.3, seed = 10*i-1)$Xnew
  
  # copula-EM
  for (j in 1:length(nIterEM)){
    a = Sys.time()
    est = impute_mixedgc(X_obs, runiter = nIterEM[j])
    b = Sys.time()
    
    error = cal_mae_scaled(xhat = est$Ximp, xobs = X_obs, x = X_true)
    res[i,1,j,1] = mean(error[1:5])
    res[i,2,j,1] = mean(error[6:10])
    res[i,3,j,1] = mean(error[11:15])
    res[i,4,j,1] = norm(Sigma - est$R, type = 'F')/norm(Sigma, type = 'F')
    res[i,5,j,1] = difftime(b,a,units = 'secs')
  }
  
  # sbgcop
  for (j in 1:length(nIterSbgcop)){
    a = Sys.time()
    est = sbgcop.mcmc(Y=X_obs, verb = FALSE, nsamp = nIterSbgcop[j])
    b = Sys.time()
    xhat = est$Y.pmean
    for (m in 6:10) xhat[,m] = trunc.rating(xhat[,m],0,1)
    for (m in 11:15) xhat[,m] = trunc.rating(xhat[,m],1,5)
    
    error = cal_mae_scaled(xhat = xhat, xobs = X_obs, x = X_true, round = TRUE)
    res[i,1,j,2] = mean(error[1:5])
    res[i,2,j,2] = mean(error[6:10])
    res[i,3,j,2] = mean(error[11:15])
    res[i,4,j,2] = mean(apply(est$C.psamp, 3, function(x){norm(x-Sigma, type = 'F')}))/norm(Sigma, type = 'F')
    res[i,5,j,2] = difftime(b,a,units = 'secs')
  }
  
  print(paste('finish iteration ',i))
}

#apply(res, c(2,3,4), mean)
#save(res, file = "Res_figure4.RData") # store the results if needed
