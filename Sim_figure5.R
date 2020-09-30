library(mixedgcImp)
library(clusterGeneration)
library(mvtnorm)
library(missMDA)
library(missForest)
library(xpcaR) # available at https://gitlab.com/xpca/xpcar
# change the working directory to "Missing-Value-Imputation-for-Mixed-Data-via-Gaussian-Copula"
source('fun_auxi.R')


ratio = 0.1 * 1:5
nrep = 100 
n = 2000
p = 15

res = array(0, dim = c(nrep,5,3,4), 
            dimnames = list(NULL, NULL,
                            c('continuous', 'binary', 'ordinal'),
                            c('Copula-EM', 'missForest', 'imputeFAMD', 'xPCA')))

for (i in 1:nrep){
  # generate data 
  gen = generateDataSim(n=n,p=p,seed=10*i-1)
  X_true = gen$X
  
  for (j in 1:length(ratio)){
    # mask
    X_obs = mask(X_true, ratio = ratio[j], seed = 10*i-1)$Xnew
    
    # copula EM
    est = impute_mixedgc(X_obs)
    error = cal_mae_scaled(xhat = est$Ximp, xobs = X_obs, x = X_true)
    res[i,j,1,1] = mean(error[1:5])
    res[i,j,2,1] = mean(error[6:10])
    res[i,j,3,1] = mean(error[11:15])
    
    # missForest
    x = as.data.frame(X_obs)
    for (m in 6:15) x[,m] = as.factor(x[,m])
    est = missForest(xmis = x, verbose = FALSE)
    error = cal_mae_scaled(xhat = est$ximp, xobs = X_obs, x = X_true)
    res[i,j,1,2]  = mean(error[1:5])
    res[i,j,2,2]  = mean(error[6:10])
    res[i,j,3,2]  = mean(error[11:15])
    
    # imputeFAMD
    x = as.data.frame(X_obs)
    for (m in 6:15) x[,m] = as.factor(X_obs[,m])
    est = try(imputeFAMD(x, ncp = 6))
    c = 1
    while(class(est) == 'try-error'){
      est = try(imputeFAMD(x, ncp = 6-c))
      c = c+1
    }
    xhat = est$completeObs
    for (m in 6:15) xhat[,m] = map_ordinal_FAMD(xhat[,m])
    error = cal_mae_scaled(xhat = xhat, xobs = X_obs, xtrue = X_true)
    res[i,j,1,3] = mean(error[1:5])
    res[i,j,2,3] = mean(error[6:10])
    res[i,j,3,3] = mean(error[11:15])
    
    # xPCA
    est = try(xpca(X_obs, rank = 3))
    c = 1
    while(class(est) == 'try-error'){
      est = try(xpca(X_obs, rank = 3-c))
      c = c+1
    }
    xhat = est$fittedEsts
    for (m in 6:10) xhat[,m] = trunc.rating(xhat[,m],0,1)
    for (m in 11:15) xhat[,m] = trunc.rating(xhat[,m],1,5)
    
    error = cal_mae_scaled(xhat = xhat, xobs = X_obs, xtrue = X_true, round = TRUE)
    res[i,j,1,4] = mean(error[1:5])
    res[i,j,2,4] = mean(error[6:10])
    res[i,j,3,4] = mean(error[11:15])
  }

  print(paste('finish iteration ',i))
}

#apply(res, c(2,3,4), mean)
save(res, file = "Res_figure5.RData") # store the results if needed
