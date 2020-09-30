# functions supporting evaluation on datasets LEV, ESL, GBSG, TIPS
evalOrdinal = function(x, ratio, rep = 20, parallel = FALSE, verbose = FALSE, rankxPCA=1, rankFAMD=5){
  x = as.matrix(x)
  p = dim(x)[2]
  n = dim(x)[1]
  if (parallel){
    Err = list()
    Err = foreach(i = 1:rep)%dopar%{
      erri = array(0, dim = c(4,2))
      xobs = mask(x, ratio = ratio, seed.start = i*5)$Xnew
      
      # Copula
      est = impute_mixedgc(xobs)
      e = cal_mae_scaled(est$Ximp, xobs, x)
      erri[1,1] = e[p]
      erri[1,2] = mean(e[1:(p-1)])
      
      # missForest
      z = as.data.frame(xobs)
      for (s in 1:p) z[,s] = as.factor(xobs[,s])
      est = missForest(z)
      e = cal_mae_scaled(est$ximp, xobs, x)
      erri[2,1] = e[p]
      erri[2,2] = mean(e[1:(p-1)])
      
      # xPCA
      est = xpca(xobs, rank = rankxPCA)
      xhat = round(est$fittedEsts)
      for (s in 1:p){
        xhat[,s] = trunc.rating(xhat[,s], xmin=min(xobs[,s], na.rm = TRUE), xmax=max(xobs[,s], na.rm = TRUE))
      }
      e = cal_mae_scaled(xhat, xobs, x)
      erri[3,1] = e[p]
      erri[3,2] = mean(e[1:(p-1)])
      
      # imputeFAMD, ncp =5 for ESL and LEV
      ncp = rankFAMD
      est = try(imputeFAMD(z, ncp = ncp))
      c = 1
      while(class(est) == 'try-error'){
        est = try(imputeFAMD(z, ncp = ncp-c))
        c = c+1
      }
      xhat = est$completeObs
      for (s in 1:p) xhat[,s] = map_ordinal_FAMD(xhat[,s])
      erri[4,1] = e[p]
      erri[4,2] = mean(e[1:(p-1)])

      erri
    }
    err = array(0, dim = c(rep,3,2))
    for (i in 1:rep) err[i,,] = Err[[i]]
  }else{
    err = array(0, dim = c(rep,4,2))
    for (i in 1:rep){
      xobs = mask(x, ratio = ratio, seed.start = i*5)$Xnew
      
      # Copula EM
      est = impute_mixedgc(xobs)
      e = cal_mae_scaled(est$Ximp, xobs, x)
      err[i,1,1] = e[p]
      err[i,1,2] = mean(e[1:(p-1)])
      
      # missForest
      z = as.data.frame(xobs)
      for (s in 1:p) z[,s] = as.factor(xobs[,s])
      est = missForest(z)
      e = cal_mae_scaled(est$ximp, xobs, x)
      err[i,2,1] = e[p]
      err[i,2,2] = mean(e[1:(p-1)])
      
      # xPCA
      est = xpca(xobs, rank = rankxPCA)
      xhat = est$fittedEsts
      for (s in 1:p){
        xhat[,s] = trunc.rating(xhat[,s], xmin=min(xobs[,s], na.rm = TRUE), xmax=max(xobs[,s], na.rm = TRUE))
      }
      e = cal_mae_scaled(xhat, xobs, x)
      err[i,3,1] = e[p]
      err[i,3,2] = mean(e[1:(p-1)])
      
      # imputeFAMD
      ncp = rankFAMD
      est = try(imputeFAMD(z, ncp = ncp))
      c = 1
      while(class(est) == 'try-error'){
        est = try(imputeFAMD(z, ncp = ncp-c))
        c = c+1
      }
      xhat = est$completeObs
      for (s in 1:p) xhat[,s] = map_ordinal_FAMD(xhat[,s])
      e = cal_mae_scaled(xhat, xobs, x)
      err[i,4,1] = e[p]
      err[i,4,2] = mean(e[1:(p-1)])
      
      if (verbose) print(paste('finish iteration',i))
    }
  }
  err
}

evalMixed = function(x, ratio, rep = 20, nlevels=20, parallel = FALSE, verbose = FALSE, rankxPCA=1, rankFAMD=5){
  x = as.matrix(x)
  p = dim(x)[2]
  n = dim(x)[1]
  if (parallel){
    Err = list()
    Err = foreach(i = 1:rep) %dopar% {
      erri = array(0, dim = c(4,2))
      xobs = mask(x, ratio = ratio, seed.start = i*5)$Xnew
      l = levels.percol(xobs)
      ind.ordinal = which(l <= nlevels)
      
      # Copula
      est = impute_mixedgc(xobs)
      erri[1,] = cal_error_mixed(est$Ximp, xobs, x, ordinal = ind.ordinal)
      
      # missForest
      z = as.data.frame(xobs)
      for (s in ind.ordinal) z[,s] = as.factor(xobs[,s])
      est = missForest(z)
      erri[2,] = cal_error_mixed(est$ximp, xobs, x, ordinal = ind.ordinal)
      
      # xPCA
      est = xpca(xobs, rank = rankxPCA)
      xhat = est$fittedEsts
      for (s in ind.ordinal){
        xhat[,s] = trunc.rating(xhat[,s], xmin=min(xobs[,s], na.rm = TRUE), xmax=max(xobs[,s], na.rm = TRUE))
      }
      erri[3,] = cal_error_mixed(xhat, xobs, x, ordinal = ind.ordinal)
      
      
      # imputeFAMD
      ncp = rankFAMD
      est = try(imputeFAMD(z, ncp = ncp))
      c = 1
      while(class(est) == 'try-error'){
        est = try(imputeFAMD(z, ncp = ncp-c))
        c = c+1
      }
      xhat = est$completeObs
      for (s in ind.ordinal) xhat[,s] = map_ordinal_FAMD(xhat[,s])
      erri[4,] = cal_error_mixed(xhat, xobs,x, ordinal = ind.ordinal)
    
      erri
    }
    err = array(0, dim = c(rep,4,2))
    for (i in 1:rep) err[i,,] = Err[[i]]
  }else{
    err = array(0, dim = c(rep,4,2))
    for (i in 1:rep){
      xobs = mask(x, ratio = ratio, seed.start = i*5)$Xnew
      l = levels.percol(xobs)
      ind.ordinal = which(l <= nlevels)
      
      # Copula
      est = impute_mixedgc(xobs)
      err[i,1,] = cal_error_mixed(est$Ximp, xobs, x, ordinal = ind.ordinal)
      
      # missForest
      z = as.data.frame(xobs)
      for (s in ind.ordinal) z[,s] = as.factor(xobs[,s])
      est = missForest(z)
      err[i,2,] = cal_error_mixed(est$ximp, xobs, x, ordinal = ind.ordinal)
      
      # xPCA
      est = try(xpca(xobs, rank = rankxPCA))
      c = 1
      while(class(est) == 'try-error'){
        est = try(xpca(xobs,rank = rankxPCA-c))
        c = c+1
      }
      xhat = est$fittedEsts
      for (s in ind.ordinal){
        xhat[,s] = trunc.rating(xhat[,s], xmin=min(xobs[,s], na.rm = TRUE), xmax=max(xobs[,s], na.rm = TRUE))
      }
      err[i,3,] = cal_error_mixed(xhat, xobs, x, ordinal = ind.ordinal)
      
      
      # imputeFAMD
      ncp = rankFAMD
      est = try(imputeFAMD(z, ncp = ncp))
      c = 1
      while(class(est) == 'try-error'){
        est = try(imputeFAMD(z, ncp = ncp-c))
        c = c+1
      }
      xhat = est$completeObs
      for (s in ind.ordinal) xhat[,s] = map_ordinal_FAMD(xhat[,s])
      err[i,4,] = cal_error_mixed(xhat, xobs,x, ordinal = ind.ordinal)
      
      if (verbose) print(paste('finish iteration',i))
    }
  }
  
  err
}

cal_error_mixed = function(xhat, xobs, x, ordinal){
  err = numeric(2)
  p = dim(xhat)[2]

  e = cal_mae_scaled(xhat, xobs, x)
  err[1] = mean(e[ordinal])
  err[2] = mean(e[setdiff(1:p,ordinal)])
  err
}


# evaluation using sbgcop
evalOrdinalSbgcop = function(x, iter, ratio, rep = 20){
  x = as.matrix(x)
  p = dim(x)[2]
  n = dim(x)[1]
  
  Err = foreach(i = 1:rep, .combine = cbind)%dopar%{
    erri = numeric(2)
    xobs = mask(x, ratio = ratio, seed.start = i*5)$Xnew
    
    # Copula
    est = sbgcop.mcmc(xobs, nsamp = iter, verbose=FALSE)
    xhat = est$Y.pmean
    for (s in 1:p){
      xhat[,s] = trunc.rating(xhat[,s], xmin=min(xobs[,s], na.rm = TRUE), xmax=max(xobs[,s], na.rm = TRUE))
    }
    
    e = cal_mae_scaled(xhat, xobs, x)
    erri[1] = e[p]
    erri[2] = mean(e[1:(p-1)])
    
    erri
  }
  
  Err
}

evalMixedSbgcop = function(x, iter, ratio, rep = 20, nlevels=20){
  x = as.matrix(x)
  p = dim(x)[2]
  n = dim(x)[1]
  
  Err = foreach(i = 1:rep, .combine = cbind) %dopar% {
    xobs = mask(x, ratio = ratio, seed.start = i*5)$Xnew
    l = levels.percol(xobs)
    ind.ordinal = which(l <= nlevels)
    
    est = sbgcop.mcmc(xobs, nsamp = iter, verbose = FALSE)
    xhat = est$Y.pmean
    for (s in ind.ordinal){
      xhat[,s] = trunc.rating(xhat[,s], xmin=min(xobs[,s], na.rm = TRUE), xmax=max(xobs[,s], na.rm = TRUE))
    }
    
    erri = cal_error_mixed(xhat, xobs, x, ordinal = ind.ordinal)
    erri
  }
  
  Err
}


