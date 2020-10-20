mask = function(X, ratio, seed.start = 1){
  n = dim(X)[1]
  p = dim(X)[2]
  loc.mis.old = which(is.na(X))
  num = n*p - length(loc.mis.old)
  seed = seed.start
  i= 0
  empty_row = TRUE
  while(empty_row){
    set.seed(seed)
    l = sample(1:num, num*ratio)
    
    Xnew = X
    if (length(loc.mis.old)>0){
      Xnew[-loc.mis.old][l] = NA
      loc = setdiff(which(is.na(Xnew)), loc.mis.old)
    } else{
      Xnew[l] = NA
      loc = l
    }
    
    empty_row = any(apply(Xnew, 1, function(x){sum(!is.na(x))}) == 0)
    i=i+1
    seed=seed+1
    if (i > 100) stop('cannot produce masking without empty row')
  }
  
  loc.obs = setdiff(which(!is.na(X)), loc)
  return(list(Xnew = Xnew, locMis = loc, seed=seed-1, locObs = loc.obs))
}

generateDataSim = function(n,p,seed){
  set.seed(seed)
  Sigma = rcorrmatrix(d=p)
  Z_true = rmvnorm(n=n, mean = rep(0,p),  sigma = Sigma)
  X_true = Z_true
  p1 = p/3
  p2 = p1*2
  # 5 continuous
  X_true[,1:p1] = qexp(pnorm(Z_true[,1:p1]), rate = 1/3)
  # 5 binary
  for (m in (p1+1):p2) X_true[,m] = continuous2ordinal(X_true[,m],k=2)
  # 5 ordinal with 5 levels
  for (m in (p2+1):p) X_true[,m] = continuous2ordinal(X_true[,m],k=5)
  list(X=X_true, Sigma=Sigma)
}

trunc.rating = function(x, xmin = 1, xmax = 5){
  x = round(x)
  x[which(x<xmin)] = xmin
  x[which(x>xmax)] = xmax
  x
}

levels.percol = function(X){
  apply(X,2,function(x){length(unique(x[!is.na(x)]))})
}

map_ordinal_FAMD = function(x){
  x = as.character(x)
  l = length(x)
  y = unlist(strsplit(x, split = '_', fixed = TRUE))
  as.numeric(y[(1:l)*2])
}

## cross-validation for imputing missing values
imputeCV = function(xdata, rank, func, nfold = 5, levels=20, seed=1, verbose=FALSE){
  xdata = as.matrix(xdata)
  # generate location
  ind = which(!is.na(xdata))
  size = round(length(ind)/nfold)
  set.seed(seed)
  ind = sample(ind, length(ind))
  loc.list = list()
  for (i in 1:nfold){
    start = (i-1)*size + 1
    if (i<nfold) end = i*size else end = length(ind)
    loc.list[[i]] = ind[start:end]
  }
  
  error_cv = numeric(nfold)
  for (i in 1:nfold){
    # generate training and testing data
    Xall = xdata
    Xobs = Xall
    Xobs[loc.list[[i]]] = NA
    ind.emptyrow = which(apply(Xobs, 1, function(x){sum(!is.na(x))}) == 0)
    if (length(ind.emptyrow)>0){
      Xall = xdata[-ind.emptyrow,]
      Xobs = Xobs[-ind.emptyrow,]
    }
    ind.disc = which(apply(Xobs, 2, function(x){length(unique(x))} ) < levels)
    
    # apply algorithm
    if (func == 'xpca'){
      require(xpcaR)
      est = try(xpca(Xobs, rank = rank))
      if (class(est) == "try-error"){
        stop(paste("Use smaller rank: xPCA fails for rank ",rank, sep = ""))
      }
      xhat = est$fittedEsts
      for (s in 1:length(ind.disc)){
        xhat[,ind.disc[s]] = trunc.rating(xhat[,ind.disc[s]], 
                                          xmin = max(Xobs[,ind.disc[s]], na.rm = TRUE), 
                                          xmax = min(Xobs[,ind.disc[s]], na.rm = TRUE))
      } 
      error_cv[i] = mean(cal_mae_scaled(xhat = xhat, xobs = Xobs, x = Xall))
    }
    
    if (func == 'imputeFAMD'){
      require(missMDA)
      Xobs = as.data.frame(Xobs)
      for (s in 1:length(ind.disc)) Xobs[,ind.disc[s]] = as.factor(Xobs[,ind.disc[s]])
      est = try(imputeFAMD(Xobs, ncp = rank))
      #est = imputeFAMD(Xobs, ncp = rank)
      if (class(est) == "try-error"){
        stop(paste("Use smaller rank: imputeFAMD fails for rank ",rank, sep = ""))
      }
      xhat = est$completeObs
      for (s in 1:length(ind.disc)) xhat[,ind.disc[s]] = map_ordinal_FAMD(xhat[,ind.disc[s]])
      error_cv[i] = mean(cal_mae_scaled(xhat = xhat, xobs = Xobs, x = Xall))
    }
    
    if (verbose) print(paste("finish fold ",i, sep = ""))
  }
  error_cv
}

