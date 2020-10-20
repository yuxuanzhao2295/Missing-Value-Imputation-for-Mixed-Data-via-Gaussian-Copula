# Examples:
# set up working directory in your place
setwd("/Users/yuxuan/Documents/GitHub/Missing-Value-Imputation-for-Mixed-Data-via-Gaussian-Copula")
source('fun_auxi.R')
load("dataset/GBSG.RData")
load("dataset/cal500exp.RData")
# GBSG
gbsg_cv_xpca = array(0, dim = c(4,5)) 
for (i in 1:4) gbsg_cv_xpca[i,] = imputeCV(GBSG, rank = i, func = 'xpca') 
apply(gbsg_cv_xpca, 1, mean) # rank 2 returns best performance 

gbsg_cv_FAMD = array(0, dim = c(4,5)) 
for (i in 1:4) gbsg_cv_FAMD[i,] = imputeCV(GBSG, rank = i, func = 'imputeFAMD')
apply(gbsg_cv_FAMD, 1, mean) # rank 2 returns best performance 


# cal500exp
cal500_cv_xpca = array(0, dim = c(5,5))
for (i in 1:5){
  cal500_cv_xpca[i,] = imputeCV(cal500exp, rank = i, func = 'xpca')
  print(paste('finished ', i))
}
apply(cal500_cv_xpca, 1, mean) # rank 1 returns best performance 

# important to remove colnames with "_"
data = cal500exp
colnames(data) = NULL
cal500_cv_FAMD = array(0, dim = c(3,5))
for (i in 1:3){
  cal500_cv_FAMD[i,] = imputeCV(data, rank = 13+i, func = 'imputeFAMD')
  print(paste('finished ', i))
}
apply(cal500_cv_FAMD, 1, mean)  # rank 15 returns best performance 


