library(mixedgcImp) 
library(missForest) # required for implementing missForest
library(missMDA) # required for implementing imputeFAMD
library(xpcaR) # required for implementing xpcaR

# change the working directory to "Missing-Value-Imputation-for-Mixed-Data-via-Gaussian-Copula"
source('fun_auxi.R')
source('fun_eval.R')
load('dataset/ESL.RData')
load('dataset/LEV.RData')
load('dataset/GBSG.RData')
load('dataset/TIPS.RData')


resLEV = evalOrdinal(LEV, ratio = 0.3, rep = 100, verbose = TRUE, rankxPCA=1, rankFAMD=5)
resESL = evalOrdinal(ESL, ratio = 0.3, rep = 100, verbose = TRUE, rankxPCA=1, rankFAMD=5)
resGBSG = evalMixed(GBSG, ratio = 0.3, rep = 100, verbose = TRUE, rankxPCA=2, rankFAMD=2)
resTIPS = evalMixed(TIPS, ratio = 0.3, rep = 100, verbose = TRUE, rankxPCA=2, rankFAMD=6)
apply(resLEV, c(2,3), mean)
apply(resESL, c(2,3), mean)
apply(resGBSG, c(2,3), mean)
apply(resTIPS, c(2,3), mean)
