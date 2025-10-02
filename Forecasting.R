library(data.table)
library(doParallel)
library(optimParallel)
library(MASS)
library(matrixcalc)
library(systemfit)
library(sstvars)
library(vars)
source('LMtest.R')
source('misc.R')
ny = 3 #number of dependent variables
In = diag(ny)
set.seed(123)
minB = 0.7
B1 = matrix(0.1, ncol = ny, nrow = ny)
diag(B1) = minB
B2 = 0.1*diag(ny)-B1
B3 = -0.7*In
S <- array(c((B1+B2), (B1+B3), B1), dim = c(ny, ny, 3))
bound_jsr_G(S, epsilon = 0.01, adaptive_eps = TRUE, ncores = 2, print_progress = TRUE)
B1 = cbind(1, B1)
B2 = cbind(-2, B2)
B3 = cbind(2, B3)
T1 =  1000#sample size
c1 = -2 ##threshold of the first regime
c2 = 0.5 ##second threshold when m = 3 
m = 3
ytot = tsDyn::TVAR.sim(cbind(B1, B2, B3), Thresh = c(c1, c2), nthresh = 2, n = T1, lag = 1, include = 'const', mTh = 1)
ntest = 200 #testing sample dimension
window = T1-ntest
test = tail(ytot, ntest)
predVAR = matrix(NA, nrow = ntest, ncol = ny)
predVTAR2 = matrix(NA, nrow = ntest, ncol = ny)
predVTARST = matrix(NA, nrow = ntest, ncol = ny)
i = 1
for(i in 1:ntest){
  y = ytot[i:(window+i-1),]
  modvar = vars::VAR(y, p =1, type = 'const')
  BVAR = matrix(0, ncol = ny, nrow = (ny+1))
  for(j in 1:ny){
    BVAR[,j] = coef(modvar)[[j]][,1]
  }
  object1 = list(st = y[1:(window-1),1], y = y[2:window,], Data = cbind(y[1:(window-1),],1), m = 1, B = BVAR,
                 residuals = rbind(rep(0, ny),residuals(modvar)))
  LMm1 = LMTEST(object1)
  LMadjm1 = FTEST(LMm1, n = ny, m = 1, nX = ncol(object1$Data), iT = window-1) 
  if(LMadjm1$pvalue >=0.05){
    predVTAR[i,] = rbindlist(lapply(predict(modvar, n.ahead = 1)[[1]], as.data.table))$fcst
  }else if(LMadjm1$pvalue <0.05){
    modtv = tsDyn::TVAR(y, lag = 1, include = 'const', model = 'TAR')
    predVTAR2[i,] = predict(modtv, n.ahead = 1)
    tlist <- lapply(coef(modtv), t)
    Bhat = do.call(cbind, tlist)
    object3 = list(st = y[1:(window-1),1], y = y[2:window,], Data = cbind(1,y[1:(window-1),]), m = 1, B = Bhat,
                   residuals = resid(modtv))
    LMm2 = LMTEST(object3, method = 'TVAR')
    LMadjm2 = FTEST(LMm2, n = ny, m = 1, nX = ncol(object3$Data), iT = window-1) 
    if(LMadjm2$pvalue < 0.05){
      modtv = tsDyn::TVAR(y, lag = 1, include = 'const', model = 'TAR', nthresh = 2)
    }
    predVTARST[i,] = predict(modtv, n.ahead = 1)
  }
  predVAR[i,] = rbindlist(lapply(predict(modvar, n.ahead = 1)[[1]], as.data.table))$fcst
  
  print(i)
}


####Computing losses####
MSE1 = matrix(NA, ncol = 3, nrow = ntest)
MSE2 = matrix(NA, ncol = 3, nrow = ntest)
MSE3 = matrix(NA, ncol = 3, nrow = ntest)

MSE1[,1] = (test[,1]-predVAR[,1])^2
MSE1[,2] = (test[,1]-predVTAR2[,1])^2
MSE1[,3] = (test[,1]-predVTARST[,1])^2
colMeans(MSE1)

MSE2[,1] = (test[,2]-predVAR[,2])^2
MSE2[,2] = (test[,2]-predVTAR2[,2])^2
MSE2[,3] = (test[,2]-predVTARST[,2])^2
colMeans(MSE2)

MSE3[,1] = (test[,3]-predVAR[,3])^2
MSE3[,2] = (test[,3]-predVTAR2[,3])^2
MSE3[,3] = (test[,3]-predVTARST[,3])^2
colMeans(MSE3)

MSEs = cbind(colMeans(MSE1), colMeans(MSE2), colMeans(MSE3))
rowMeans(MSEs)

lossE = matrix(NA, ncol = 3, nrow = ntest)

for (i in 1:ntest){
  lossE[i,1] = euc.dist(as.matrix(test[i,]),as.matrix(predVTAR2[i,]))
  lossE[i,2] = euc.dist(as.matrix(test[i,]),as.matrix(predVTARST[i,]))
  lossE[i,3] = euc.dist(as.matrix(test[i,]),as.matrix(predVAR[i,]))
  print(i)
}
colMeans(lossE)

####Computing the GW test####
GWSE = matrix(ncol = 2, nrow = 2)
for(i in 1:2){
  GWSE[i,1] = multiGW(lossE[,i], lossE[,3], tau = 1, alternative = 'two.sided',T = nrow(lossE),3)$statistic
  GWSE[i,2] = multiGW(lossE[,i], lossE[,3], tau = 1, alternative = 'two.sided',T = nrow(lossE),3)$p.value
}
rownames(GWSE) = c('VTAR2', 'VTARST')
colnames(GWSE) = c('Statistic', 'pvalue')
