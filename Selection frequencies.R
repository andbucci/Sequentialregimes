library(MASS)
library(matrixcalc)
library(systemfit)
library(vars)
source('LMtest.R')
source('misc.R')
####Table 2 - Simulation from a VLSTAR with m = 2####
ny = 3 #number of dependent variables
set.seed(123)
minB = 0.7
B1 = matrix(0.1, ncol = ny, nrow = ny)
diag(B1) = minB
B2 = 0.2*diag(ny)-B1
jsr((B1+B2), B1)###stationarity condition, joint radius lower than 1
T1 =  400#sample size
c1 = 0 ##threshold of the first regime
gamma1 = 2 ##gamma of the first regime
Nsim = 1000 #number of simulations
m = 2
pvalLMm1 = rep(NA, Nsim)
pvalLMadm1 = rep(NA, Nsim)
pvalWilkm1 = rep(NA, Nsim)
pvalLMm2 = rep(NA, Nsim)
pvalLMadm2 = rep(NA, Nsim)
pvalWilkm2 = rep(NA, Nsim)
pvalLMm3 = rep(NA, Nsim)
pvalLMadm3 = rep(NA, Nsim)
pvalWilkm3 = rep(NA, Nsim)
filename = paste('SelFreqVLSTAR', T1, 'B',sub("\\.", "", minB), 'n', ny,'m',m, '.RData', sep = '')
i=1
for(i in 1:Nsim){
  set.seed(1234+i)
  simulation <- VLSTAR.sim(n = ny, nsim = T1, B = cbind(B1, B2), c = c1, gamma = gamma1, m = 2)
  y = simulation$sim
  time <- 1:T1
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 2)) # 
  plot(time, y[,1], type = "l", col = "black", ylim = range(y[,1]), ylab = "y1", xlab = "")
  abline(h = c1, col = 'grey', lty = 3)
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,1], pch = 2, col = "blue")
    } else if (simulation$st[j] > c1) {
      points(time[j], y[j,1], pch = 15, col = "red") 
    }
  }
  plot(time, y[,2], type = "l", col = "black", ylim = range(y[,2]), ylab = "y2", xlab = "")
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,2], pch = 2, col = "blue")
    } else if (simulation$st[j] > c1) {
      points(time[j], y[j,2], pch = 15, col = "red") 
    }
  }  
  plot(time, y[,3], type = "l", col = "black", ylim = range(y[,3]), ylab = "y3", xlab = "")
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,3], pch = 2, col = "blue")
    } else if (simulation$st[j] > c1) {
      points(time[j], y[j,3], pch = 15, col = "red") 
    }
  } 
  par(mfrow = c(1, 1))
  #Test H0: m=1
  modvar = vars::VAR(y, p =1, type = 'const')
  BVAR = matrix(0, ncol = ny, nrow = (ny+1))
  for(j in 1:ny){
    BVAR[,j] = coef(modvar)[[j]][,1]
  }
  object1 = list(st = y[1:(T1-1),1], y = y[2:T1,], Data = cbind(y[1:(T1-1),],1), m = 1, B = BVAR,
                 residuals = rbind(rep(0, ny),residuals(modvar)))
  LMm1 = LMTEST(object1)
  LMadjm1 = FTEST(LMm1, n = ny, m = 1, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm1 = wilks(object1)
  pvalLMm1[i] = LMm1$pvalue
  pvalLMadm1[i] = LMadjm1$pvalue
  pvalWilkm1[i] = Wilkm1$pvalue
  #Test H0: m=2
  stam2 = starting(y, st = simulation$st, n.combi = 20, ncores = 8)
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NLS', n.iter = 20,
               starting = stam2,
               ncores = 6, constant = T, maxgamma = 20)
  object2 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = cbind(1,y[1:(T1-1),]), m = 2,
                gamma = mod$Gammac[,1], c = mod$Gammac[,2], residuals = resid(mod),
                Gtilde = mod$Gtilde, B =mod$B)
  LMm2 = LMTEST(object2)
  LMadjm2 = FTEST(LMm2, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm2 = wilks(object2)
  pvalLMm2[i] = LMm2$pvalue
  pvalLMadm2[i] = LMadjm2$pvalue
  pvalWilkm2[i] = Wilkm2$pvalue
  #save.image(filename)
  if(pvalLMm1[i]>0.10){cat('m = 1')
  }else if(pvalLMm1[i]<0.10 & pvalLMm2[i]>0.10){
      cat('m = 2')}else if(pvalLMm1[i]<0.10 & pvalLMm2[i]<0.10){cat('m = 3')}
  print(i)
}

selfreq = function(data, alpha = 0.05){
  data$m = NA
  data$m[data$LMm1 > alpha] = 1
  data$m[data$LMm1 < alpha & data$LMm2>alpha] = 2
  data$m[data$LMm1 < alpha & data$LMm2<alpha] = 3
  return(table(data$m)/Nsim)
}

LMtest = data.frame(LMm1 = pvalLMm1, LMm2 = pvalLMm2)
selfreq(LMtest, 0.10)*100
selfreq(LMtest, 0.05)*100
selfreq(LMtest, 0.01)*100

LMadtest = data.frame(LMm1 = pvalLMadm1, LMm2 = pvalLMadm2)
selfreq(LMadtest, 0.10)*100
selfreq(LMadtest, 0.05)*100
selfreq(LMadtest, 0.01)*100

Wilks = data.frame(LMm1 = pvalWilkm1, LMm2 = pvalWilkm2)
selfreq(Wilks, 0.10)*100
selfreq(Wilks, 0.05)*100
selfreq(Wilks, 0.01)*100

save.image(filename)

####Table 3 - Simulation from a VLSTAR with m = 3####
library(sparsevar)
library(sstvars)
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
B1 = rbind(1, B1)
B2 = rbind(-2, B2)
B3 = rbind(2, B3)
T1 =  1000#sample size
c1 = -2 ##threshold of the first regime
c2 = 0.5 ##second threshold when m = 3 
gamma1 = 2 ##gamma of the first regime
gamma2 = 40 ##second gamma when m = 3
Nsim = 1000 #number of simulations
source('LMtest.R')
source('misc.R')
m = 3
pvalLMm1 = rep(NA, Nsim)
pvalLMadm1 = rep(NA, Nsim)
pvalWilkm1 = rep(NA, Nsim)
pvalLMm2 = rep(NA, Nsim)
pvalLMadm2 = rep(NA, Nsim)
pvalWilkm2 = rep(NA, Nsim)
pvalLMm3 = rep(NA, Nsim)
pvalLMadm3 = rep(NA, Nsim)
pvalWilkm3 = rep(NA, Nsim)
filename = paste('SelFreqVLSTAR', T1, 'B',sub("\\.", "", minB), 'n', ny,'m',m, '.RData', sep = '')
i=1
for(i in 1:Nsim){
  set.seed(1234+i)
  simulation <- VLSTAR.sim(n = ny, nsim = T1, B = cbind(B1, B2, B3), c = c(c1, c2), 
                           gamma = c(gamma1, gamma2), m = 3, const = T)
  y = simulation$sim
  time <- 1:T1
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 2)) # 
  plot(time, y[,1], type = "l", col = "black", ylim = range(y[,1]), ylab = "y1", xlab = "")
  abline(h = c(c1, c2), col = 'grey', lty = 3)
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,1], pch = 2, col = "blue")
    } else if (simulation$st[j] >= c1 && simulation$st[j] <= c2) {
      points(time[j],y[j,1], pch = 16, col = "green")
    } else if (simulation$st[j] > c2) {
      points(time[j], y[j,1], pch = 15, col = "red") 
    }
  }
  plot(time, y[,2], type = "l", col = "black", ylim = range(y[,2]), ylab = "y2", xlab = "")
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,2], pch = 2, col = "blue")
    } else if (simulation$st[j] >= c1 && simulation$st[j] <= c2) {
      points(time[j],y[j,2], pch = 16, col = "green")
    } else if (simulation$st[j] > c2) {
      points(time[j], y[j,2], pch = 15, col = "red") 
    }
  }
  plot(time, y[,3], type = "l", col = "black", ylim = range(y[,3]), ylab = "y3", xlab = "")
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,3], pch = 2, col = "blue")
    } else if (simulation$st[j] >= c1 && simulation$st[j] <= c2) {
      points(time[j],y[j,3], pch = 16, col = "green")
    } else if (simulation$st[j] > c2) {
      points(time[j], y[j,3], pch = 15, col = "red") 
    }
  }
  par(mfrow = c(1, 1))
  #Test H0: m=1
  modvar = vars::VAR(y, p =1, type = 'const')
  BVAR = matrix(0, ncol = ny, nrow = (ny+1))
  for(j in 1:ny){
    BVAR[,j] = coef(modvar)[[j]][,1]
  }
  object1 = list(st = y[1:(T1-1),1], y = y[2:T1,], Data = cbind(y[1:(T1-1),],1), m = 1, B = BVAR,
                 residuals = rbind(rep(0, ny),residuals(modvar)))
  LMm1 = LMTEST(object1)
  LMadjm1 = FTEST(LMm1, n = ny, m = 1, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm1 = wilks(object1)
  pvalLMm1[i] = LMm1$pvalue
  pvalLMadm1[i] = LMadjm1$pvalue
  pvalWilkm1[i] = Wilkm1$pvalue
  #Test H0: m=2
  stam2 = starting(y, st = simulation$st, n.combi = 20, ncores = 8)
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NLS', n.iter = 20,
               starting = stam2,
               ncores = 6, constant = T, maxgamma = 50)
  object2 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = cbind(1,y[1:(T1-1),]), m = 2,
                 gamma = mod$Gammac[,1], c = mod$Gammac[,2], residuals = resid(mod),
                 Gtilde = mod$Gtilde, B =mod$B)
  LMm2 = LMTEST(object2)
  print(LMm2$pvalue)
  LMadjm2 = FTEST(LMm2, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm2 = wilks(object2)
  pvalLMm2[i] = LMm2$pvalue
  pvalLMadm2[i] = LMadjm2$pvalue
  pvalWilkm2[i] = Wilkm2$pvalue
  save.image(filename)
  if(pvalLMm1[i]>0.10){cat('m = 1 \n')
  }else if(pvalLMm1[i]<0.10 & pvalLMm2[i]>0.10){
    cat('m = 2 \n')}else if(pvalLMm1[i]<0.10 & pvalLMm2[i]<0.10){cat('m = 3 \n')}
  print(i)
}

selfreq = function(data, alpha = 0.05){
  data$m = NA
  data$m[data$LMm1 > alpha] = 1
  data$m[data$LMm1 < alpha & data$LMm2>alpha] = 2
  data$m[data$LMm1 < alpha & data$LMm2<alpha] = 3
  return(table(data$m)/Nsim)
}

LMtest = data.frame(LMm1 = pvalLMm1, LMm2 = pvalLMm2)
selfreq(LMtest, 0.10)*100
selfreq(LMtest, 0.05)*100
selfreq(LMtest, 0.01)*100

LMadtest = data.frame(LMm1 = pvalLMadm1, LMm2 = pvalLMadm2)
selfreq(LMadtest, 0.10)*100
selfreq(LMadtest, 0.05)*100
selfreq(LMadtest, 0.01)*100

Wilks = data.frame(LMm1 = pvalWilkm1, LMm2 = pvalWilkm2)
selfreq(Wilks, 0.10)*100
selfreq(Wilks, 0.05)*100
selfreq(Wilks, 0.01)*100

save.image(filename)

####Table 4 - Simulation from a VTAR with m = 2####

library(doParallel)
library(MASS)
library(tsDyn)
library(zoo)
ny = 3
In = diag(ny)
set.seed(123)
minB = 0.7
B1 = matrix(0.1, ncol = ny, nrow = ny)
diag(B1) = minB
B2 = 0.2*diag(ny)-B1
T1 = 1000
c1 = 0
Nsim = 1000
m=2
source('LMtest.R')
source('misc.R')
#True m = 2
pvalLMm1 = rep(NA, Nsim)
pvalLMadm1 = rep(NA, Nsim)
pvalWilkm1 = rep(NA, Nsim)
pvalLMm2 = rep(NA, Nsim)
pvalLMadm2 = rep(NA, Nsim)
pvalWilkm2 = rep(NA, Nsim)
pvalLMm3 = rep(NA, Nsim)
pvalLMadm3 = rep(NA, Nsim)
pvalWilkm3 = rep(NA, Nsim)
i=1
filename = paste('SelFreqVTAR', T1, 'B',sub("\\.", "", minB), 'n', ny,'m',m, '.RData', sep = '')
for(i in 1:Nsim){
  set.seed(1234+i)
 # simulation <- TVAR.sim(n = ny, nsim = T1, phi = list(B1, B2), c = c1, m = 2)
  y = tsDyn::TVAR.sim(cbind(B1, B2), Thresh = c1, nthresh = 1, n = (T1+1), lag = 1, include = 'none', mTh = 1)
  simulation = list(y = y[-1,], st = y[1:T1,1])
  y = y[-1,]
  #y = simulation$sim
  time <- 1:T1
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 2)) # 
  plot(time, simulation$y[,1], type = "l", col = "black", ylim = range(y[,1]), ylab = "y1", xlab = "")
  abline(h = c1, col = 'grey', lty = 3)
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], simulation$y[j,1], pch = 2, col = "blue")
    } else if (simulation$st[j] > c1) {
      points(time[j], simulation$y[j,1], pch = 15, col = "red") 
    }
  }
  plot(time, simulation$y[,2], type = "l", col = "black", ylim = range(y[,2]), ylab = "y2", xlab = "")
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], simulation$y[j,2], pch = 2, col = "blue")
    } else if (simulation$st[j] > c1) {
      points(time[j], simulation$y[j,2], pch = 15, col = "red") 
    }
  }  
  plot(time, simulation$y[,3], type = "l", col = "black", ylim = range(y[,3]), ylab = "y3", xlab = "")
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], simulation$y[j,3], pch = 2, col = "blue")
    } else if (simulation$st[j] > c1) {
      points(time[j], simulation$y[j,3], pch = 15, col = "red") 
    }
  } 
  par(mfrow = c(1, 1))
  #Test H0: m=1
  modvar = vars::VAR(y, p =1, type = 'const')
  BVAR = matrix(0, ncol = ny, nrow = (ny+1))
  for(j in 1:ny){
    BVAR[,j] = coef(modvar)[[j]][,1]
  }
  object1 = list(st = y[1:(T1-1),1], y = y[2:T1,], Data = cbind(y[1:(T1-1),],1), m = 1, B = BVAR,
                 residuals = rbind(rep(0, ny),residuals(modvar)))
  LMm1 = LMTEST(object1)
  LMadjm1 = FTEST(LMm1, n = ny, m = 1, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm1 = wilks(object1)
  pvalLMm1[i] = LMm1$pvalue
  pvalLMadm1[i] = LMadjm1$pvalue
  pvalWilkm1[i] = Wilkm1$pvalue
  #Test H0: m=2
  stam2 = starting(y, st = simulation$st, n.combi = 20, ncores = 8)
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NLS', n.iter = 20,
               starting = stam2,
               ncores = 6, constant = T, maxgamma = 20)
  object2 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = cbind(1,y[1:(T1-1),]), m = 2,
                 gamma = mod$Gammac[,1], c = mod$Gammac[,2], residuals = resid(mod),
                 Gtilde = mod$Gtilde, B =mod$B)
  LMm2 = LMTEST(object2)
  LMadjm2 = FTEST(LMm2, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm2 = wilks(object2)
  pvalLMm2[i] = LMm2$pvalue
  pvalLMadm2[i] = LMadjm2$pvalue
  pvalWilkm2[i] = Wilkm2$pvalue
  #Test H0: m=2 ST APPROACH
  modtv = tsDyn::TVAR(y, lag = 1, include = 'const', model = 'TAR')
  tlist <- lapply(coef(modtv), t)
  Bhat = do.call(cbind, tlist)
  object3 = list(st = y[1:(T1-1),1], y = y[2:T1,], Data = cbind(1,y[1:(T1-1),]), m = 1, B = Bhat,
                 residuals = resid(modtv))
  LMm3 = LMTEST(object3, method = 'TVAR')
  print(LMm3$pval)
  LMadjm3 = FTEST(LMm3, n = ny, m = 1, nX = ncol(object3$Data), iT = T1-1) 
  Wilkm3 = wilks(object3, method = 'TVAR')
  pvalLMm3[i] = LMm3$pvalue
  pvalLMadm3[i] = LMadjm3$pvalue
  pvalWilkm3[i] = Wilkm3$pvalue
  save.image(filename)
  if(pvalLMm1[i]>0.10){cat('m = 1 \n')
  }else if(pvalLMm1[i]<0.10 & pvalLMm2[i]>0.10){
    cat('m = 2 \n')}else if(pvalLMm1[i]<0.10 & pvalLMm2[i]<0.10){cat('m = 3 \n')}
  print(i)
}

selfreq = function(data, alpha = 0.05){
  data$m = NA
  data$m[data$LMm1 > alpha] = 1
  data$m[data$LMm1 < alpha & data$LMm2>alpha] = 2
  data$m[data$LMm1 < alpha & data$LMm2<alpha] = 3
  return(table(data$m)/Nsim)
}

LMtest = data.frame(LMm1 = pvalLMm1, LMm2 = pvalLMm2)
selfreq(LMtest, 0.10)*100
selfreq(LMtest, 0.05)*100
selfreq(LMtest, 0.01)*100

LMadtest = data.frame(LMm1 = pvalLMadm1, LMm2 = pvalLMadm2)
selfreq(LMadtest, 0.10)*100
selfreq(LMadtest, 0.05)*100
selfreq(LMadtest, 0.01)*100

Wilks = data.frame(LMm1 = pvalWilkm1, LMm2 = pvalWilkm2)
selfreq(Wilks, 0.10)*100
selfreq(Wilks, 0.05)*100
selfreq(Wilks, 0.01)*100

LMtestST = data.frame(LMm1 = pvalLMm1, LMm2 = pvalLMm3)
selfreq(LMtestST, 0.10)*100
selfreq(LMtestST, 0.05)*100
selfreq(LMtestST, 0.01)*100

LMadtestST = data.frame(LMm1 = pvalLMadm1, LMm2 = pvalLMadm3)
selfreq(LMadtestST, 0.10)*100
selfreq(LMadtestST, 0.05)*100
selfreq(LMadtestST, 0.01)*100

WilksST = data.frame(LMm1 = pvalWilkm1, LMm2 = pvalWilkm3)
selfreq(WilksST, 0.10)*100
selfreq(WilksST, 0.05)*100
selfreq(WilksST, 0.01)*100

####Table 5 - Simulation from a VTAR with m = 3####
library(sparsevar)
library(sstvars)
library(doParallel)
library(MASS)
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
T1 =  400#sample size
c1 = -2 ##threshold of the first regime
c2 = 0.5 ##second threshold when m = 3 
Nsim = 1000 #number of simulations
m = 3
source('LMtest.R')
source('misc.R')
#True m = 2
pvalLMm1 = rep(NA, Nsim)
pvalLMadm1 = rep(NA, Nsim)
pvalWilkm1 = rep(NA, Nsim)
pvalLMm2 = rep(NA, Nsim)
pvalLMadm2 = rep(NA, Nsim)
pvalWilkm2 = rep(NA, Nsim)
pvalLMm3 = rep(NA, Nsim)
pvalLMadm3 = rep(NA, Nsim)
pvalWilkm3 = rep(NA, Nsim)
filename = paste('SelFreqVTAR', T1, 'B',sub("\\.", "", minB), 'n', ny,'m',m, '.RData', sep = '')
i=1
for(i in 1:Nsim){
  set.seed(1234+i)
  y = tsDyn::TVAR.sim(cbind(B1, B2, B3), Thresh = c(c1, c2), nthresh = 2, n = (T1+1), lag = 1, include = 'const', mTh = 1)
  simulation = list(y = y[-1,], st = y[1:T1,1])
  y = y[-1,]
  time <- 1:T1
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 2)) # 
  plot(time, y[,1], type = "l", col = "black", ylim = range(y[,1]), ylab = "y1", xlab = "")
  abline(h = c(c1, c2), col = 'grey', lty = 3)
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,1], pch = 2, col = "blue")
    } else if (simulation$st[j] >= c1 && simulation$st[j] <= c2) {
      points(time[j],y[j,1], pch = 16, col = "green")
    } else if (simulation$st[j] > c2) {
      points(time[j], y[j,1], pch = 15, col = "red") 
    }
  }
  plot(time, y[,2], type = "l", col = "black", ylim = range(y[,2]), ylab = "y2", xlab = "")
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,2], pch = 2, col = "blue")
    } else if (simulation$st[j] >= c1 && simulation$st[j] <= c2) {
      points(time[j],y[j,2], pch = 16, col = "green")
    } else if (simulation$st[j] > c2) {
      points(time[j], y[j,2], pch = 15, col = "red") 
    }
  }
  plot(time, y[,3], type = "l", col = "black", ylim = range(y[,3]), ylab = "y3", xlab = "")
  for (j in 1:T1) {
    if (simulation$st[j] < c1) {
      points(time[j], y[j,3], pch = 2, col = "blue")
    } else if (simulation$st[j] >= c1 && simulation$st[j] <= c2) {
      points(time[j],y[j,3], pch = 16, col = "green")
    } else if (simulation$st[j] > c2) {
      points(time[j], y[j,3], pch = 15, col = "red") 
    }
  }
  par(mfrow = c(1, 1))
  #Test H0: m=1
  modvar = vars::VAR(y, p =1, type = 'const')
  BVAR = matrix(0, ncol = ny, nrow = (ny+1))
  for(j in 1:ny){
    BVAR[,j] = coef(modvar)[[j]][,1]
  }
  object1 = list(st = y[1:(T1-1),1], y = y[2:T1,], Data = cbind(y[1:(T1-1),],1), m = 1, B = BVAR,
                 residuals = rbind(rep(0, ny),residuals(modvar)))
  LMm1 = LMTEST(object1)
  LMadjm1 = FTEST(LMm1, n = ny, m = 1, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm1 = wilks(object1)
  pvalLMm1[i] = LMm1$pvalue
  pvalLMadm1[i] = LMadjm1$pvalue
  pvalWilkm1[i] = Wilkm1$pvalue
  #Test H0: m=2
  stam2 = starting(y, st = simulation$st, n.combi = 20, ncores = 8)
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NLS', n.iter = 20,
               starting = stam2,
               ncores = 6, constant = T, maxgamma = 50)
  object2 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = cbind(1,y[1:(T1-1),]), m = 2,
                 gamma = mod$Gammac[,1], c = mod$Gammac[,2], residuals = resid(mod),
                 Gtilde = mod$Gtilde, B =mod$B)
  LMm2 = LMTEST(object2)
  print(LMm2$pval)
  LMadjm2 = FTEST(LMm2, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm2 = wilks(object2)
  pvalLMm2[i] = LMm2$pvalue
  pvalLMadm2[i] = LMadjm2$pvalue
  pvalWilkm2[i] = Wilkm2$pvalue
  #Test H0: m=2 ST APPROACH
  modtv = tsDyn::TVAR(y, lag = 1, include = 'const', model = 'TAR')
  tlist <- lapply(coef(modtv), t)
  Bhat = do.call(cbind, tlist)
  object3 = list(st = y[1:(T1-1),1], y = y[2:T1,], Data = cbind(1,y[1:(T1-1),]), m = 1, B = Bhat,
                 residuals = resid(modtv))
  LMm3 = LMTEST(object3, method = 'TVAR')
  print(LMm3$pval)
  LMadjm3 = FTEST(LMm3, n = ny, m = 1, nX = ncol(object3$Data), iT = T1-1) 
  Wilkm3 = wilks(object3, method = 'TVAR')
  pvalLMm3[i] = LMm3$pvalue
  pvalLMadm3[i] = LMadjm3$pvalue
  pvalWilkm3[i] = Wilkm3$pvalue
  save.image(filename)
  if(pvalLMm1[i]>0.10){cat('m = 1 \n')
  }else if(pvalLMm1[i]<0.10 & pvalLMm2[i]>0.10){
    cat('m = 2 \n')}else if(pvalLMm1[i]<0.10 & pvalLMm2[i]<0.10){cat('m = 3 \n')}
  print(i)
}

selfreq = function(data, alpha = 0.05){
  data$m = NA
  data$m[data$LMm1 > alpha] = 1
  data$m[data$LMm1 < alpha & data$LMm2>alpha] = 2
  data$m[data$LMm1 < alpha & data$LMm2<alpha] = 3
  return(table(data$m)/Nsim)
}

LMtest = data.frame(LMm1 = pvalLMm1, LMm2 = pvalLMm2)
selfreq(LMtest, 0.10)*100
selfreq(LMtest, 0.05)*100
selfreq(LMtest, 0.01)*100

LMadtest = data.frame(LMm1 = pvalLMadm1, LMm2 = pvalLMadm2)
selfreq(LMadtest, 0.10)*100
selfreq(LMadtest, 0.05)*100
selfreq(LMadtest, 0.01)*100

Wilks = data.frame(LMm1 = pvalWilkm1, LMm2 = pvalWilkm2)
selfreq(Wilks, 0.10)*100
selfreq(Wilks, 0.05)*100
selfreq(Wilks, 0.01)*100

LMtestST = data.frame(LMm1 = pvalLMm1, LMm2 = pvalLMm3)
selfreq(LMtestST, 0.10)*100
selfreq(LMtestST, 0.05)*100
selfreq(LMtestST, 0.01)*100

LMadtestST = data.frame(LMm1 = pvalLMadm1, LMm2 = pvalLMadm3)
selfreq(LMadtestST, 0.10)*100
selfreq(LMadtestST, 0.05)*100
selfreq(LMadtestST, 0.01)*100

WilksST = data.frame(LMm1 = pvalWilkm1, LMm2 = pvalWilkm3)
selfreq(WilksST, 0.10)*100
selfreq(WilksST, 0.05)*100
selfreq(WilksST, 0.01)*100
