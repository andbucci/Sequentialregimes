####Simulate from VLSTAR####
library(MASS)
library(matrixcalc)
library(systemfit)
library(tsDyn)
library(doParallel)
library(optimParallel)

source('LMtest.R')
source('misc.R')

ny = 5 #number of dependent variables
m = 1
set.seed(123)
minB = 0.6
B1 = matrix(0.1, ncol = ny, nrow = ny)
diag(B1) = minB
T1 =  600#sample size
Nsim = 1000 #number of simulations

#Table 1 - Empirical size simulation from m = 2
emp10 = emp05 = emp01 = rep(NA, Nsim)
emp10ad = emp05ad = emp01ad = rep(NA, Nsim)
emp10w = emp05w = emp01w = rep(NA, Nsim)

LMlist = rep(NA, Nsim)
LMlistad = rep(NA, Nsim)
Wilklist = rep(NA, Nsim)
gammamat = matrix(ncol = ny, nrow = Nsim)
cmat = matrix(ncol = ny, nrow = Nsim)
i=1
filename = paste('Empsize', T1, 'B',sub("\\.", "", minB), 'n', ny,'m',m, '.RData', sep = '')
for(i in 1:Nsim){
  set.seed(1234+i)
  y = VAR.sim(B = B1, lag= 1, include = 'none', n = (T1+1))
  modvar = vars::VAR(y, p =1, type = 'none')
  BVAR = matrix(0, ncol = ny, nrow = ny)
  for(j in 1:ny){
    BVAR[,j] = coef(modvar)[[j]][,1]
  }
  object = list(st = y[1:(T1-1),1], y = y[2:T1,], Data = y[1:(T1-1),], m = 1, B = BVAR,
                 residuals = rbind(rep(0, ny),residuals(modvar)))
  plot.ts(y,main="Time Series", panel=my.ts.panel)
  LM1 = LMTEST(object)
  LMadj1 = FTEST(LM1, n = ny, m = 1, nX = ncol(object$Data), iT = T1-2) 
  Wilk = wilks(object)
  LMlist[i] = LM1$LM
  LMlistad[i] = LMadj1$LM
  Wilklist[i] = Wilk$Wilks
  
  emp10[i] = LM1$pvalue <= 0.10
  emp05[i] = LM1$pvalue <= 0.05
  emp01[i] = LM1$pvalue <= 0.01
  
  emp10ad[i] = LMadj1$pvalue <= 0.10
  emp05ad[i] = LMadj1$pvalue <= 0.05
  emp01ad[i] = LMadj1$pvalue <= 0.01
  
  emp10w[i] = Wilk$pvalue <= 0.10
  emp05w[i] = Wilk$pvalue <= 0.05
  emp01w[i] = Wilk$pvalue <= 0.01
  #save.image(filename)
  print(i)
  print(paste('Partial size: ', sum(emp10[1:i])*100/i, ', p-val:', round(LM1$pvalue,3), sep = ''))
}
empsize10 = sum(emp10)*100/Nsim
empsize05 = sum(emp05)*100/Nsim
empsize01 = sum(emp01)*100/Nsim

empsize10ad = sum(emp10ad)*100/Nsim
empsize05ad = sum(emp05ad)*100/Nsim
empsize01ad = sum(emp01ad)*100/Nsim

empsize10w = sum(emp10w)*100/Nsim
empsize05w = sum(emp05w)*100/Nsim
empsize01w = sum(emp01w)*100/Nsim

emptest = c(empsize10,empsize05,empsize01,empsize10ad,empsize05ad,empsize01ad,
            empsize10w,empsize05w,empsize01w)
emptest
save.image(filename)
