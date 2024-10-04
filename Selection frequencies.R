library(MASS)
library(matrixcalc)
library(systemfit)
ny = 3 #number of dependent variables
In = diag(ny)
set.seed(123)
B1 = matrix(0.1, nrow = ny, ncol = ny)
diag(B1) = runif(ny, min = 0.5, max = 0.8) ###rho
B2 = -B1
B3 = -0.7*In
T =  400#sample size
c1 = 2 ##threshold of the first regime
c2 = 4 ##second threshold when m = 3 
gamma1 = 2 ##gamma of the first regime
gamma2 = 1 ##second gamma when m = 3
Nsim = 1000 #number of simulations
source('LMtest.R')
source('misc.R')
####Table 3 - Simulation from a VLSTAR with m = 2####
pvalLMm1 = rep(NA, Nsim)
pvalLMadm1 = rep(NA, Nsim)
pvalWilkm1 = rep(NA, Nsim)
pvalRaom1 = rep(NA, Nsim)
pvalLMm2 = rep(NA, Nsim)
pvalLMadm2 = rep(NA, Nsim)
pvalWilkm2 = rep(NA, Nsim)
pvalRaom2 = rep(NA, Nsim)
pvalLMm3 = rep(NA, Nsim)
pvalLMadm3 = rep(NA, Nsim)
pvalWilkm3 = rep(NA, Nsim)
pvalRaom3 = rep(NA, Nsim)
i=1
for(i in 1:Nsim){
  set.seed(1234+i)
  simulation <- VLSTAR.sim(n = ny, nsim = T1, B = cbind(B1, B2), c = c1, gamma = gamma1, m = 2,
                           rho = as.matrix(c(B1[1,1], B2[1,1])))
  y = simulation$sim
  plot.ts(cbind(y, simulation$st))
  #Test H0: m=1
  modvar = vars::VAR(y, p =1, type = 'none')
  BVAR = matrix(0, ncol = ny, nrow = ny)
  for(j in 1:ny){
    BVAR[,j] = coef(modvar)[[j]][,1]
  }
  object1 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = y[2:T1,], m = 1, B = BVAR,
                 residuals = rbind(rep(0, ny),residuals(modvar)))
  LMm1 = LMTEST(object1)
  LMadjm1 = FTEST(LMm1, n = ny, m = 1, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm1 = wilks(object1)
  Raom1 = RAO(object1)
  pvalLMm1[i] = LMm1$pvalue
  pvalLMadm1[i] = LMadjm1$pvalue
  pvalWilkm1[i] = Wilkm1$pvalue
  pvalRaom1[i] = Raom1$pval
  #Test H0: m=2
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NLS', n.iter = 1, constant = F,
               starting = list(cbind(rep(1.5, ny), rep(1.5, ny))), ncores = 6)
  object2 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = y[2:T1,], m = 2, gamma = mod$Gammac[,1],
                 c = mod$Gammac[,2], residuals = rbind(rep(0, ny),residuals(mod)), Gtilde = mod$Gtilde, B = t(mod$Bhat))
  LMm2 = LMTEST(object2)
  LMadjm2 = FTEST(LMm2, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm2 = wilks(object2)
  Raom2 = RAO(object2)
  pvalLMm2[i] = LMm2$pvalue
  pvalLMadm2[i] = LMadjm2$pvalue
  pvalWilkm2[i] = Wilkm2$pvalue
  pvalRaom2[i] = Raom2$pval
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
selfreq(LMtest, 0.10)
selfreq(LMtest, 0.05)
selfreq(LMtest, 0.01)

LMadtest = data.frame(LMm1 = pvalLMadm1, LMm2 = pvalLMadm2)
selfreq(LMadtest, 0.10)
selfreq(LMadtest, 0.05)
selfreq(LMadtest, 0.01)

Wilks = data.frame(LMm1 = pvalWilkm1, LMm2 = pvalWilkm2)
selfreq(Wilks, 0.10)
selfreq(Wilks, 0.05)
selfreq(Wilks, 0.01)

Rao = data.frame(LMm1 = pvalRaom1, LMm2 = pvalRaom2)
selfreq(Rao, 0.10)
selfreq(Rao, 0.05)
selfreq(Rao, 0.01)


####Table 6 - Simulation from a VTAR with m = 2####
ny = 3
In = diag(ny)
set.seed(123)
B1 = matrix(0.1, nrow = ny, ncol = ny)
diag(B1) = runif(ny, min = 0.3, max = 0.5)
B2 = -B1
B3 = 0.3*In
T1 = 400
c1 = 2
c2 = 4
gamma1 = 2
gamma2 = 1
Nsim = 1000
source('LMtest.R')
source('misc.R')
#True m = 2
pvalLMm1 = rep(NA, Nsim)
pvalLMadm1 = rep(NA, Nsim)
pvalWilkm1 = rep(NA, Nsim)
pvalRaom1 = rep(NA, Nsim)
pvalLMm2 = rep(NA, Nsim)
pvalLMadm2 = rep(NA, Nsim)
pvalWilkm2 = rep(NA, Nsim)
pvalRaom2 = rep(NA, Nsim)
pvalLMm3 = rep(NA, Nsim)
pvalLMadm3 = rep(NA, Nsim)
pvalWilkm3 = rep(NA, Nsim)
pvalRaom3 = rep(NA, Nsim)
i=1
for(i in 1:Nsim){
  set.seed(1234+i)
  simulation <- TVAR.sim(n = ny, nsim = T1, phi = list(B1, B2), c = c1, m = 2)
  y = simulation$sim
  #Test H0: m=1
  modvar = vars::VAR(y, p =1, type = 'none')
  BVAR = matrix(0, ncol = ny, nrow = ny)
  for(j in 1:ny){
    BVAR[,j] = coef(modvar)[[j]][,1]
  }
  object1 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = y[2:T1,], m = 1, B = BVAR,
                 residuals = rbind(rep(0, ny),residuals(modvar)))
  LMm1 = LMTEST(object1)
  LMadjm1 = FTEST(LMm1, n = ny, m = 1, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm1 = wilks(object1)
  Raom1 = RAO(object1)
  pvalLMm1[i] = LMm1$pvalue
  pvalLMadm1[i] = LMadjm1$pvalue
  pvalWilkm1[i] = Wilkm1$pvalue
  pvalRaom1[i] = Raom1$pval
  #Test H0: m=2
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NLS', n.iter = 1, constant = F,
               starting = list(cbind(rep(1.5, ny), rep(1.5, ny))), ncores = 6)
  object2 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = y[2:T1,], m = 2, gamma = mod$Gammac[,1],
                 c = mod$Gammac[,2], residuals = rbind(rep(0, ny),residuals(mod)), Gtilde = mod$Gtilde, B = t(mod$Bhat))
  
  LMm2 = LMTEST(object2)
  LMadjm2 = FTEST(LMm2, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
  Wilkm2 = wilks(object2)
  Raom2 = RAO(object2)
  pvalLMm2[i] = LMm2$pvalue
  pvalLMadm2[i] = LMadjm2$pvalue
  pvalWilkm2[i] = Wilkm2$pvalue
  pvalRaom2[i] = Raom2$pval
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
selfreq(LMtest, 0.10)
selfreq(LMtest, 0.05)
selfreq(LMtest, 0.01)

LMadtest = data.frame(LMm1 = pvalLMadm1, LMm2 = pvalLMadm2)
selfreq(LMadtest, 0.10)
selfreq(LMadtest, 0.05)
selfreq(LMadtest, 0.01)

Wilks = data.frame(LMm1 = pvalWilkm1, LMm2 = pvalWilkm2)
selfreq(Wilks, 0.10)
selfreq(Wilks, 0.05)
selfreq(Wilks, 0.01)

Rao = data.frame(LMm1 = pvalRaom1, LMm2 = pvalRaom2)
selfreq(Rao, 0.10)
selfreq(Rao, 0.05)
selfreq(Rao, 0.01)
