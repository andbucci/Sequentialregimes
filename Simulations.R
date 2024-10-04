####Simulate from VLSTAR####
library(MASS)
library(matrixcalc)
library(systemfit)
library(tsDyn)

ny = 3 #number of dependent variables
In = diag(ny)
set.seed(123)
B1 = matrix(0.1, nrow = ny, ncol = ny)
diag(B1) = runif(ny, min = 0.5, max = 0.8) ###rho
B2 = -B1
B3 = -0.7*In
T1 =  400#sample size
c1 = 2 ##threshold of the first regime
c2 = 4 ##second threshold when m = 3 
gamma1 = 2 ##gamma of the first regime
gamma2 = 1 ##second gamma when m = 3
Nsim = 1000 #number of simulations
source('LMtest.R')
source('misc.R')

#Table 1 - Empirical size simulation from m = 2
emp10 = emp05 = emp01 = rep(NA, Nsim)

emp10ad = emp05ad = emp01ad = rep(NA, Nsim)

emp10w = emp05w = emp01w = rep(NA, Nsim)

LMlist = rep(NA, Nsim)
LMlistad = rep(NA, Nsim)
Wilklist = rep(NA, Nsim)

for(m in 1:Nsim){
  set.seed(1234+m)
  simulation <- VLSTAR.sim(n = ny, nsim = T1, B = cbind(B1, B2), c = c1, gamma = gamma1, m = 2)
  y = simulation$sim
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NLS', n.iter = 1, 
               starting = list(cbind(rep(gamma1, ny), rep(c1, ny))), ncores = 6, constant = F)
  object = list(st = simulation$st[2:T1], y = y[2:T1,], Data = y[2:T1,], m = 2,
                gamma = mod$Gammac[,1], c = mod$Gammac[,2], residuals = resid(mod),
                Gtilde = mod$Gtilde, B = t(mod$Bhat))
  plot.ts(y,main="Time Series", panel=my.ts.panel)
  LM1 = LMTEST(object)
  LMadj1 = FTEST(LM1, n = ny, m = 2, nX = ncol(object$Data), iT = T1-1) 
  Wilk = wilks(object)
  LMlist[m] = LM1$LM
  LMlistad[m] = LMadj1$LM
  Wilklist[m] = Wilk$Wilks
  
  emp10[m] = LM1$pvalue <= 0.10
  emp05[m] = LM1$pvalue <= 0.05
  emp01[m] = LM1$pvalue <= 0.01
  
  emp10ad[m] = LMadj1$pvalue <= 0.10
  emp05ad[m] = LMadj1$pvalue <= 0.05
  emp01ad[m] = LMadj1$pvalue <= 0.01
  
  emp10w[m] = Wilk$pvalue <= 0.10
  emp05w[m] = Wilk$pvalue <= 0.05
  emp01w[m] = Wilk$pvalue <= 0.01
  
  print(m)
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

####Table 2 - Empirical power with m = 3####
emp10a = emp05a = emp01a = rep(NA, Nsim)

emp10b = emp05b = emp01b = rep(NA, Nsim)

emp10ad_a = emp05ad_a = emp01ad_a = rep(NA, Nsim)

emp10ad_b = emp05ad_b = emp01ad_b = rep(NA, Nsim)

emp10w = emp05w = emp01w = rep(NA, Nsim)


LMlist_a = rep(NA, Nsim)
LMlist_b = rep(NA, Nsim)
LMlistad_a = rep(NA, Nsim)
LMlistad_b = rep(NA, Nsim)

Wilklist_pot  = rep(NA, Nsim)

for(m in 1:Nsim){
  set.seed(1234+m)
  simulation <- VLSTAR.sim(n = ny, nsim = T1, B = cbind(B1, B2, B3), c = c(c1, c2), 
                           gamma = c(gamma1, gamma2), m = 3)
  y = simulation$sim
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NLS', n.iter = 1, constant = F,
                starting = list(cbind(rep(1.5, ny), rep(1.5, ny))), ncores = 6)
  object1 = list(st = simulation$st[2:T1], y = y[2:T1,], Data = y[2:T1,], m = 2, gamma = mod$Gammac[,1],
                 c = mod$Gammac[,2], residuals = rbind(rep(0, ny),residuals(mod)),
                 Gtilde = mod$Gtilde, B = t(mod$Bhat))
  LM1a = LMTEST(object1)
  LMadj1a = FTEST(LM1a, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
  Wilk = wilks(object1)
  LMlist_a[m] = LM1a$LM
  LMlist_b[m] = LM1b$LM
  LMlistad_a[m] = LMadj1a$LM
  LMlistad_b[m] = LMadj1b$LM
  Wilklist_pot[m] = Wilk$Wilks
  
  emp10a[m] = LM1a$pvalue <= 0.10
  emp05a[m] = LM1a$pvalue <= 0.05
  emp01a[m] = LM1a$pvalue <= 0.01
  
  emp10ad_a[m] = LMadj1a$pvalue <= 0.10
  emp05ad_a[m] = LMadj1a$pvalue <= 0.05
  emp01ad_a[m] = LMadj1a$pvalue <= 0.01
  
  emp10w[m] = Wilk$pvalue <= 0.10
  emp05w[m] = Wilk$pvalue <= 0.05
  emp01w[m] = Wilk$pvalue <= 0.01
  
  print(m)
}

Pot10 = sum(emp10a)*100/Nsim
Pot05 = sum(emp05a)*100/Nsim
Pot01 = sum(emp01a)*100/Nsim

Pot10ad = sum(emp10ad_a)*100/Nsim
Pot05ad = sum(emp05ad_a)*100/Nsim
Pot01ad = sum(emp01ad_a)*100/Nsim

Pot10w = sum(emp10w)*100/Nsim
Pot05w = sum(emp05w)*100/Nsim
Pot01w = sum(emp01w)*100/Nsim

print(c(Pot05, Pot05ad, Pot05w))



####Simulate from TVAR####
library(MASS)
library(matrixcalc)
library(systemfit)
library(starvars)
ny = 3
In = diag(ny)
B1 = 0.4*In
B2 = -0.4*In
B3 = -0.7*In
T = 100
c1 = 2
c2 = 4
gamma1 = 2
gamma2 = 1
Nsim = 1000
source('LMtest.R')
source('misc.R')

####Table 4 - Simulation from a VTAR with m=2####
emp10 = rep(NA, Nsim)
emp05 = rep(NA, Nsim)
emp01 = rep(NA, Nsim)

emp10ad = rep(NA, Nsim)
emp05ad = rep(NA, Nsim)
emp01ad = rep(NA, Nsim)

emp10w = rep(NA, Nsim)
emp05w = rep(NA, Nsim)
emp01w = rep(NA, Nsim)

LMlist = rep(NA, Nsim)
LMlistad = rep(NA, Nsim)
Wilklist = rep(NA, Nsim)

for(m in 1:Nsim){
  set.seed(1234+m)
  simulation <- TVAR.sim(n = ny, nsim = T, phi = list(B1, B2), c = c1, m = 2)
  y = simulation$sim
  mod = VLSTAR(y, p = 1, st = simulation$st, method = 'NL', n.iter = 1, 
              starting = list(cbind(rep(gamma1, ny), rep(c1, ny))))
  object = list(st = simulation$st[2:T], y = y[2:T,], Data = y[2:T,], m = 2, gamma = rep(gamma1, ny),
                c = rep(c1,ny), residuals = rbind(rep(0, ny), mod$residuals), Gtilde = mod$Gtilde, B = cbind(B1, B2))
  plot.ts(y,main="Time Series", panel=my.ts.panel)
  LM1 = LMTEST(object)
  LMadj1 = FTEST(LM1, n = ny, m = 2, nX = ncol(object$Data), iT = 199) 
  Wilk = wilks(object, method = 'LVSTAR')
  LMlist[m] = LM1$LM
  LMlistad[m] = LMadj1$LM
  Wilklist[m] = Wilk$Wilks
  emp10[m] = LM1$pvalue <= 0.10
  emp05[m] = LM1$pvalue <= 0.05
  emp01[m] = LM1$pvalue <= 0.01
  
  emp10ad[m] = LMadj1$pvalue <= 0.10
  emp05ad[m] = LMadj1$pvalue <= 0.05
  emp01ad[m] = LMadj1$pvalue <= 0.01
  
  emp10w[m] = Wilk$pvalue <= 0.10
  emp05w[m] = Wilk$pvalue <= 0.05
  emp01w[m] = Wilk$pvalue <= 0.01
  print(m)
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

emptest = c(empsize10,empsize05,empsize01,empsize10ad,empsize05ad,empsize01ad,empsize10w,empsize05w,empsize01w)
emptest
