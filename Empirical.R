library(MASS)
library(MTS)
library(doParallel)
library(optimParallel)
library(plotly)
library(starvars)
library(tsDyn)
library(tseries)
library(urca)
library(xts)
library(zoo)
source('LMtest.R')
source('misc.R')
####Empirical application on interest rates - Section 6.1####
data <- read.csv2('TsayData.csv')
data$Date <- as.Date(data$Date, format = '%d/%m/%Y')
data[,2:7] <- apply(data[,2:7], 2, as.numeric)
data = data[3:nrow(data),]

###Figure 1
tindex <- seq(as.Date("1953/06/01", format = '%Y/%m/%d'), as.Date("2022/09/30", format = '%Y/%m/%d'), "months")

p1 <- plot_ly(x = tindex, y = data$g3m, type = 'scatter', mode = 'lines', 
              fill = "tozeroy", name = '3-month bill')#fillcolor
p2 <- plot_ly(x = tindex, y = data$g3Y, type="scatter", mode="lines",  
              fill = "tozeroy",name = '3-year treasury')
p3 <- plot_ly(x = tindex, y = data$spreadavg, type="scatter", mode="lines",  
              fill = "tozeroy",name = 'Moving average of spread')
fig1 <- subplot(p1, p2,p3, nrows = 3)
fig1

lvl <- log(data[, 2:3])
jotest = ca.jo(lvl[2:nrow(lvl),], ecdet = 'none')
summary(jotest)

####Table 10####
mY = data[2:nrow(data),4:5]
mX = data[1:(nrow(data)-1),4:5]
st = data$spreadavg[2:nrow(data)]


####LM, LMresc and Wilks
set.seed(3010)
T1 = nrow(mY)
ny = ncol(mY)
modvar = vars::VAR(mY, p =1, type = 'const')
AIC1 = AIC(modvar)
BVAR = matrix(0, ncol = ny, nrow = (ny+1))
for(j in 1:ny){
  BVAR[,j] = coef(modvar)[[j]][,1]
}
object0 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(mY[1:(nrow(mY)-1),],1), m = 1, B = BVAR,
               residuals = rbind(rep(0, ny),residuals(modvar)))
LM0 = LMTEST(object0)
LMadj0 = FTEST(LM0, n = ny, m = 1, nX = ncol(object0$Data), iT = T1-1) 
Wilk0 <- wilks(object0)
Panel0 = cbind(rbind(LM0$LM, LMadj0$LM, Wilk0$Wilks), rbind(LM0$pval, LMadj0$pval, Wilk0$pval))
Panel0

####Test vs m = 3
m = 2
stam2 = starting(mY, st = st, n.combi = 20, ncores = 8)
mod = VLSTAR(mY, p = 1, st = st, method = 'NLS', n.iter = 20,
             starting = stam2,
             ncores = 6, constant = T, maxgamma = 50)
AIC2 = sum(mod$AIC)
object1 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1,mY[1:(T1-1),]), m = 2,
               gamma = mod$Gammac[,1], c = mod$Gammac[,2], residuals = resid(mod),
               Gtilde = mod$Gtilde, B =mod$B)
LM1 = LMTEST(object1)
LMadj1 = FTEST(LM1, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
Wilk1 <- wilks(object1)
##ST approach
modtv = tsDyn::TVAR(mY, lag = 1, include = 'const', model = 'TAR', thVar = st)
AIC2ST = AIC(modtv)
tlist <- lapply(coef(modtv), t)
Bhat = do.call(cbind, tlist)
objectST = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1,mY[1:(T1-1),]), m = 1, B = Bhat,
               residuals = resid(modtv))
LM1ST = LMTEST(objectST, method = 'TVAR')
LMadj1ST = FTEST(LM1ST, n = ny, m = 1, nX = ncol(objectST$Data), iT = T1-1) 
Wilk1ST = wilks(objectST, method = 'TVAR')
Panel1 = cbind(rbind(LM1$LM, LMadj1$LM, Wilk1$Wilks), rbind(LM1$pval, LMadj1$pval, Wilk1$pval), 
               rbind(LM1ST$LM, LMadj1ST$LM, Wilk1ST$Wilks), rbind(LM1ST$pval, LMadj1ST$pval, Wilk1ST$pval))
Panel1

##Estimate m = 3 (for AIC)#
m = 3
stam3 = starting(mY, st = st, n.combi = 20, ncores = 8, m = 3)
mod2 = VLSTAR(mY, p = 1, st = st, method = 'NLS', n.iter = 100,
             starting = stam3, m = 3,
             ncores = 6, constant = T, maxgamma = 50)
modtv2 = tsDyn::TVAR(mY, lag = 1, include = 'const', model = 'TAR', thVar = st, nthresh = 2, trim = 0.05)
AIC3 = sum(mod2$AIC)
AIC3ST = AIC(modtv2)

#Table 11
c(AIC1, AIC2, AIC3)
c(AIC1, AIC2ST, AIC3ST)

##Figure with regimes identified through ST
time = data$Date[-1]
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2)) # 
plot(time, mY[,1], type = "l", col = "black", ylim = range(mY[,1]), ylab = "3-month", xlab = "")
for (j in 1:T1) {
  if (st[j] < getTh(modtv2)[1]) {
    points(time[j], mY[j,1], pch = 2, col = "blue")
  } else if (st[j] > getTh(modtv2)[1] & st[j] <= getTh(modtv2)[2]) {
    points(time[j], mY[j,1], pch = 15, col = "red") 
  } else if(st[j] > getTh(modtv2)[2]){
    points(time[j], mY[j,1], pch = 15, col = "grey")
  }
    
}
plot(time, mY[,2], type = "l", col = "black", ylim = range(mY[,2]), ylab = "3-year", xlab = "")
for (j in 1:T1) {
  if (st[j] < getTh(modtv2)[1]) {
    points(time[j], mY[j,2], pch = 2, col = "blue")
  } else if (st[j] > getTh(modtv2)[1] & st[j] <= getTh(modtv2)[2]) {
    points(time[j], mY[j,2], pch = 15, col = "red") 
  } else if(st[j] > getTh(modtv2)[2]){
    points(time[j], mY[j,2], pch = 15, col = "grey")
  }
}  
par(mfrow = c(1, 1))

##VECM
modVECM <- VECM(lvl, lag = 1, r = 1, include = "none", estim = "ML")
resVECM = residuals(modVECM)
#MSE VECM vs MSE VLSTAR/VTAR
sum(diag(resVECM %*% t(resVECM)))/T1
sum(diag(resVLSTAR3 %*% t(resVLSTAR3)))/T1
sum(diag(resVTAR3 %*% t(resVTAR3)))/T1

####Empirical application on river flows data - Section 6.2####
###Figure 2 and Table 12
data("ice.river")
source('LMtest.R')
source('misc.R')

tindex <- seq(as.Date("1972/1/1"), as.Date("1974/12/31"), "days")
p1 <- plot_ly(x = tindex, y = flow.jok, type = 'scatter', mode = 'lines', 
             fill = "tozeroy", name = 'Jökulsá')#fillcolor
p2 <- plot_ly(x = tindex, y = flow.vat, type = 'scatter', mode = 'lines', 
              fill = "tozeroy", name = 'Vatndalsá')#fillcolor
p3 <- plot_ly(x = tindex, y = temp, type = 'scatter', mode = 'lines', 
              fill = "tozeroy", name = 'Temperature')#fillcolor
p4 <- plot_ly(x = tindex, y = prec, type = 'scatter', mode = 'lines', 
              fill = "tozeroy", name = 'Precipitation')#fillcolor
fig2 <- subplot(p1, p2, p3, p4, nrows = 4)
fig2

####Table 12 - Panel A####
##st = temperature##
mY = ice.river[2:nrow(ice.river),1:2]
mX = ice.river[1:(nrow(ice.river)-1),1:2]
st = temp[1:(nrow(ice.river)-1)]

#LM, LMresc and Wilks##
set.seed(3010)
T1 = nrow(mY)
ny = ncol(mY)
modvar = vars::VAR(mY, p =1, type = 'const')
AIC1 = AIC(modvar)
BVAR = matrix(0, ncol = ny, nrow = (ny+1))
for(j in 1:ny){
  BVAR[,j] = coef(modvar)[[j]][,1]
}
object0 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(mY[1:(nrow(mY)-1),],1), m = 1, B = BVAR,
               residuals = rbind(rep(0, ny),residuals(modvar)))
LM0 = LMTEST(object0)
starvars::VLSTARjoint(mY, st = st)
LMadj0 = FTEST(LM0, n = ny, m = 1, nX = ncol(object0$Data), iT = T1-1) 
Wilk0 <- wilks(object0)
Panel0temp = cbind(rbind(LM0$LM, LMadj0$LM, Wilk0$Wilks), rbind(LM0$pval, LMadj0$pval, Wilk0$pval))
Panel0temp


####Test vs m = 3
m = 2
stam2 = starting(mY, st = st, n.combi = 30, ncores = 8)
mod = VLSTAR(mY, p = 1, st = st, method = 'NLS', n.iter = 200,
             starting = stam2,
             ncores = 6, constant = T, maxgamma = 50)
AIC2 = sum(mod$AIC)
object1 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1,mY[1:(T1-1),]), m = 2,
               gamma = mod$Gammac[,1], c = mod$Gammac[,2], residuals = resid(mod),
               Gtilde = mod$Gtilde, B = mod$B)
LM1 = LMTEST(object1)
LMadj1 = FTEST(LM1, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
Wilk1 <- wilks(object1)
##ST approach
modtv = tsDyn::TVAR(mY, lag = 1, include = 'const', model = 'TAR', thVar = temp[2:(nrow(ice.river))])
AIC2ST = AIC(modtv)
tlist <- lapply(coef(modtv), t)
Bhat = do.call(cbind, tlist)
objectST = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1,mY[1:(T1-1),]), m = 1, B = Bhat,
                residuals = resid(modtv))
LM1ST = LMTEST(objectST, method = 'TVAR')
LMadj1ST = FTEST(LM1ST, n = ny, m = 1, nX = ncol(objectST$Data), iT = T1-1) 
Wilk1ST = wilks(objectST, method = 'TVAR')
Panel1temp = cbind(rbind(LM1$LM, LMadj1$LM, Wilk1$Wilks), rbind(LM1$pval, LMadj1$pval, Wilk1$pval), 
               rbind(LM1ST$LM, LMadj1ST$LM, Wilk1ST$Wilks), rbind(LM1ST$pval, LMadj1ST$pval, Wilk1ST$pval))
Panel1temp

##Estimate m = 3 (for AIC)#
m = 3
stam3 = starting(mY, st = st, n.combi = 20, ncores = 8, m = 3)
mod2 = VLSTAR(mY, p = 1, st = st, method = 'NLS', n.iter = 100,
              starting = stam3, m = 3,
              ncores = 6, constant = T, maxgamma = 50)
modtv2 = tsDyn::TVAR(mY, lag = 1, include = 'const', model = 'TAR', thVar = st, nthresh = 2, trim = 0.05)
AIC3 = sum(mod2$AIC)
AIC3ST = AIC(modtv2)

#Table 13
c(AIC1, AIC2, AIC3)
c(AIC1, AIC2ST, AIC3ST)

time = tindex
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2)) # 
plot(time, mY[,1], type = "l", col = "black", ylim = range(mY[,1]), ylab = "Jökulsá", xlab = "")
for (j in 1:T1) {
  if (st[j] < getTh(modtv2)[1]) {
    points(time[j], mY[j,1], pch = 2, col = "blue")
  } else if (st[j] > getTh(modtv2)[1] & st[j] <= getTh(modtv2)[2]) {
    points(time[j], mY[j,1], pch = 15, col = "red") 
  } else if(st[j] > getTh(modtv2)[2]){
    points(time[j], mY[j,1], pch = 15, col = "grey")
  }
  
}
plot(time, mY[,2], type = "l", col = "black", ylim = range(mY[,2]), ylab = "Vatndalsá", xlab = "")
for (j in 1:T1) {
  if (st[j] < getTh(modtv2)[1]) {
    points(time[j], mY[j,2], pch = 2, col = "blue")
  } else if (st[j] > getTh(modtv2)[1] & st[j] <= getTh(modtv2)[2]) {
    points(time[j], mY[j,2], pch = 15, col = "red") 
  } else if(st[j] > getTh(modtv2)[2]){
    points(time[j], mY[j,2], pch = 15, col = "grey")
  }
}  
par(mfrow = c(1, 1))

####Table 12 - Panel B####
##st = precipitation##
mY = ice.river[2:nrow(ice.river),1:2]
mX = ice.river[1:(nrow(ice.river)-1),1:2]
st = prec[1:(nrow(ice.river)-1)]

##LM, LMresc and Wilks##
set.seed(3010)
T1 = nrow(mY)
ny = ncol(mY)
modvar = vars::VAR(mY, p =1, type = 'const')
AIC1 = AIC(modvar)
BVAR = matrix(0, ncol = ny, nrow = (ny+1))
for(j in 1:ny){
  BVAR[,j] = coef(modvar)[[j]][,1]
}
object0 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(mY[1:(nrow(mY)-1),],1), m = 1, B = BVAR,
               residuals = rbind(rep(0, ny),residuals(modvar)))
LM0 = LMTEST(object0)
starvars::VLSTARjoint(mY, st = st)
LMadj0 = FTEST(LM0, n = ny, m = 1, nX = ncol(object0$Data), iT = T1-1) 
Wilk0 <- wilks(object0)
Panel0prec = cbind(rbind(LM0$LM, LMadj0$LM, Wilk0$Wilks), rbind(LM0$pval, LMadj0$pval, Wilk0$pval))
Panel0prec

####Test vs m = 3
m = 2
stam2 = starting(mY, st = st, n.combi = 20, ncores = 8)
mod = VLSTAR(mY, p = 1, st = st, method = 'NLS', n.iter = 50,
             starting = stam2,
             ncores = 6, constant = T, maxgamma = 50)
AIC2 = sum(mod$AIC)
object1 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1,mY[1:(T1-1),]), m = 2,
               gamma = mod$Gammac[,1], c = mod$Gammac[,2], residuals = resid(mod),
               Gtilde = mod$Gtilde, B =mod$B)
LM1 = LMTEST(object1)
LMadj1 = FTEST(LM1, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
Wilk1 <- wilks(object1)
##ST approach
modtv = tsDyn::TVAR(mY, lag = 1, include = 'const', model = 'TAR', thVar = prec[2:(nrow(ice.river))])
AIC2ST = AIC(modtv)
tlist <- lapply(coef(modtv), t)
Bhat = do.call(cbind, tlist)
objectST = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1,mY[1:(T1-1),]), m = 1, B = Bhat,
                residuals = resid(modtv))
LM1ST = LMTEST(objectST, method = 'TVAR')
LMadj1ST = FTEST(LM1ST, n = ny, m = 1, nX = ncol(objectST$Data), iT = T1-1) 
Wilk1ST = wilks(objectST, method = 'TVAR')
Panel1prec = cbind(rbind(LM1$LM, LMadj1$LM, Wilk1$Wilks), rbind(LM1$pval, LMadj1$pval, Wilk1$pval), 
               rbind(LM1ST$LM, LMadj1ST$LM, Wilk1ST$Wilks), rbind(LM1ST$pval, LMadj1ST$pval, Wilk1ST$pval))
Panel1prec

##Estimate m = 3 (for AIC)#
m = 3
stam3 = starting(mY, st = st, n.combi = 20, ncores = 8, m = 3)
mod2 = VLSTAR(mY, p = 1, st = st, method = 'NLS', n.iter = 100,
              starting = stam3, m = 3,
              ncores = 6, constant = T, maxgamma = 50)
modtv2 = tsDyn::TVAR(mY, lag = 1, include = 'const', model = 'TAR', thVar = st, nthresh = 2, trim = 0.05)
AIC3 = sum(mod2$AIC)
AIC3ST = AIC(modtv2)

#Table 13
c(AIC1, AIC2, AIC3)
c(AIC1, AIC2ST, AIC3ST)

time = tindex[-1]
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2)) # 
plot(time, mY[,1], type = "l", col = "black", ylim = range(mY[,1]), ylab = "Jökulsá", xlab = "")
for (j in 1:T1) {
  if (st[j] < getTh(modtv2)[1]) {
    points(time[j], mY[j,1], pch = 2, col = "blue")
  } else if (st[j] > getTh(modtv2)[1] & st[j] <= getTh(modtv2)[2]) {
    points(time[j], mY[j,1], pch = 15, col = "red") 
  } else if(st[j] > getTh(modtv2)[2]){
    points(time[j], mY[j,1], pch = 15, col = "grey")
  }
  
}
plot(time, mY[,2], type = "l", col = "black", ylim = range(mY[,2]), ylab = "Vatndalsá", xlab = "")
for (j in 1:T1) {
  if (st[j] < getTh(modtv2)[1]) {
    points(time[j], mY[j,2], pch = 2, col = "blue")
  } else if (st[j] > getTh(modtv2)[1] & st[j] <= getTh(modtv2)[2]) {
    points(time[j], mY[j,2], pch = 15, col = "red") 
  } else if(st[j] > getTh(modtv2)[2]){
    points(time[j], mY[j,2], pch = 15, col = "grey")
  }
}  
par(mfrow = c(1, 1))
