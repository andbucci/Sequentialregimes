library(MASS)
library(doParallel)
library(optimParallel)
library(plotly)
library(starvars)
library(tsDyn)
library(tseries)
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

####Table 10####
mY = data[2:nrow(data),4:5]
mX = data[1:(nrow(data)-1),4:5]
st = data$spreadavg[2:nrow(data)]

####LM, LMresc and Wilks
set.seed(3010)
T1 = nrow(mY)
ny = ncol(mY)
modvar = vars::VAR(mY, p =1, type = 'const')
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
Panel0 = cbind(rbind(LM0$LM, LMadj0$LM, Wilk0$Wilks), rbind(LM0$pval, LMadj0$pval, Wilk0$pval))
Panel0

####Test vs m = 3
data = list(mY = mY, X1 = cbind(1,mX), st = st)
m = 2
nparamH0 = (ncol(mX)+1)*m*2
stam2 = starting(mY, st = st, n.combi = 20, ncores = 8)
ti1 <- c(stam2[[1]][,1],stam2[[1]][,2], rep(0.3, nparamH0))
MLm2 = MLiter(ti = ti1, ll = logl2, data = data, m = 2, epsilon = 1e-3)
object1 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1, mY[1:(T1-1),]), m = 2, gamma = MLm2$ti[1:2],
               c = MLm2$ti[3:4], residuals = t(MLm2$residuals[,-1]), 
               Gtilde = MLm2$Gtilde[-1], B = matrix(MLm2$ti[5:length(ti1)], nrow = 3, ncol=4, byrow = T))
LM1 = LMTEST(object1)
LMadj1 = FTEST(LM1, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
Wilk1 <- wilks(object1)
Panel1 = cbind(rbind(LM1$LM, LMadj1$LM, Wilk1$Wilks), rbind(LM1$pval, LMadj1$pval, Wilk1$pval))
Panel1

###Test vs m = 4
m = 3
nparamH0 = (ncol(mX)+1)*m*2
stam3 = starting(mY, st = st, n.combi = 20, ncores = 8, m = 3)
ti2 <- c(stam3[[1]][,1], stam3[[2]][,1], stam3[[1]][,2], stam3[[2]][,2], 
         rep(0.3, nparamH0))
MLm3 = MLiter(ti = ti2, ll = logl2m3, data = data, m = 3, epsilon = 1e-3)
object2 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1, mY[1:(T1-1),]), m = 3, gamma = MLm3$ti[1:4],
               c = MLm3$ti[5:8], residuals = t(MLm3$residuals[,-1]), 
               Gtilde = MLm3$Gtilde[-1], B = matrix(MLm3$ti[9:length(ti2)], nrow = 3, ncol=6, byrow = T))
LM2 = LMTEST(object2)
LMadj2 = FTEST(LM2, n = ny, m = 3, nX = ncol(object2$Data), iT = T1-1) 
Wilk2 <- wilks(object2)
Panel2 = cbind(rbind(LM2$LM, LMadj2$LM, Wilk2$Wilks), rbind(LM2$pval, LMadj2$pval, Wilk2$pval))
Panel2

####Empirical application on river flows data - Section 6.2####
###Figure 2 and Table 11
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

####Table 11 - Panel A
####st = temperature####
mY = ice.river[2:nrow(ice.river),1:2]
mX = ice.river[1:(nrow(ice.river)-1),1:2]
st = temp[1:(nrow(ice.river)-1)]

####LM, LMresc and Wilks####
set.seed(3010)
T1 = nrow(mY)
ny = ncol(mY)
modvar = vars::VAR(mY, p =1, type = 'const')
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
Panel0 = cbind(rbind(LM0$LM, LMadj0$LM, Wilk0$Wilks), rbind(LM0$pval, LMadj0$pval, Wilk0$pval))
Panel0


####Test vs m = 3
data = list(mY = mY, X1 = cbind(1,mX), st = st)
m = 2
nparamH0 = (ncol(mX)+1)*m*2

stam2 = starting(mY, st = st, n.combi = 20, ncores = 8)
ti1 <- c(stam2[[1]][,1],stam2[[1]][,2], rep(0.3, nparamH0))
MLm2 = MLiter(ti = ti1, ll = logl2, data = data, m = 2, epsilon = 1e-3)
object1 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1, mY[1:(T1-1),]), m = 2, gamma = MLm2$ti[1:2],
               c = MLm2$ti[3:4], residuals = t(MLm2$residuals[,-1]), 
               Gtilde = MLm2$Gtilde[-1], B = matrix(MLm2$ti[5:length(ti1)], nrow = 3, ncol=4, byrow = T))
LM1 = LMTEST(object1)
LMadj1 = FTEST(LM1, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
Wilk1 <- wilks(object1)
Panel1 = cbind(rbind(LM1$LM, LMadj1$LM, Wilk1$Wilks), rbind(LM1$pval, LMadj1$pval, Wilk1$pval))
Panel1

###Test vs m = 4
m=3
nparamH0 = (ncol(mX)+1)*m*2
stam3 = starting(mY, st = st, n.combi = 20, ncores = 8, m =3)
ti2 <- c(stam3[[1]][,1], stam3[[2]][,1], stam3[[1]][,2], stam3[[2]][,2], 
         rep(0.3, nparamH0))
MLm3 = MLiter(ti = ti2, ll = logl2m3, data = data, m = 3, epsilon = 1e-3)
object2 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1, mY[1:(T1-1),]), m = 3, gamma = MLm3$ti[1:4],
               c = MLm3$ti[5:8], residuals = t(MLm3$residuals[,-1]), 
               Gtilde = MLm3$Gtilde[-1], B = matrix(MLm3$ti[9:length(ti2)], nrow = 3, ncol=6, byrow = T))
LM2 = LMTEST(object2)
LMadj2 = FTEST(LM2, n = ny, m = 3, nX = ncol(object2$Data), iT = T1-1) 
Wilk2 <- wilks(object2)
Panel2 = cbind(rbind(LM2$LM, LMadj2$LM, Wilk2$Wilks), rbind(LM2$pval, LMadj2$pval, Wilk2$pval))
Panel2


####st = precipitation####
mY = ice.river[2:nrow(ice.river),1:2]
mX = ice.river[1:(nrow(ice.river)-1),1:2]
st = prec[1:(nrow(ice.river)-1)]

####LM, LMresc and Wilks####
set.seed(3010)
T1 = nrow(mY)
ny = ncol(mY)
modvar = vars::VAR(mY, p =1, type = 'const')
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
Panel0 = cbind(rbind(LM0$LM, LMadj0$LM, Wilk0$Wilks), rbind(LM0$pval, LMadj0$pval, Wilk0$pval))
Panel0


####Test vs m = 3
data = list(mY = mY, X1 = cbind(1,mX), st = st)
m = 2
nparamH0 = (ncol(mX)+1)*m*2
stam2 = starting(mY, st = st, n.combi = 20, ncores = 8)
ti1 <- c(stam2[[1]][,1],stam2[[1]][,2], rep(0.3, nparamH0))
MLm2 = MLiter(ti = ti1, ll = logl2, data = data, m = 2, epsilon = 1e-3)
object1 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1, mY[1:(T1-1),]), m = 2, gamma = MLm2$ti[1:2],
               c = MLm2$ti[3:4], residuals = t(MLm2$residuals[,-1]), 
               Gtilde = MLm2$Gtilde[-1], B = matrix(MLm2$ti[5:length(ti1)], nrow = 3, ncol=4, byrow = T))
LM1 = LMTEST(object1)
LMadj1 = FTEST(LM1, n = ny, m = 2, nX = ncol(object1$Data), iT = T1-1) 
Wilk1 <- wilks(object1)
Panel1 = cbind(rbind(LM1$LM, LMadj1$LM, Wilk1$Wilks), rbind(LM1$pval, LMadj1$pval, Wilk1$pval))
Panel1

###Test vs m = 4
data = list(mY = mY, X1 = cbind(1,mX), st = st)
m=3
nparamH0 = (ncol(mX)+1)*m*2
stam3 = starting(mY, st = st, n.combi = 30, ncores = 8, m = 3)
ti2 <- c(stam3[[1]][,1], stam3[[2]][,1], stam3[[1]][,2], stam3[[2]][,2], 
         rep(0.1, nparamH0))
MLm3 = MLiter(ti = ti2, ll = logl2m3, data = data, m = 3, epsilon = 1e-3)
object2 = list(st = st[2:T1], y = mY[2:T1,], Data = cbind(1, mY[1:(T1-1),]), m = 3, gamma = MLm3$ti[1:4],
               c = MLm3$ti[5:8], residuals = t(MLm3$residuals[,-1]), 
               Gtilde = MLm3$Gtilde[-1], B = matrix(MLm3$ti[9:length(ti2)], nrow = 3, ncol=6, byrow = T))
LM2 = LMTEST(object2)
LMadj2 = FTEST(LM2, n = ny, m = 3, nX = ncol(object2$Data), iT = T1-1) 
Wilk2 <- wilks(object2)
Panel2 = cbind(rbind(LM2$LM, LMadj2$LM, Wilk2$Wilks), rbind(LM2$pval, LMadj2$pval, Wilk2$pval))
Panel2
