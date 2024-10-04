DerGFunc <- function(Psit,BX,ip,im,ms,mgamma,mc)
  # compute the current derivative of Psi_t =
  # ( I, G_t^1(s_t|gamma,c), ... )'
  # input:
  # Psit, mp by p, the current Psit
  # BX, mp by 1, vector, BX = t(mB)%*%Xt
  # ms, p by m-1, vector of the current transition variable s_t
  # mgamma, p by m-1
  # mc, p by m-1
  # output:
  # dPsi_t/ddelta, 2*p*p*(m-1) vector
{
  DPsit = Psit - Psit*Psit
  DPsit1 = DPsit %*% (ms-mc)
  DPsit2 = -DPsit %*% mgamma
  tDP = matrix(0,(ip+1),(im*ip))
  RET = NULL
  for(iter in 1:(im-1)){
    for(jter in 1:ip){
      tDP[iter,(iter*ip+jter)] = DPsit1[(iter*ip+jter),iter]
      RET = c(RET, c(tDP%*%BX))
      tDP[iter,(iter*ip+jter)] = 0
      tDP[iter,(iter*ip+jter)] = DPsit2[(iter*ip+jter),iter]
      RET = c(RET, c(tDP%*%BX))
      tDP[iter,(iter*ip+jter)] = 0
    }
  }
  return(RET)
}

LMTEST <- function(object, d = 3){
  mY = as.matrix(object$y)
  mX = as.matrix(object$Data)
  st <- object$st
  x <- as.matrix(object$Data)
  ncolx <- ncol(x)
  nrowx <- nrow(x)
  ncoly <- ncol(object$y)
  m <- object$m
  im = m
  ZZ <- matrix(nrow = nrowx, ncol = ncolx*3)
  for (i in 1:nrowx){
    xst1 <-  as.matrix(x[i,]*st[i])
    xst2 <- as.matrix(x[i,]*st[i]^2)
    xst3 <- as.matrix(x[i,]*st[i]^3)
    ZZ[i,] <- cbind(xst1, xst2, xst3)
  }
  mZ = ZZ[,(1:(d*ncolx))]
  #mY, mX, mZ must be matrices!
  #returns the LM test statistic and the degree of freedom
  iT = dim(mY)[1]
  ip = dim(mY)[2]
  iDF = dim(mZ)[2]*ip
  #ee <- object$residual[2:nrow(object$residuals),]
  ee <- object$residual
  Omega <- (t(ee) %*% ee)/(nrowx-ncoly)
  It <- diag(nrowx)
  ###Derivative of Gt
  if(object$m >1){
    dimgit = dim(object$Gtilde[[1]])[1]
    mK = NULL
    for(t in 1:nrowx){
      Psit = object$Gtilde[[t]]
      BX = t(object$B)%*%x[t,]
      temp = c(Psit)
      dim(temp) = c(1,im*ip*ip)
      Kt = c(x[t,]%*%temp)
      Kt = c(Kt,DerGFunc(Psit, BX = BX, ip = ncoly, ms = rep(st[t], ncoly), im = object$m,
                         mgamma = matrix(object$gamma, nrow = ncoly), mc =  matrix(object$gamma, nrow = ncoly)))
      mK = rbind(mK,Kt)
    }
    mU = svd(mK)$u
    #mE = mY - mX%*%solve(t(mX)%*%mX)%*%t(mX)%*%mY
    mE = ee - (mU%*%t(mU))%*%ee
    mXX = cbind(mK, mZ)
  }else{
    mE = mY - mX%*%MASS::ginv(t(mX)%*%mX,tol = 1e-20)%*%t(mX)%*%mY
    mXX = cbind(mX, mZ)
    mk = mE - mXX%*%MASS::ginv(t(mXX)%*%mXX,tol = 1e-20)%*%t(mXX)%*%mE
  }
  #mE = ee
  RSS0 = t(mE)%*%mE
  mK = mE - mXX%*%ginv(t(mXX)%*%mXX)%*%t(mXX)%*%mE
  RSS1 = t(mK)%*%mK
  dTR = sum(diag(solve(RSS0)%*%RSS1))
  LM = iT*(ip-dTR)
  pval = 1-pchisq(LM,df=iDF)
  return(list(LM = LM, pvalue= pval, df = iDF))
}


FTEST <- function(LM, iT, m, n, nX)
{
  iDF1 = LM[[3]]
  iK = 2*n*m + m*nX
  iDF2 = n*(iT-iK)
  FT = LM[[1]]*(iT-iK)/(iT*LM[[3]])
  pval = 1-pf(FT,df1=iDF1,df2=iDF2)
  return(list(LM = FT, pvalue = pval, df1 = iDF1, df2 = iDF2))
}

wilks <- function(object, d = 3, method = 'LVSTAR'){
  mY = as.matrix(object$y)
  mX = as.matrix(object$Data)
  st <- object$st
  x <- as.matrix(object$Data)
  ncolx <- ncol(x)
  nrowx <- nrow(x)
  ncoly <- ncol(object$y)
  m <- object$m
  im = m
  ZZ <- matrix(nrow = nrowx, ncol = ncolx*3)
  for (i in 1:nrowx){
    xst1 <-  as.matrix(x[i,]*st[i])
    xst2 <- as.matrix(x[i,]*st[i]^2)
    xst3 <- as.matrix(x[i,]*st[i]^3)
    ZZ[i,] <- cbind(xst1, xst2, xst3)
  }
  mZ = ZZ[,(1:(d*ncolx))]
  #mY, mX, mZ must be matrices!
  #returns the LM test statistic and the degree of freedom
  iT = dim(mY)[1]
  ip = dim(mY)[2]
  iz = dim(mZ)[2]
  ix = ncolx
  iDF = iz*ip
  #ee <- object$residual[2:nrow(ee),]
  ee <- object$residual
  Omega <- (t(ee) %*% ee)/(nrowx-ncoly)
  It <- diag(nrowx)
  if(method == 'LVSTAR'){
    if(object$m > 1){
      ###Derivative of Gt
      dimgit = dim(object$Gtilde[[1]])[1]
      mK = NULL
      for(t in 1:nrowx){
        Psit = object$Gtilde[[t]]
        BX = t(object$B)%*%x[t,]
        temp = c(Psit)
        dim(temp) = c(1,im*ip*ip)
        Kt = c(x[t,]%*%temp)
        Kt = c(Kt,DerGFunc(Psit, BX = BX, ip = ncoly, ms = rep(st[t], ncoly), im = object$m,
                           mgamma = matrix(object$gamma, nrow = ncoly), mc =  matrix(object$gamma, nrow = ncoly)))
        mK = rbind(mK,Kt)
      }
      mU = svd(mK)$u
      #mE = mY - mX%*%solve(t(mX)%*%mX)%*%t(mX)%*%mY
      mE = ee - (mU%*%t(mU))%*%ee
      mXX = cbind(mK, mZ)
    }else{
      mE = mY - mX%*%solve(t(mX)%*%mX,tol = 1e-20)%*%t(mX)%*%mY
      mXX = cbind(mX, mZ)
      mK = mE - mXX%*%solve(t(mXX)%*%mXX,tol = 1e-20)%*%t(mXX)%*%mE
    }
  } else if(method == 'TVAR'){
    mE = ee
    mXX = cbind(mX, mZ)
  }
  #mE = ee
  RSS0 = t(mE)%*%mE
  mK1 = mE - mXX%*%ginv(t(mXX)%*%mXX)%*%t(mXX)%*%mE
  #mK1 = mE - mK%*%solve(t(mK)%*%mK)%*%t(mK)%*%mE
  RSS1 = t(mK1)%*%mK1
  R0 = svd(RSS0)$d
  R1 = svd(RSS1)$d
  Lambda = sum(log(R1))-sum(log(R0))
  Lambda = Lambda * ((ip+iz+1)*.5 + ix - iT)
  pval = 1-pchisq(Lambda,df=iDF)
  return(list(Wilks = Lambda, pvalue = pval, df = iDF))
}
