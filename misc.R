#'Multivariate log-likelihood
#'
#'Log-likelihood to be optimized in the ML method and used to check convergence in both methods
#'@keywords internal
#'
loglike <- function(param, data){
  y = data$y
  Omegahat = data$Omegahat
  gamma <- param[1:((data$m-1)*ncol(y))]
  c <- param[(ncol(y)*(data$m-1)+1):length(param)]
  gamma1 <- matrix(gamma, ncol = (data$m-1))
  c1 <- matrix(c, ncol = (data$m-1))
  glog <- rep(0, ncol(y))
  GT <- list()
  dify <- matrix(ncol = 1, nrow = nrow(y))
  for (z in 1:nrow(y)){
    for(t in 1:(data$m-1)){
      for (o in 1:ncol(y)){
        glog[o] <- (1L+exp(-gamma1[o,t]*(data$st[z]-c1[o,t])))^(-1)}
      if(data$singlecgamma == TRUE){
        GT[[t]] <- diag(rep(glog[1], ncol(y)))
      }else{
        GT[[t]] <- diag(glog)
      }
    }
    Gtilde <- t(cbind(diag(ncol(y)), do.call(cbind,GT)))
    dify[z] <-  t(y[z, ] - t(Gtilde)%*%t(data$BB)%*%data$x[z,])%*%MASS::ginv(Omegahat)%*%(y[z, ] - t(Gtilde)%*%t(data$BB)%*%data$x[z,])
  }
  sumdif <- sum(dify)
  logll <- -(nrow(y)*log(det(Omegahat))/2L) - sumdif/2L  - (nrow(y)*ncol(y)/2L)*log(2L*pi)##Normal distribution assumed
  return(-logll)
}

#'Sum of squared error
#'
#'Sum of squared error to be optimized in the NLS method
#'@keywords internal
#'
SSQ <- function(param, data){
  y <- data$y
  gamma <- param[1:((data$m-1)*ncol(y))]
  c <- param[(ncol(y)*(data$m-1)+1):length(param)]
  gamma1 <- matrix(gamma, ncol = (data$m-1))
  c1 <- matrix(c, ncol = (data$m-1))
  glog <- rep(0, ncol(y))
  GT <- list()
  Gtilde <- matrix(ncol = ncol(y), nrow = (ncol(y)+ncol(y)*(data$m-1)))
  dify <- matrix(ncol = 1, nrow = nrow(y))
  for (z in 1:nrow(y)){
    for(t in 1:(data$m-1)){
      for (o in 1:ncol(y)){
        glog[o] <- (1L+exp(-gamma1[o,t]*(data$st[z]-c1[o,t])))^(-1)}
      if(data$singlecgamma == TRUE){
        GT[[t]] <- diag(rep(glog[1], ncol(y)))
      }else{
        GT[[t]] <- diag(glog)
      }
    }
    Gtilde <- t(cbind(diag(ncol(y)), do.call(cbind,GT)))
    dify[z] <-  t(y[z, ] - t(Gtilde)%*%t(data$BB)%*%data$x[z,])%*%(y[z, ] - t(Gtilde)%*%t(data$BB)%*%data$x[z,])
  }
  sumdif <- sum(dify)
  return(sumdif)
}


SSQsingle <- function(param, data){
  y <- data$y
  gamma <- param[1:((data$m-1)*ncol(y))]
  c <- param[(ncol(y)*(data$m-1)+1):length(param)]
  gamma1 <- matrix(gamma, ncol = (data$m-1))
  c1 <- matrix(c, ncol = (data$m-1))
  glog <- rep(0, ncol(y))
  GT <- list()
  Gtilde <- matrix(ncol = ncol(y), nrow = (ncol(y)+ncol(y)*(data$m-1)))
  dify <- matrix(ncol = 1, nrow = nrow(y))
  for (z in 1:nrow(y)){
    for(t in 1:(data$m-1)){
      for (o in 1:ncol(y)){
        glog[o] <- (1L+exp(-gamma1[o,t]*(data$st[z]-c1[o,t])))^(-1)}
      if(data$singlecgamma == TRUE){
        GT[[t]] <- diag(rep(glog[o], ncol(y)))
      }else{
        GT[[t]] <- diag(glog)
      }
    }
    Gtilde <- t(cbind(diag(ncol(y)), do.call(cbind,GT)))
    dify[z] <-  t(y[z, ] - t(Gtilde)%*%t(data$BB)%*%data$x[z,])%*%(y[z, ] - t(Gtilde)%*%t(data$BB)%*%data$x[z,])
  }
  sumdif <- sum(dify)
  return(sumdif)
}

VLSTAR <- function(y, exo = NULL, p = 1,
                   m = 2, st = NULL, constant = TRUE, starting = NULL,
                   method = c('ML', 'NLS'), n.iter = 500,
                   singlecgamma = FALSE,
                   epsilon = 10^(-3), ncores = NULL, maxgamma = 20){
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if(is.null(ncores)){
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      ncores <- 2L
    } else {
      ncores <- parallel::detectCores()
    }}
  y <- as.matrix(y)
  x <- exo
  method <- match.arg(method)
  ##Checks and warnings
  if (anyNA(y))
    stop("\nNAs in y.\n")
  if(m < 2)
    stop('The number of regimes should be greater than one.')
  if(is.null(st))
    stop('The transition variable must be supplied.')
  if(is.null(x)){
    if(length(y[,1]) != length(st))
      stop('The length of the variables does not match!')
  }else{
    if(length(y[,1]) != length(as.matrix(x[,1])) | length(st) != length(as.matrix(x[,1])) | length(y[,1]) != length(st))
      stop('The length of the variables does not match!')
  }
  if(is.null(starting)){
    stop('Starting values should be provided.')
  }
  if(!is.list(starting)){
    stop('Starting c and gamma should be put in a list.')
  }
  if(is.null(p)){
    stop('Please, specify a valid lag order.')
  }
  if (ncol(y) < 2)
    stop("The matrix 'y' should contain at least two variables. For univariate analysis consider lstar() function in this package.\n")
  if (is.null(colnames(y))) {
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    warning(paste("No column names supplied in y, using:",
                  paste(colnames(y), collapse = ", "), ", instead.\n"))
  }
  colnames(y) <- make.names(colnames(y))
  ##Definition of dimensions, creating variable x with constant
  yt <- zoo::zoo(y)
  if(p == 0){
    ylag = NULL
  }else{
    ylag <- stats::lag(yt, -(1:p))
    ylag <- as.matrix(ylag)
    if(p>1){
      lagg <- p-1
      ylag <- ylag[-(1:lagg),]
    }
  }
  y <- y[-c(1:p), ]
  ncoly <- ncol(y)
  if(singlecgamma == TRUE){
    starting1 <- list()
    for(i in 1:(m-1)){
      starting1[[i]] <- cbind(rep(starting[[i]][,1], ncoly), rep(starting[[i]][,2], ncoly))
    }
    starting = starting1
  }
  if(!is.null(starting)){
    if(length(starting)!= (m-1)){
      stop('The length of the list of initial values should be equal to m-1.')
    }else{
      if(any(unlist(lapply(starting, ncol))!=2) | any(unlist(lapply(starting, nrow))!=ncoly)){
        stop('Each element of the starting argument should have two columns and n rows.')
      }
      
    }
  }
  ncolylag <- ncoly*p
  nrowy <- nrow(y)
  ncolx1 <- ncol(x)
  const <- rep(1, nrowy)
  if (constant == TRUE){
    if(!is.null(exo)){
      x1a <- as.matrix(x[-c(1:p),])
      x <- as.matrix(cbind(const,ylag,x1a))
    }else{
      x <- as.matrix(cbind(const,ylag))
    }
    ncolx <- ncol(x)
  }  else{if(!is.null(exo)){
    x1a <- as.matrix(x[-c(1:p),])
    x <- as.matrix(cbind(ylag,x1a))
  }
    x <- as.matrix(ylag)
  }
  nrowx <- nrow(x)
  st <- st[(1+p):length(st)]
  ncolx <- ncol(x)
  ny <- ifelse(singlecgamma == TRUE, 1, ncoly)
  q <- ncol(x)-ncolylag
  PARAM <- starting
  
  ####Estimating VLSTAR model####
  
  ##Estimating initial values to be used in the iterative algorithm
  In <- diag(ncoly)
  glog <- matrix(ncol=ny, nrow = nrowy)
  GT <- list()
  Gtilde <- list()
  kro <- list()
  for (i in 1:nrowx){
    for (t in 1:(m-1)){
      for (j in 1 : ny){
        glog[i,j] <- (1L+exp(-PARAM[[t]][j,1]*(st[i]-PARAM[[t]][j,2])))^(-1)
      }
      if(singlecgamma == TRUE){
        Gt <- diag(rep(glog[i,1], ncoly))
      }else{
        Gt <- diag(glog[i,])
      }
      GT[[t]] <- Gt
    }
    Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
    kro[[i]] <- kronecker(Gtilde[[i]], x[i,])
  }
  M <- t(do.call("cbind", kro))
  Y <- matrixcalc::vec(t(y))
  Bhat <- MASS::ginv(t(M)%*%M)%*%t(M)%*%Y
  BB <- ks::invvec(Bhat, ncol = (m*ncoly), nrow = (ncolylag + q)) ##Estimated coefficients
  ##Calculating the estimated covariance matrix
  #resi <- list()
  Ehat1 <- matrix(NA, ncol = ncoly, nrow = nrowy)
  for (i in 1:nrowx){
    #resi[[i]] <- y[i, ] - t(Gtilde[[i]])%*%t(BB)%*%x[i,]
    Ehat1[i,] <- t(y[i, ] - t(Gtilde[[i]])%*%t(BB)%*%x[i,])
  }  #Estimated covariance matrix
  Omegahat <- (t(Ehat1)%*%Ehat1)/(nrowy)
  
  data = list(y = y, x = x, m = m, BB = BB, Omegahat = Omegahat, st = st, singlecgamma = singlecgamma)
  
  #Inizialization of iter
  iter <- 0
  ll0 <- 10^(4)
  loglik1 <- NULL
  omega <- list()
  PARAM1 <- list()
  for(t in 1:(m-1)){
    PARAM1[[t]] <- as.data.frame(PARAM[[t]])
  }
  param <- do.call(rbind,PARAM1)
  param <- matrixcalc::vec(as.matrix(param)) ##the parameters c and gamma are vectorized in order to be passed in the optimParallel function
  err <- 10^5
  
  ##Actual iteration to estimate the coefficients
  # cl <- parallel::makeCluster(ncores, outfile="")     # set the number of processor cores
  # parallel::setDefaultCluster(cl=cl)
  if(requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(ncores, outfile="")
    parallel::setDefaultCluster(cl=cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  if(method == 'ML'){
    #NLS Estimation of Bhat and Omegahat to be used in the first iteration of maximum likelihood
    message('Maximum likelihood estimation\n')
    #Convergence algorithm
    #1. Maximum likelihood estimation of gamma and c with NLS estimates of Bhat and Omegahat
    #2. Maximum likelihood estimation of Bhat with new values of gamma and c
    #3. Convergence check
    #4. 1-2-3 until convergence
    while (iter < n.iter & err > epsilon){
      Sys.sleep(0)
      iter <- iter+1
      #Parameters
      low1 <- replicate(ny, 0)
      min1 <- replicate(ncoly, min(st))
      up1 <- replicate(ncoly, maxgamma)
      max1 <- replicate(ncoly, max(st))
      #1.Maximum likelihood estimation of gamma and c
      param1 <- optimParallel::optimParallel(par = as.vector(param), fn = loglike, lower = c(low1, min1),
                                             upper = c(up1, max1),
                                             data = data, parallel = list(cl = cl, forward = FALSE, loginfo = FALSE))
      #parallel::stopCluster(cl)
      if(singlecgamma == TRUE){
        cgam1 <- matrix(rep(param1$par, ncoly), ncol = 2L, nrow = (ncoly*(m-1)))
      }else{
        cgam1 <- matrix(param1$par, ncol = 2L, nrow = (ncoly*(m-1)))
      }
      
      #2.Maximum likelihood estimation of Bhat with new values of gamma and c
      glog <- rep(0, ncol(y))
      GT <- list()
      Gtilde <- list()
      XX <- list()
      XYOG <- list()
      kro <- list()
      PsiOmegaPsi <- list()
      for (i in 1:nrowx){
        for(t in 1:(m-1)){
          for (j in 1 : ny){
            glog[j] <- (1+exp(-cgam1[j,1]*(st[i]-cgam1[j,2])))^(-1)
          }
          if(singlecgamma == TRUE){
            GT[[t]] <- diag(rep(glog[1], ncoly))
          }else{
            GT[[t]] <- diag(glog)
          }
        }
        Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
        XX[[i]] <- x[i,] %*%t(x[i,])
        XYOG[[i]] <- matrixcalc::vec((x[i, ]%*%t(y[i,]))%*%MASS::ginv(Omegahat)%*%t(Gtilde[[i]]))
        PsiOmegaPsi[[i]] <- Gtilde[[i]]%*%MASS::ginv(Omegahat)%*%t(Gtilde[[i]])
        kro[[i]] <- kronecker(PsiOmegaPsi[[i]], XX[[i]])
      }
      xyog <- Reduce(`+`, XYOG)/nrowy
      kroxx <- Reduce(`+`, kro)/nrowy
      Bhat <- t(t(xyog)%*%kroxx)
      BB <- ks::invvec(Bhat, ncol = (m*ncoly), nrow = (ncolylag + q))
      resi <- list()
      fitte <- matrix(nrow = nrowy, ncol = ncoly)
      Ehat1 <- matrix(NA, ncol = ncoly, nrow = nrowy)
      for (i in 1:nrowx){
        resi[[i]] <- y[i, ] - t(Gtilde[[i]])%*%t(BB)%*%x[i,]
        fitte[i,] <- t(t(Gtilde[[i]])%*%t(BB)%*%x[i,])
      }
      Ehat1 <- t(do.call("cbind", resi))
      Omegahat <- (t(Ehat1)%*%Ehat1)/(nrowy)
      data$Omegahat = Omegahat
      data$BB = BB
      #3. Convergence check
      ll1 <- loglike(param = param1$par, data)
      if((ll1-ll0)<0){
        param <- as.matrix(matrixcalc::vec(cgam1))
      }
      err <- abs(ll1 - ll0)
      ll0 <- ll1
      loglik1[iter] <- ll1
      
      message(paste("iteration", iter, "complete\n"))
      
      cat(paste('Log-likelihood:', round(ll0,3), '\n'))
      
      if (err<epsilon | iter == n.iter) message('Converged\n')}
  } else{
    message('NLS estimation\n')
    
    #Inizialization of iter
    #Convergence algorithm
    #1. NLS estimation of gamma and c with NLS estimates of Bhat and Omegahat
    #2. NLS of Bhat with new values of gamma and c
    #3. Convergence check
    #4. 1-2-3 until convergence
    while (iter < n.iter & err > epsilon){
      Sys.sleep(0)
      iter <- iter+1L
      #Parameters
      low1 <- replicate(ny, 0)
      min1 <- replicate(ncoly, min(st))
      up1 <- replicate(ncoly, maxgamma)
      max1 <- replicate(ncoly, max(st))
      #1.Maximum likelihood estimation of gamma and c
      param1 <- optimParallel::optimParallel(par = as.vector(param), fn = SSQ, lower = c(low1, min1),
                                             upper = c(up1, max1),
                                             data = data, parallel = list(cl = cl, forward = FALSE, loginfo = FALSE))
      #parallel::stopCluster(cl)
      if(singlecgamma == TRUE){
        cgam1 <- matrix(param1$par, ncol = 2L, nrow = (ncoly*(m-1))) 
      }else{
        cgam1 <- matrix(param1$par, ncol = 2L, nrow = (ncoly*(m-1)))  
      }
      
      #2.NLS estimation of Bhat with new values of gamma and c
      glog <- rep(0, ncol(y))
      GT <- list()
      Gtilde <- list()
      kro <- list()
      for (i in 1:nrowx){
        for(t in 1:(m-1)){
          for (j in 1 : ny){
            glog[j] <- (1+exp(-cgam1[((t-1)*ny + j),1]*(st[i]-cgam1[((t-1)*ny + j),2])))^(-1)
          }
          if(singlecgamma == TRUE){
            GT[[t]] <- diag(rep(glog[1], ncoly))
          }else{
            GT[[t]] <- diag(glog)
          }
        }
        Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
        kro[[i]] <- kronecker(Gtilde[[i]], x[i,])
      }
      Bhat <- MASS::ginv(t(t(do.call("cbind", kro)))%*%(t(do.call("cbind", kro))))%*%t(t(do.call("cbind", kro)))%*%(matrixcalc::vec(t(y)))
      BB <- ks::invvec(Bhat, ncol = (m*ncoly), nrow = (ncolylag + q))
      resi <- list()
      resiresi <- list()
      fitte <- matrix(nrow = nrowy, ncol = ncoly)
      for (o in 1:nrowx){
        resi[[o]] <- y[o, ] - t(Gtilde[[o]])%*%t(BB)%*%x[o,]
        resiresi[[o]] <- resi[[o]]%*%t(resi[[o]])
        fitte[o,] <- t(t(Gtilde[[o]])%*%t(BB)%*%x[o,])
      }
      Ehat1 <- Reduce("+", resiresi)
      Omegahat <- Ehat1/(nrowy-1L)
      
      data$Omegahat = Omegahat
      data$BB = BB
      
      #3. Convergence check
      ll1 <- loglike(param = param1$par, data)
      if((ll1-ll0)<0){
        param <- as.matrix(matrixcalc::vec(cgam1))
      }
      err <- abs(ll1 - ll0)
      ll0 <- ll1
      loglik1[iter] <- ll1
      
      message(paste("iteration", iter, "complete\n"))
      
      cat(paste('Log-likelihood:',round(ll0, 3), '\n'))
      
      if(iter > 50){
        if(loglik1[[iter]] == loglik1[[(iter-2)]]){
          iter <- n.iter
        }}
      
      if (err<epsilon | iter == n.iter) message('Converged\n')}
  }
  
  
  #Calculating residuals and estimating standard errors
  residuals1 <- t(do.call("cbind", resi))
  varhat <- diag(Omegahat)
  bb1 <- BB[,1:ncoly]
  bb2 <- list()
  for(t in 1:(m-1)){
    bb2[[t]] <- as.data.frame(BB[,(ncoly*(t-1)+1):(ncoly*t)] +
                                BB[,(ncoly*(t)+1):(ncoly*(t+1))])
  }
  bb4 <- as.matrix(do.call(rbind,bb2))
  BBhat <- rbind(bb1, bb4) ##Estimates of the coefficients
  
  ##Calculating standard errors, t-test and p-values
  colnames(cgam1) <- c('gamma', 'c')
  covbb <- matrix(nrow = m*ncolx, ncol = ncoly)
  ttest <- matrix(nrow = m*ncolx, ncol = ncoly)
  pval <- matrix(nrow = m*ncolx, ncol = ncoly)
  ee <- matrix(nrow = m*ncolx, ncol = ncoly)
  for (j in 1 : ncoly){
    covbb[,j] <- sqrt(diag(MASS::ginv(t(x) %*%x))*varhat[j])
    ttest[,j] <- BBhat[,j]/covbb[,j]
    pval[,j] <- 2*pt(abs(ttest[,j]),df=(nrowy*m-m*(ncolx)), lower.tail = FALSE)
  }
  ##Significances used for the summary function
  signifi <- matrix('',nrow = nrow(pval), ncol = ncol(pval))
  for(i in 1:nrow(pval)){
    for(j in 1:ncol(pval)){
      if(pval[i,j] < 0.01)
      {signifi[i,j] <- '***'}
      else if(pval[i,j] <=0.05 & pval[i,j]>0.01){
        signifi[i,j] <- '**'
      }
      else if(pval[i,j] <= 0.1 & pval[i,j]>0.05){
        signifi[i,j] <- '*'
      }
    }
  }
  
  ##Calculating univariate log-likelihoods, AIC and BIC criteria
  residui <- residuals1
  omega1 <- diag(Omegahat)
  k <- nrow(BBhat)
  ll2 <- NULL
  AIC1 <- NULL
  BIC1 <- NULL
  for (l in 1:ncoly){
    ll2[l] <- -(nrow(y)/2)*log(omega1[l]) - (t(residui[,l])%*%residui[,l])/(2*omega1[l])
    AIC1[l] <- 2*k - 2*ll2[l]
    BIC1[l] <- -2*ll2[l] + k*log(nrowy)
  }
  names1 <- list()
  for(j in 1:m){
    names1[[j]] <- as.data.frame(paste(colnames(x), ' m_', j, sep = ''))
  }
  names1 <- as.matrix(do.call(rbind,names1))
  rownames(BBhat) <- names1
  colnames(BBhat) <- colnames(y)
  regime_matrices <- list()
  for (regime in 1:m) {
         start_row <- (regime - 1) * ncol(x) + 1
         end_row <- start_row + ncol(x) - 1
         regime_matrices[[regime]] <- BBhat[start_row:end_row, , drop = FALSE]
  }
  BB <- do.call(cbind, regime_matrices)
  modeldata <- list(y, x)
  fitte <- fitte[!is.na(fitte[,1]),]
  results <- list(BBhat, covbb, ttest, pval, cgam1, Omegahat, fitte, residuals1, ll1, ll2, AIC1, BIC1, Gtilde, modeldata, BB, m, p,
                  st, y, exo, constant, method, singlecgamma)
  names(results) <- c('Bhat','StDev', 'ttest', 'pval', 'Gammac', 'Omega', 'fitted', 'residuals', 'MultiLL', 'LL', 'AIC',
                      'BIC', 'Gtilde', 'Data', 'B', 'm', 'p', 'st', 'yoriginal', 'exo', 'constant', 'method', 'singlecgamma')
  class(results) = 'VLSTAR'
  return(results)
}


starting <- function(y, exo = NULL, p = 1,
                     m = 2, st = NULL, constant = TRUE,
                     n.combi = NULL, ncores = 2,
                     singlecgamma = FALSE){
  y <- as.matrix(y)
  x <- exo
  if (anyNA(y))
    stop("\nNAs in y.\n")
  if(m < 2)
    stop('The number of regimes should be greater than one.')
  if(is.null(st))
    stop('The transition variable must be supplied.')
  if(is.null(x)){
    if(length(y[,1]) != length(st))
      stop('The length of the variables does not match!')
  }else{
    if(length(y[,1]) != length(as.matrix(x[,1])) | length(st) != length(as.matrix(x[,1])) | length(y[,1]) != length(st))
      stop('The length of the variables does not match!')
  }
  if(is.null(n.combi)){
    stop('A number of combinations should be provided.')
  }
  
  if(is.null(p) || p < 1){
    stop('Please, specify a valid lag order.')
  }
  if (ncol(y) < 2)
    stop("The matrix 'y' should contain at least two variables. For univariate analysis consider lstar() function in this package.\n")
  if (is.null(colnames(y))) {
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    warning(paste("No column names supplied in y, using:",
                  paste(colnames(y), collapse = ", "), ", instead.\n"))
  }
  colnames(y) <- make.names(colnames(y))
  ##Definition of dimensions, creating variable x with constant
  yt <- zoo(y)
  ylag <- as.matrix(stats::lag(yt, -(1:p)))
  if(p>1){
    lagg <- p-1
    ylag <- ylag[-(1:lagg),]
  }
  y <- y[-c(1:p), ]
  ncoly <- ncol(y)
  In <- diag(ncoly)
  ncolylag <- ncoly*p
  nrowy <- nrow(y)
  ncolx1 <- ncol(x)
  const <- rep(1, nrowy)
  if (constant == TRUE){
    if(!is.null(exo)){
      x1a <- as.matrix(x[-c(1:p),])
      x <- as.matrix(cbind(const,ylag,x1a))
    }else{
      x <- as.matrix(cbind(const,ylag))
    }
    ncolx <- ncol(x)
  } else{
    if(!is.null(exo)){
      x1a <- as.matrix(x[-c(1:p),])
      x <- as.matrix(cbind(ylag,x1a))
    }else{
      x <- as.matrix(ylag)
    }
  }
  nrowx <- nrow(x)
  st <- st[(1+p):length(st)]
  ncolx <- ncol(x)
  param.init <- list()
  # Starting values for c and gamma
  if (singlecgamma == TRUE) {
    param.init$gamma <- 2L
    param.init$c <- mean(y)
  } else {
    param.init$gamma <- rep(2L, ncoly)
    param.init$c <- colMeans(y)
  }
  
  ny <- ifelse(singlecgamma == TRUE, 1, ncoly)
  COMBI <- list()
  GAMMA <- list()
  CJ <- list()
  for (t in 1:(m - 1)) {
    if (singlecgamma == TRUE) {
      gamma <- seq(from = 2L, to = 20, length.out = n.combi)
      cj <- seq(from = rangey[1] * 1.1, to = rangey[2], length.out = n.combi)
      GAMMA[[t]] <- gamma
      CJ[[t]] <- cj
      COMBI[[t]] <- expand.grid(CJ[[t]], GAMMA[[t]])
    } else {
      GAMMA[[t]] <- matrix(ncol = ny, nrow = n.combi)
      CJ[[t]] <- matrix(ncol = ny, nrow = n.combi)
      rangey <- matrix(nrow = ny, ncol = 2)
      combi <- list()
      for (j in 1:ny) {
        rangey[j, ] <- range(y[, j])
        CJ[[t]][, j] <- seq(from = rangey[j, 1] * 1.1, to = rangey[j, 2], length.out = n.combi)*1/(1+exp(-t))
        GAMMA[[t]][, j] <- seq(from = 2L, to = 20, length.out = n.combi)*(1+exp(-t))
        combi[[j]] <- expand.grid(CJ[[t]][, j], GAMMA[[t]][, j])
      }
      COMBI[[t]] <- combi
    }
  }
  
  #NLS for each combination of c and gamma
  
  message(paste('Searching optimal c and gamma among', n.combi*n.combi, 'combinations\n'))
  cl <- parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = n.combi*n.combi, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  l = 1
  ssq <- foreach(l = 1:(n.combi*n.combi), .errorhandling='pass', .combine = rbind,.options.snow = opts) %dopar% {
    Sys.sleep(0.1)
    glog <- matrix(ncol=ny, nrow = nrowy)
    GT <- list()
    Gtilde <- list()
    kro <- list()
    for (i in 1:nrowx){
      for(t in 1:(m-1)){
        if(singlecgamma == TRUE){
          glog[i,j] <- (1+exp(-COMBI[[t]][l,2]*(st[i]-COMBI[[t]][l,1])))^(-1)
          Gt <- diag(rep(glog[i,1], ncoly))
        }else{
          for (j in 1 : ny){
            glog[i,j] <- (1+exp(-COMBI[[t]][[j]][l,2]*(st[i]-COMBI[[t]][[j]][l,1])))^(-1)
          }
          Gt <- diag(glog[i,])
        }
        GT[[t]] <- Gt
      }
      Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
      kro[[i]] <- kronecker(Gtilde[[i]], x[i,])
    }
    M <- t(do.call("cbind", kro))
    Y <- matrixcalc::vec(t(y))
    Bhat <- MASS::ginv(t(M)%*%M)%*%t(M)%*%Y
    BB <- ks::invvec(Bhat, ncol = (m*ncoly), nrow = nrow(Bhat)/(m*ncoly))
    Ehat <- matrix(NA, ncol = ncoly, nrow = nrowy)
    for (o in 1:nrowx){
      Ehat[o,] <- t(y[o, ] - t(Gtilde[[o]])%*%t(BB)%*%x[o,])
    }
    SSQ <- colSums(Ehat^2)
    return(SSQ)
  }
  close(pb)
  stopCluster(cl)
  #c and gamma minimazing the sum of squared residuals for each equation
  if(singlecgamma == TRUE){
    cgamma <- vector(length = 2L)
  }else{
    cgamma <- matrix(nrow = ny, ncol = 2)
  }
  
  CGAMMA <- list()
  for (t in 1:(m-1)){
    if(singlecgamma == T){
      cgamma <-  as.matrix(COMBI[[t]][which.min(rowSums(ssq)),])
    } else{
      for(j in 1:ny){
        cgamma[j,] <- as.matrix(COMBI[[t]][[j]][which.min(ssq[,j]),])
        #cgamma[j,] <- as.matrix(COMBI[[t]][[j]][which.min(rowSums(ssq)),])
      }
    }
    CGAMMA[[t]] <- cgamma
  }
  
  
  #Definition of c0 gamma0
  PARAM <- list()
  for (t in 1:(m-1)){
    cj <- as.matrix(CGAMMA[[t]][,1])
    gamma <- as.matrix(CGAMMA[[t]][,2])
    PARAM[[t]] <- cbind(gamma, cj)
  }
  results <- PARAM
  class(results) = 'startingVLSTAR'
  return(results)
}

TVAR.simexo <- function(n, nsim, m, phi, c){
  In = diag(n)
  nsimu = nsim+100
  y = matrix(0, nrow = nsimu, ncol = n)
  epsilon = MASS::mvrnorm(nsimu, mu = rep(0, n), Sigma = In)
  y[1,] = epsilon[1,]
  epsilonst = rnorm(nsimu)
  st = rep(0, nsimu)
  st[1] = epsilonst[1]
  for(t in 2:nsimu){
    st[t] = 0.95 * st[(t-1)] + epsilonst[t]
  }
  check = matrix(0, ncol = (m-1), nrow = nsimu)
  for(l in 1:(m-1)){
    check[(st<c[l]),l] = 1
  }
  check1 = rep(0, nsimu)
  if(m == 3){
    check1[(check[,1] != check[,2])] = 1
  }
  checkoppo = +(!check)
  for(i in 2:nsimu){
    if(m == 2){
      y[i,] =  (phi[[1]]%*% as.matrix(y[(i-1),]))*check[i,1] + (phi[[2]]%*% as.matrix(y[(i-1),]))*checkoppo[i,1] + epsilon[i,]
    } else if(m == 3){
      y[i,] =  (phi[[1]]%*% as.matrix(y[(i-1),]))*check[i,1] + 
        (phi[[2]]%*% as.matrix(y[(i-1),]))*check1[i] + (phi[[3]]%*% as.matrix(y[(i-1),]))*checkoppo[i,2] + epsilon[i,]
    }
  }
  return(list(sim = tail(y,nsim), st = tail(st, nsim), residuals = tail(epsilon, nsim)))
}

TVAR.sim <- function(n, nsim, m, phi, c){
  In = diag(n)
  nsimu = nsim + 100
  y = matrix(0, nrow = nsimu, ncol = n)
  epsilon = MASS::mvrnorm(nsimu, mu = rep(0, n), Sigma = In)
  y[1, ] = epsilon[1, ]
  
  st = rep(0, nsimu)
  st[1] = y[1, 1]
  
  check = matrix(0, nrow = nsimu, ncol = m - 1)
  if (m == 3) {
    check1 = rep(0, nsimu)
  }
  
  for (i in 2:nsimu) {
    st[i] = y[i - 1, 1]
    
    for (l in 1:(m - 1)) {
      check[i, l] = as.integer(st[i] < c[l])
    }
    
    if (m == 2) {
      check_oppo = 1 - check[i, 1]
      y[i, ] = (phi[[1]] %*% y[i - 1, ]) * check[i, 1] +
        (phi[[2]] %*% y[i - 1, ]) * check_oppo + epsilon[i, ]
    } else if (m == 3) {
      check_oppo = 1 - check[i, 2]
      check1[i] = as.integer(check[i, 1] != check[i, 2])
      y[i, ] = (phi[[1]] %*% y[i - 1, ]) * check[i, 1] +
        (phi[[2]] %*% y[i - 1, ]) * check1[i] +
        (phi[[3]] %*% y[i - 1, ]) * check_oppo + epsilon[i, ]
    }
  }
  
  return(list(sim = tail(y, nsim), st = tail(st, nsim), residuals = tail(epsilon, nsim)))
}

VLSTAR.sim <- function(n, nsim, m, B, c, gamma, const = FALSE){
  In = diag(n)
  nsimu = nsim+100
  y = matrix(0, nrow = nsimu, ncol = n)
  epsilon = MASS::mvrnorm(nsimu, mu = rep(0, n), Sigma = In)
  
  # Initializing y[1,] and st[1]
  y[1,] = epsilon[1,]
  st = rep(0, nsimu)
  st[1] = y[1,1]
  
  Psit = list()
  Psitpre = list()
  Itmp = diag(0.1, ncol = n, nrow = n)
  Psit[[1]] = rbind(In, matrix(rep(Itmp, (m-1)), ncol = n, nrow = n*(m-1), byrow = TRUE))
  Psitpre[[1]] <- rbind(In, matrix(rep(Itmp, (m-2)), ncol = n, nrow = n*(m-2), byrow = TRUE))
  
  for(i in 2:nsimu){
    # Updating st
    st[i] = y[(i-1), 1]
    
    # Computing Gtilde on the new st
    Gtilde = In
    for(l in 1:(m-1)){
      G <- 1/(1 + exp(-gamma[l]*(st[i]-c[l])))
      Gt = diag(G, ncol = n, nrow = n)
      Gtilde = rbind(Gtilde, Gt)
    }
    
    # Updating y
    if(const == FALSE){
      y[i,] =  t(Gtilde) %*% t(B)%*%y[(i-1),] + epsilon[i,]
    }else if(const == TRUE){
      y[i,] =  t(Gtilde) %*% t(B)%*%c(1,y[(i-1),]) + epsilon[i,]
    }
    
    Psit[[i]] <- Gtilde
    Psitpre[[i]] <- Gtilde[1:(n*(m-1)),]
    }

  return(list(sim = tail(y, nsim), Psit = tail(Psit, nsim), 
              Psitpre = tail(Psitpre, nsim), st = tail(st, nsim), residuals = tail(epsilon, nsim)))
}

VLSTAR.simexo <- function(n, nsim, m, B, c, gamma, rho = 0.95){
  In = diag(n)
  y = matrix(0, nrow = nsim, ncol = n)
  epsilon = mvrnorm(nsim, mu = rep(0, n), Sigma = In)
  y[1,] = epsilon[1,]
  epsilonst = rnorm(nsim)
  st = rep(0, nsim)
  st[1] = epsilonst[1]
  for(t in 2:nsim){
    st[t] = rho * st[(t-1)] + epsilonst[t]
  }
  Psit = list()
  Psitpre = list()
  Itmp = diag(0.1, ncol = n, nrow = n)
  Psit[[1]] = rbind(In, matrix(rep(Itmp, (m-1)), ncol = n, nrow = n*(m-1), byrow = TRUE))
  Psitpre[[1]] <- rbind(In, matrix(rep(Itmp, (m-2)), ncol = n, nrow = n*(m-2), byrow = TRUE))
  for(i in 2:nsim){
    Gtilde = In
    for(l in 1:(m-1)){
      G <- 1/(1 + exp(-gamma[l]*(st[(i-1)]-c[l])))
      Gt = diag(G, ncol = n, nrow = n)
      Gtilde = rbind(Gtilde, Gt)
    }
    y[i,] =  t(Gtilde) %*% t(B)%*%y[(i-1),] + epsilon[i,]
    Psit[[i]] <- Gtilde
    Psitpre[[i]] <- Gtilde[1:(n*(m-1)),]
  }
  return(list(sim = y, Psit = Psit, Psitpre = Psitpre, st = st, residuals = epsilon))
}


lognreg = function(ti, data){
  mY = data$mY        # Matrice delle variabili dipendenti
  X1 = as.matrix(data$X1)  # Matrice delle variabili indipendenti
  st = data$st        # Variabile di stato
  T1 = nrow(mY)       # Numero di osservazioni
  n = ncol(mY)        # Numero di variabili dipendenti
  nx = ncol(X1)
  m = data$m          # Numero di regimi
  
  gamma_start = 1
  gamma_end = n*(m-1)
  c_start = gamma_end + 1
  c_end = n*(m-1)*2
  beta_start = c_end + 1
  beta_end = beta_start-1+nx*n
  psi_start = beta_end + 1
  psi_end = psi_start-1+nx*n*(m-1)
  
  gammalist = split_vector(ti[gamma_start:gamma_end], n)
  clist = split_vector(ti[c_start:c_end], n)
  psilist = split_to_matrices(ti[psi_start:psi_end], n, nx)
  beta = matrix(ti[beta_start:beta_end], ncol = nx, nrow = n, byrow = T)
  
  e1 = matrix(0, nrow = T1, ncol = n)  # Matrice degli errori
  
  for(t in 1:T1){
    reg_effect = rep(0, n)
    for(k in 1:(m-1)){
      gt = matrix(0, nrow = n, ncol = n)
      for(i in 1:n){
        gt[i,i] = 1/(1 + exp(-gammalist[[k]][i] * (st[t] - clist[[k]][i])))
      }
      reg_effect = gt%*%psilist[[k]] %*% X1[t,] + reg_effect
    }
    e1[t,] = as.matrix(mY[t,] - beta%*%X1[t,] - reg_effect)
  }
  
  # Calcolo della matrice di covarianza degli errori
  o1 = t(e1) %*% e1 / T1
  io1 = MASS::ginv(o1)
  
  log_det_o1 = log(det(o1))
  
  # Calcolo del termine quadratico
  quadratic_form = sum((io1 %*% t(e1))*t(e1))
  
  # Calcolo della log-verosimiglianza
  log_likelihood = -T1 * n * 0.5 * log(2 * pi) - 0.5 * T1 * log_det_o1 - 0.5 * quadratic_form
  
  return(-log_likelihood)  # Negativo per la minimizzazione
}



my.ts.panel <- function(x, col = col, bg = bg, pch = pch, type = type,  vpos=8.75, ...){
  lines(x, col = col, bg = bg, pch = pch, type = type, ...)
  abline(v=which(object$st>object$c), col = 'grey')}

jsr <- function(A, B, max_prod_length = 5, norm_type = "spectral") {
  # Check if matrices are square and of the same size
  if (!all(dim(A) == dim(B)) || nrow(A) != ncol(A)) {
    stop("Matrices A and B must be square and of the same dimensions.")
  }
  
  n <- nrow(A)
  
  # List to store all products
  products <- list()
  
  # Initialize with individual matrices
  products[[1]] <- A
  products[[2]] <- B
  
  # Generate all possible products up to max_prod_length
  for (k in 2:max_prod_length) {
    new_products <- list()
    for (prod in products) {
      new_products <- c(new_products, list(prod %*% A, prod %*% B))
    }
    products <- c(products, new_products)
  }
  
  # Function to compute the norm of a matrix
  compute_norm <- function(mat) {
    if (norm_type == "spectral") {
      return(max(abs(eigen(mat, only.values = TRUE)$values)))
    } else {
      return(norm(mat, type = norm_type))
    }
  }
  
  # Compute norms of all products
  norms <- sapply(products, compute_norm)
  
  # Estimate the joint spectral radius
  jsr_estimate <- max(norms)^(1 / max_prod_length)
  
  return(jsr_estimate)
}


euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


multiGW <- function (l1, l2, T, tau, method = c("HAC", "NeweyWest", "Andrews", 
                                                "LumleyHeagerty"), alternative = c("two.sided", "less", 
                                                                                   "greater")) 
{
  if (is.matrix(l1) && ncol(l1) > 2) 
    stop("multivariate time series not allowed")
  if (is.matrix(l2) && ncol(l2) > 2) 
    stop("multivariate time series not allowed")
  if (tau < 1) 
    stop("Predictive Horizon must to be a positive integer")
  if (length(l1) != length(l2)) 
    stop("size of l1 and l2 difier")
  alternative <- match.arg(alternative)
  #DNAME <- deparse(substitute(x))
  dif = l1 - l2
  q = length(dif)
  m = T - q
  n = T - tau - m + 1
  delta = mean(dif)
  mod <- lm(dif ~ 0 + rep(1, q))
  if (tau == 1) {
    re = summary(mod)
    STATISTIC = re$coefficients[1, 3]
    if (alternative == "two.sided") 
      PVAL <- 2 * pnorm(-abs(STATISTIC))
    else if (alternative == "less") 
      PVAL <- round(pnorm(STATISTIC), 4)
    else if (alternative == "greater") 
      PVAL <- round(pnorm(STATISTIC, lower.tail = FALSE), 
                    4)
    names(STATISTIC) <- "Normal Standad"
    METHOD <- "Standard Statistic Simple Regression Estimator"
  }
  if (tau > 1) {
    if (method == "HAC") {
      METHOD <- "HAC Covariance matrix Estimation"
      ds = sqrt(vcovHAC(mod)[1, 1])
    }
    if (method == "NeweyWest") {
      METHOD <- "Newey-West HAC Covariance matrix Estimation"
      ds = sqrt(NeweyWest(mod, tau)[1, 1])
    }
    if (method == "LumleyHeagerty") {
      METHOD <- "Lumley HAC Covariance matrix Estimation"
      ds = sqrt(weave(mod)[1, 1])
    }
    if (method == "Andrews") {
      METHOD <- "kernel-based HAC Covariance matrix Estimator"
      ds = sqrt(kernHAC(mod)[1, 1])
    }
    STATISTIC = delta/ds
    if (alternative == "two.sided") 
      PVAL <- 2 * pnorm(-abs(STATISTIC))
    else if (alternative == "less") 
      PVAL <- pnorm(STATISTIC)
    else if (alternative == "greater") 
      PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
    names(STATISTIC) <- "Normal Standar"
  }
  structure(list(statistic = STATISTIC, alternative = alternative, 
                 p.value = PVAL, method = METHOD))
}

