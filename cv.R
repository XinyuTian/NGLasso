cv <- function(dat, k = 10, type = "brier", lambdaseq=NULL, lambda2seq=0.1, fitype = NULL) {
  
  if(k < 3) stop("'less than 3'-fold crossvalidation not supported")
  nobs <- nrow(dat$y)
  P <- ncol(dat$x) - 1
  Q <- ncol(dat$y)
  if(is.null(lambdaseq)) {
    lambdamax <- getlambdamax(dat=dat)
    lambdaseq <- getlambdaseq(lambdamax)
  }
  n1 <- length(lambdaseq)
  n2 <- length(lambda2seq)
  coef.init <- matrix(1,nrow=Q,ncol=P+1)
  ## a small function giving the norm of difference
  Brier <- function(y, mu, weights = rep(1, nrow(y)), ...){
    weighted.mean(rowMeans((y - mu)^2), weights)
  }
  
  ## permutation until the segmentation of samples are eligible
  sane.permutation <- 0
  while(sane.permutation < 1){
    permutation <- sample(seq(nobs))
    boundaries <- numeric(k-1)
    for(i in seq(1, length(boundaries))){
      boundaries[i] <- ceiling(nobs * i / k)
    }
    boundaries <- c(0, boundaries, nobs)
    cvinds <- list(); length(cvinds) <- k
    for(i in seq(k)){
      cvinds[[i]] <- permutation[seq((boundaries[i] + 1), boundaries[i+1])]
    }
    
    ## test if all columns of y have samples in each training set
    sane.flag <- numeric(k)
    for(i in seq(k)){
      daty.test <- dat$y[-cvinds[[i]], ]
      daty.test <- cbind(daty.test, 1 - rowSums(daty.test))
      if(all(colSums(daty.test) >= 1)) sane.flag[i] <- T
    }
    if(all(sane.flag == T)) sane.permutation <- 1  else sane.permutation <- sane.permutation-1
    if(sane.permutation < -20) stop("data cannot be divided")
  }
  
  newdata <- list(); length(newdata) <- k
  for(j in seq(k))  newdata[[j]] <- list(y = dat$y[cvinds[[j]], ], x = dat$x[cvinds[[j]], ])
  
#  cvlist <- list(); length(cvlist) <- n2
  dfmax <- list(); length(dfmax) <- n2
  dfmin <- list(); length(dfmin) <- n2
  mdev <- list(); length(mdev) <- n2
  mbrierscore <- list(); length(mbrierscore) <- n2
  sdev <- list(); length(sdev) <- n2
  sbrierscore <- list(); length(sbrierscore) <- n2
  
  for (l in 1:n2) {
    lambda2 <- lambda2seq[l]
    # initialization
    i <- 1
    lambda1 <- lambdaseq[1]
    df <- 0
#    cvlist[[l]] <- list()
    dev <-  matrix(nrow = n1, ncol = k)
    brierscore <-  matrix(nrow = n1, ncol = k)
    
    # initial values
    ctl <- 0
    oldmean <- 100
    while(i <= n1) {      
      ## a function doing the j-th cv, returns a list of coef with length = k
      cvcore <- function(j){
        cvdat <- list(y = dat$y[-cvinds[[j]], , drop = F],
                      x = dat$x[-cvinds[[j]], , drop = F],
                      Lmatrix = dat$Lmatrix, pwt = dat$pwt)
        coef <- fista(dat=cvdat, tuning=list(lambda1,lambda2), coef.init = coef.init,
                      fitype = fitype)$coef

        if (all(!is.na(coef))) coef.init <<- coef
        return(coef)  
      }
      ## now the fitting of the cv'd models:
      cvtemp <- lapply(seq(k), cvcore)
      
      ## df
#      dfseq <- lapply(cvtemp, function(coef) sum(coef[1,]!=0)-1)
#      dfseq <- unlist(unlist(dfseq))
#      dfmax[[l]][i] <- max(dfseq)
#      dfmin[[l]][i] <- min(dfseq)
#      df <- max(dfseq)
            
      logl <- c()

      for(j in seq(k)){
        nobsj <- length(cvinds[[j]])
        coefj <- cvtemp[[j]]
        if (all(!is.na(coefj))) {
          etaj <- update.eta(dat = newdata[[j]], coef = coefj, weights = rep(1,nobsj))
          muj <- update.mu(etaj)
          
          logl[j] <- loglik(y = newdata[[j]]$y, mu = muj, weights = rep(1,nobsj))
          dev[i, j] <- 2*(loglik(y = newdata[[j]]$y, mu = newdata[[j]]$y,
                                 weights = rep(1,nobsj)) - logl[j])
          brierscore[i, j] <- Brier(y = newdata[[j]]$y, mu = muj, weights = rep(1,nobsj))
        } else {dev[i, j] <- NA; brierscore[i, j] <- NA}
      }
#      print(paste("this is the ", i, "th iteration."))
#      cvlist[[l]][[i]] <- cvtemp
      newmean <- mean(brierscore[i, ], na.rm = TRUE)
      if (!is.na(newmean) & abs(oldmean-newmean) > 1e-6) {
        if ((oldmean-newmean) <= 0.005*newmean) ctl<-ctl+1 else ctl<-0 
        if (ctl > 3) {break; print(paste("break while loop at i = ", i))}
        oldmean <-newmean
      }
      i <- i+1 
      lambda1 <- lambdaseq[i]
    }
    
    dev <- na.omit(dev)
    brierscore <- na.omit(brierscore)
    mdev[[l]] <- rowMeans(dev, na.rm = TRUE)
    mbrierscore[[l]] <- rowMeans(brierscore, na.rm = TRUE)
    sdev[[l]] <- apply(dev, 1, sd, na.rm = TRUE) / sqrt(k - 1)
    sbrierscore[[l]] <- apply(brierscore, 1, sd, na.rm = TRUE) / sqrt(k - 1)
  }
  
  means <- switch(type,
                  "deviance" = mdev,
                  "brier" = mbrierscore
  )
  
  sds <- switch(type,    ## we compute the sd of the estimator for the mean, not the sample sd!!
                "deviance" = sdev,
                "brier" = sbrierscore
  )
  
  return(list(mean = means, sd = sds, type = type, dfmax = dfmax, dfmin = dfmin, 
              lambdaseq = lambdaseq, lambda2seq = lambda2seq, mdev=mdev, sdev=sdev))
  
}


