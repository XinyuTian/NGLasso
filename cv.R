
cv <- function(dat, k = 10, dfmax=NULL, crt.measure = "brier", lambdaseq=NULL, lambda2seq=0.1, 
               fitype = NULL, ratio=NULL, scl=NULL) {
  
  if(k < 3) stop("'less than 3'-fold crossvalidation not supported")
  nobs <- nrow(dat$y)
  P <- ncol(dat$x) - 1
  Q <- ncol(dat$y)
  if(is.null(dfmax)) dfmax <- P/2
  if(is.null(lambdaseq)) {
    lambdamax <- getlambdamax(dat=dat)
    lambdaseq <- getlambdaseq(lambdamax=lambdamax, ratio=ratio, scl=scl)
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
  
  df <- list(); length(df) <- n2
  means <- list(); length(means) <- n2
  sds <- list(); length(sds) <- n2
  
  for (l in 1:n2) {
    lambda2 <- lambda2seq[l]
    # initialization
    i <- 1
    lambda1 <- lambdaseq[1]
    dftemp <- 0
    oldmean <- 100

    measure <-  matrix(nrow = n1, ncol = k)
    
    while(i <= n1 & dftemp <= dfmax) {      
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
      dfseq <- lapply(cvtemp, function(coef) length(nvar(coef))-1 )
      dfseq <- unlist(unlist(dfseq))
      dftemp <- Mode(dfseq)
      df[[l]][i] <- dftemp
      logl <- c()

      for(j in seq(k)){
        nobsj <- length(cvinds[[j]])
        coefj <- cvtemp[[j]]
        if (all(!is.na(coefj))) {
          etaj <- update.eta(dat = newdata[[j]], coef = coefj, weights = rep(1,nobsj))
          muj <- update.mu(etaj)
          
          if (crt.measure == "brier") 
            measure[i, j] <- Brier(y = newdata[[j]]$y, mu = muj, weights = rep(1,nobsj))
          else if (crt.measure == "dev") {
            logl[j] <- loglik(y = newdata[[j]]$y, mu = muj, weights = rep(1,nobsj))
            measure[i, j] <- 2*(loglik(y = newdata[[j]]$y, mu = newdata[[j]]$y,
                                       weights = rep(1,nobsj)) - logl[j])
          }
        } else measure[i, j] <- NA
      }
      i <- i+1 
      lambda1 <- lambdaseq[i]
    }
    
    measure <- na.omit(measure)
    means[[l]] <- rowMeans(measure, na.rm = TRUE)
    sds[[l]] <- apply(measure, 1, sd, na.rm = TRUE) / sqrt(k)
  }
  
  return(list(mean = means, sd = sds, type = crt.measure, df = df, 
              lambdaseq = lambdaseq, lambda2seq = lambda2seq, fold=k))
  
}


