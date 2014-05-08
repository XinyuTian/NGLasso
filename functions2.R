## find the lambda.min and lambda.1se of a cv.
## if length(lambda2seq) > 1, choose the global minimum
findlambda <- function(dat, lambdaseq = NULL, lambda2seq=NULL, cv1, fitype = NULL) {
#  if (missing(cv1)) cv1 <<- cv(dat=dat, lambdaseq = lambdaseq, lambda2seq = lambda2seq, adapt = adapt)
  lambdaseq <- cv1$lambdaseq
  lambda2seq <- cv1$lambda2seq
  n2 <- length(lambda2seq)
#  if (fitype == "refit" | fitype == "adprf") {
#    meanmin <- lapply(cv1$mean, Mode)
#  } else  meanmin <- lapply(cv1$mean, min)
  meanmin <- lapply(cv1$mean, min)
  indmmin <- which.min(unlist(meanmin))
#  if (fitype == "refit" | fitype == "adprf") {
#    indmin <-  median(which(cv1$mean[[indmmin]] == meanmin[[indmmin]]))
#  } else    indmin <- which.min(cv1$mean[[indmmin]])
  indmin <- which.min(cv1$mean[[indmmin]])
  minsd <- cv1$sd [[indmmin]] [indmin]
  m1se <- meanmin[[indmmin]] + minsd
  lambda.min <- lambdaseq[indmin]
  mdiff <- cv1$mean[[indmmin]] - m1se
  ind1se <- which.min(abs(mdiff))
  lambda.1se <- lambdaseq[ind1se]
  
  return(list(lambda.min = lambda.min, lambda.1se = lambda.1se, lambda2seq=lambda2seq, indmmin = indmmin))
}

## return all coefs for each combination of lambda1 and lambda2 based on fista
## it is called in solution path
findcoef <- function(dat, lambdaseq=NULL, lambda2seq=0.1, fitype = NULL) {
  if(is.null(lambdaseq)) {
    lambdamax <- getlambdamax(dat=dat)
    lambdaseq <- getlambdaseq(lambdamax)
  }
  n1 <- length(lambdaseq)
  n2 <- length(lambda2seq)
  tempcoef <- list()
  for (l in 1:n2) {
    lambda2 <- lambda2seq[l]
    
    ## a function calls fista, returning coef
    pfista <- function(lambda1){
      coef <- fista(dat=dat, tuning=list(lambda1, lambda2), fitype = fitype)$coef
      return(coef)  
    }
    ## now the fitting of the cv'd models:
    tempcoef[[l]] <- sapply(lambdaseq, pfista, simplify = FALSE)
  }
  return(tempcoef)
}


## the refit function is wrong!!
refit <- function(dat, lambda1, lambda2, allcoef=NULL, lambdaseq=NULL, lambda2seq=NULL) {
  
  if (!is.null(allcoef)) {
    ind1 <- integer(0)
    ind2 <- integer(0)
    if (is.null(lambdaseq)) {
      lambdamax <- getlambdamax(dat=dat)
      lambdaseq <- getlambdaseq(lambdamax)
    }
    ind1 <- which(lambdaseq == lambda1) 
    if (!is.null(lambda2seq)) ind2 <- which(lambda2seq == lambda2)
    
    if(length(ind1) + length(ind2) > 0) coef <- allcoef[[ind2]][[ind1]]
  } else {
    coef <- fista(dat=dat, tuning=list(lambda1, lambda2))$coef
  }
  
  ind <- nvar(coef)
  if(ind[1]==1) ind0 <- ind[-1]-1 else ind0 <- ind-1
  newdata <- list(y=dat$y, x=dat$x[, ind], Lmatrix=Lmatrix[ind0,ind0])
  out <-multinom(cbind(1-rowSums(newdata$y), newdata$y) ~ newdata$x)
  return(summary(out))
}

## a fucntion returns the results of model selection 
## using the lambda.min or lambda.1es from findlambda
slcmod <- function (dat, flambda, lambda1seq=NULL, lambda2seq=NULL, 
                    type="1se", fitype = NULL) {
  if (!missing(flambda)) {
    lambda2seq <- flambda$lambda2seq
    lambda1seq <- switch(type,
                         "min" = flambda$lambda.min,
                         "1se" = flambda$lambda.1se
    )
  }
  n2 <- length(lambda2seq)
  lambda1seq <- unlist(lambda1seq)
  coef <- dat$coef
  out <- list()
  for (l in 1:n2) {
    lambda2 <- lambda2seq[l]
    lambda1 <- lambda1seq[l]
    tempcoef <- fista(dat=dat, tuning=list(lambda1,lambda2), fitype=fitype)$coef
    ind <- nvar(tempcoef)
    mse <- sum((coef-tempcoef)^2)/length(coef)
    out[[l]] <- list(ind=ind, coef=tempcoef, MSE=mse)
  }
  return(out)
}

## prediction
pred <- function(newdat, coef0, weights=NULL) {
  y <- newdat$y
  if (is.null(weights)) weights <- rep(1, nrow(y))
  eta <- update.eta(newdat, coef0, weights=weights)
  mu <- update.mu(eta)
  yc <- cbind(y,1-rowSums(y))
  muc <- cbind(mu,1-rowSums(mu))
  yv <- as.vector(yc)
  muv <- as.vector(muc)
  predev <- norm(as.matrix(yv-muv), type="F")
  yi <- apply(yc, 1, which.max)
  mui <- apply(muc, 1, which.max)
  accu <- sum(yi == mui)/nrow(y)
  return(c(predev, accu))
}


## exstract coef from cv.glmnet fit
## returns a list, the first is grouped-penalized coef/mse, the second is unprouped
## type could be "1se" or "min"
excoefg <- function(dat, type="1se") {
  gc <- cv.glmnet(dat$x[,-1], y=cbind(dat$y,1-rowSums(dat$y)), 
                  family="multinomial", standardize=F, type.multinomial="grouped")
  li <- switch(type,
               "min" = which(gc$lambda == gc$lambda.min),
               "1se" = which(gc$lambda == gc$lambda.1se)
  )
  g1 <- gc$glmnet.fit
  coef1 <-lapply(g1$beta, function(u) u[,li])
  coef1 <- unlist(coef1)
  coef1 <- matrix(coef1, nrow=4,byrow=T)
  alpha=g1$a0[,li]
  coef1 <- cbind(alpha, coef1)
  coef1 <- apply(coef1,1, function(u) u-coef1[4,])
  coef1 <- t(coef1[,-4])
  
  g2 <- cv.glmnet(dat$x[,-1], y=cbind(dat$y,1-rowSums(dat$y)), 
                  family="multinomial", standardize=F)
  li2 <- switch(type,
                "min" = which(g2$lambda == g2$lambda.min),
                "1se" = which(g2$lambda == g2$lambda.1se)
  )
  g3 <- g2$glmnet.fit
  coef3 <-lapply(g3$beta, function(u) u[,li2])
  coef3 <- unlist(coef3)
  coef3 <- matrix(coef3, nrow=4,byrow=T)
  alpha2=g3$a0[,li2]
  coef3 <- cbind(alpha2, coef3)
  coef3 <- apply(coef3,1, function(u) u-coef3[4,])
  coef3 <- t(coef3[,-4])
  
  out <- list(coefgrp = coef1, coefugp = coef3)
  return(out)
}

## extract coef from fista
## calls findlambda to get the lambda.min or lambda.1es
## returns a list of length n2, in each element contains (df, coef,) mse
## indmmin indicates which lambda2 to use
excoeff <- function (dat, type="1se", fitype = NULL, lambda2seq=1, dfmax=NULL) {
  if(is.null(fitype)) fitype <- "ordinary"
#  if (fitype == "refit" | fitype == "adprf") type="min"
  cv1 <- cv(dat=dat, k=5, dfmax=dfmax, fitype=fitype, lambda2seq=lambda2seq)
  flambda <- findlambda(dat = dat, cv1=cv1, fitype = fitype)
  lambda1 <- switch(type,
                       "min" = flambda$lambda.min,
                       "1se" = flambda$lambda.1se
  )
  indmmin <- flambda$indmmin
  lambda2 <- lambda2seq[indmmin]
  coef0 <- fista(dat=dat, tuning=list(lambda1,lambda2), fitype=fitype)$coef
  return(coef0)
}

## extract coef from multinom
excoefm <- function(dat) {
  out <- multinom(cbind(1-rowSums(dat$y), dat$y) ~ dat$x[,-1])
  coefm <- summary(out)$coefficients
  return(coefm)
}

## input the type of data to be simulated, simulate a group of data, gives output
outfct <- function (setting="small", lambda2seq, type="1se") {
  coef <- switch(setting,
                 "small" = crtcoef(P=20,nz=4),
                 "medium" = crtcoef(P=100,nz=10),
                 "large" = crtcoef(P=200,nz=10,coef1=c(0.5, -1.2, 1))
    )
  Lmatrix <- switch(setting,
                    "small" = crtLmat(P=20,nz=4),
                    "medium" = crtLmat(P=100,nz=10),
                    "large" = crtLmat(P=200,nz=10)
    )
  
  dat <- switch(setting,
                "small" = multinom.simdata(nobs = 200, P = 20, K = 4, coef = coef, Lmatrix = Lmatrix),
                "medium" = multinom.simdata(nobs = 200, P = 100, K = 4, coef = coef, Lmatrix = Lmatrix),
                "large" = multinom.simdata(nobs = 200, P = 200, K = 4, coef = coef, Lmatrix = Lmatrix)
  )
  predat <- switch(setting,
                "small" = multinom.simdata(nobs = 200, P = 20, K = 4, coef = coef, Lmatrix = Lmatrix),
                "medium" = multinom.simdata(nobs = 200, P = 100, K = 4, coef = coef, Lmatrix = Lmatrix),
                "large" = multinom.simdata(nobs = 200, P = 200, K = 4, coef = coef, Lmatrix = Lmatrix)
  )
  dfmax <- switch(setting,
                "small" = 10,
                "medium" = 50,
                "large" = 50
  )
  coef <- coefsdat(dat, lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  crite <- getcrt(coef=coef, coef0=dat$coef, setting=setting, predat=predat)
  return(crite)
}

## the function generates a list containing 7 coefs
coefsdat <- function(dat, lambda2seq, type=NULL, dfmax=NULL) {
  coef <- list(7)
  coef[[1]] <- excoefm(dat)
  gfit <- excoefg(dat)
  coef[[2]] <- gfit$coefugp
  coef[[3]] <- gfit$coefgrp
  coef[[4]] <- excoeff(dat=dat, lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  coef[[5]] <- excoeff(dat=dat, fitype="adapt", lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  coef[[6]] <- excoeff(dat=dat, fitype="refit", lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  coef[[7]] <- excoeff(dat=dat, fitype="adprf", lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  return(coef)
}

## input a list of coefs and the true value coef0
## generate 6 criteria ---- MSE of coefficients, FNR/FPR/FDR for variable selection,
## Brier score and prediction accuracy
getcrt <- function(coef, coef0, setting, predat) {
  noc <- length(coef0)
  mse <- lapply(coef, function(c1) sum((c1-coef0)^2) / noc) 
  ind <- lapply(coef, nvar)
  ind0 <- switch(setting,
                 "small" = c(2:5),
                 "medium" = c(2:11),
                 "large" = c(2:11)
    )
  P <- switch(setting,
              "small" = 20,
              "medium" = 100,
              "large" = 200
  )
  nz <- length(ind0)
  nTP <- lapply(ind, function(XX) sum(!is.na(match(ind0, XX))))
  nFN <- lapply(nTP, function(XX) nz - XX)
  nFP <- mapply(function(xx, yy) length(xx) - yy - 1, ind, nTP)
  nTN <- lapply(nFP, function(XX) P-nz - XX)
  predcrite <- lapply(coef, pred, newdat=predat)
  predcrite0 <- do.call(rbind, predcrite)
  brier <- predcrite0[,1]
  accu <- predcrite0[,2]
  return(list(MSE = unlist(mse), FPR = as.numeric(nFP)/(P-nz), FNR = as.numeric(nFN)/nz,
              FDR = as.numeric(nFP)/(P-nz-as.numeric(nFN)+as.numeric(nFP)), brier=brier, accuacy=accu))
}

