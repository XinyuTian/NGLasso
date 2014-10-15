## find the lambda.min and lambda.1se of a cv.
## if length(lambda2seq) > 1, choose the global minimum
findlambda <- function(dat, lambdaseq = NULL, lambda2seq=NULL, cv1, fitype = NULL) {
  #  if (missing(cv1)) cv1 <<- cv(dat=dat, lambdaseq = lambdaseq, lambda2seq = lambda2seq, adapt = adapt)
  lambdaseq <- cv1$lambdaseq
  lambda2seq <- cv1$lambda2seq
  n2 <- length(lambda2seq)
  meanmin <- lapply(cv1$mean, min)
  indl2 <- which.min(unlist(meanmin))
  cvm <- cv1$mean[[indl2]]
  indmin <- which.min(cvm)
  minsd <- cv1$sd [[indl2]] [indmin]
  mini <- meanmin[[indl2]]
  all.1se <- which(cvm < (mini + minsd))
  lambda.min <- lambdaseq[indmin]
  lambda.1se <- lambdaseq[min(all.1se)]
  
  return(list(lambda.min = lambda.min, lambda.1se = lambda.1se, 
              lambda2seq=lambda2seq, indmmin = indl2))
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
  predev <- norm(as.matrix(yc-muc), type="F")
  yi <- apply(yc, 1, which.max)
  mui <- apply(muc, 1, which.max)
  accu <- sum(yi == mui)/nrow(y)
  return(c(predev, accu))
}



## exstract coef using cv.glmnet, "grouped" or "ungrouped"
## returns a list, type=
## "both"  list: 1-coef.min, 2-coef.1se
## "min"   list: 1-coef.min
## "1se"   list: 1-coef.1se
excoefg <- function(dat, type=c("both", "min", "1se"), is.group="grouped") {
  gc <- cv.glmnet(dat$x[,-1], y=cbind(dat$y,1-rowSums(dat$y)), 
                  family="multinomial", standardize=F, type.multinomial=is.group)
  
  if (type=="both") return(list(coef.min=excoefg.coef(gc,type="min"), 
                                coef.1se=excoefg.coef(gc,type="1se")))
  else {
    out <- list(excoefg.coef(gc=gc,type=type))
    names(out) <- paste0("coef.", type)
    return(out)
  }
}

# a function nested inside excoefg
# extract coef from the result of cv.glmnet, "1se", "min" or "both"
excoefg.coef <- function(gc, type){
  li <- switch(type,
               "min" = which(gc$lambda == gc$lambda.min),
               "1se" = which(gc$lambda == gc$lambda.1se)
  )
  g1 <- gc$glmnet.fit
  Q <- length(g1$beta)
  coef1 <-lapply(g1$beta, function(u) u[,li])
  coef1 <- unlist(coef1)
  coef1 <- matrix(coef1, nrow=Q,byrow=T)
  alpha=g1$a0[,li]
  coef1 <- cbind(alpha, coef1)
  coef1 <- apply(coef1,1, function(u) u-coef1[Q,])
  coef1 <- t(coef1[,-Q])
  return(coef1)
}

## extract coef from cross validation based on fista
## calls findlambda to get the lambda.min or lambda.1es
## returns a list, type=
## "both"  list: 1-coef.min, 2-coef.1se, 3-ind
## "min"   list: 1-coef.min, 2-ind
## "1se"   list: 1-coef.1se, 2-ind
## indmmin indicates which lambda2 to use
excoeff <- function (dat, type="both", fitype = NULL, lambda2seq=1, dfmax=50) {
  if(is.null(fitype)) fitype <- "ordinary"
  cv1 <- cv(dat=dat, k=5, dfmax=dfmax, fitype=fitype, lambda2seq=lambda2seq)
  flambda <- findlambda(dat = dat, cv1=cv1, fitype = fitype)
  indmmin <- flambda$indmmin
  lambda2 <- lambda2seq[indmmin]
  coef.min <- fista(dat=dat, tuning=list(flambda$lambda.min,lambda2), fitype=fitype)$coef
  coef.1se <- fista(dat=dat, tuning=list(flambda$lambda.1se,lambda2), fitype=fitype)$coef
  
  if (type=="both") return(list(coef.min=coef.min, coef.1se=coef.1se, ind=indmmin))
  else if (type=="min") return(list(coef.min=coef.min, ind=indmmin))
  else if (type=="1se") return(list(coef.1se=coef.1se, ind=indmmin))
  else stop("'type' should be one of 'both', '1se', 'min'")
}

## extract coef from multinom
excoefm <- function(dat) {
  out <- multinom(cbind(1-rowSums(dat$y), dat$y) ~ dat$x[,-1])
  coefm <- summary(out)$coefficients
  return(coefm)
}

## input the type of data to be simulated, simulate a group of data, gives output
outfct <- function (modsize="small", coeftype="ideal", Lmattype="ideal", SNR=1, 
                    covtype="blockwise", rho=0.5, lambda2seq, type="1se") {
  P <- switch(modsize,
                  "small" = 20,
                  "medium" = 100,
                  "large" = 200
  )
  nz <- switch(modsize,
                  "small" = 4,
                  "medium" = 10,
                  "large" = 10
  )
  coef <- switch(coeftype,
                 "ideal" = crtcoef(P=P,nz=nz),
                 "simil" = crtcoef1(P=P,nz=nz),
                 "randm" = crtcoef2(P=P,nz=nz)
  )
  coef <- coef*SNR
  cov <- switch(covtype,
                 "AR" = getcov(rho=rho, P=P),
                 "blockwise" = crtBlockWiseCov(P=P,nz=nz, rho=rho)
  )
  Amatrix <- switch(Lmattype,
                    "ideal" = crtAmat(upper.lim=0, P=P,nz=nz),
                    "noise" = crtAmat(upper.lim=0.25, P=P,nz=nz),
                    "incor" = crtAmat(upper.lim=0.6, P=P,nz=nz)
  )
  dat <- multinom.simdata(nobs = 200, P = P, K = 4, coef = coef, cov=cov, Amatrix = Amatrix)
  predat <- multinom.simdata(nobs = 200, P = P, K = 4, coef = coef, cov=cov, Amatrix = Amatrix)
  dfmax <- switch(modsize,
                  "small" = 10,
                  "medium" = 50,
                  "large" = 50
  )
  coefs <- coefsdat(dat, lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  coefs.min <- coefs$coefs.min
  coefs.1se <- coefs$coefs.1se
  ind <- coefs$ind
  crt.min <- getcrt(coef=coefs.min, coef0=dat$coef, modsize=modsize, predat=predat)
  crt.1se <- getcrt(coef=coefs.1se, coef0=dat$coef, modsize=modsize, predat=predat)
  al.crt <- list(crt.min=crt.min, crt.1se=crt.1se, ind=ind)
  return(list(criteria = al.crt, coefs = coefs, dat = dat))
}

## the function generates a list containing 7 coefs
coefsdat <- function(dat, lambda2seq, type=NULL, dfmax=NULL) {
  coefs.min <- list(7)
  coefs.1se <- list(7)
  coefs.min[[1]] <- excoefm(dat)
  coefs.1se[[1]] <- coefs.min[[1]]
  gfit.ung <- excoefg(dat, type=type, is.group="ungrouped")
  coefs.min[[2]] <- gfit.ung$coef.min
  coefs.1se[[2]] <- gfit.ung$coef.1se
  gfit.grp <- excoefg(dat, type=type, is.group="grouped")
  coefs.min[[3]] <- gfit.grp$coef.min
  coefs.1se[[3]] <- gfit.grp$coef.1se
  fit1 <- excoeff(dat=dat, lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  fit2 <- excoeff(dat=dat, fitype="adapt", lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  fit3 <- excoeff(dat=dat, fitype="refit", lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  fit4 <- excoeff(dat=dat, fitype="adprf", lambda2seq=lambda2seq, type=type, dfmax=dfmax)
  coefs.min[[4]] <- fit1$coef.min
  coefs.1se[[4]] <- fit1$coef.1se
  coefs.min[[5]] <- fit2$coef.min
  coefs.1se[[5]] <- fit2$coef.1se
  coefs.min[[6]] <- fit3$coef.min
  coefs.1se[[6]] <- fit3$coef.1se
  coefs.min[[7]] <- fit4$coef.min
  coefs.1se[[7]] <- fit4$coef.1se
  ind <- c(fit1$ind, fit2$ind, fit3$ind, fit4$ind)
  return(list(coefs.min=coefs.min, coefs.1se=coefs.1se, ind=ind))  
}

## input a list of coefs and the true value coef0
## generate 6 criteria ---- MSE of coefficients, FNR/FPR/FDR for variable selection,
## Brier score and prediction accuracy
getcrt <- function(coef, coef0, modsize, predat) {
  nparam <- length(coef0)
  mse <- lapply(coef, function(c1) sum((c1-coef0)^2) / nparam) 
  ind <- lapply(coef, nvar)
  ind0 <- switch(modsize,
                 "small" = c(2:5),
                 "medium" = c(2:11),
                 "large" = c(2:11)
  )
  P <- switch(modsize,
              "small" = 20,
              "medium" = 100,
              "large" = 200
  )
  nz <- length(ind0)
  nTP <- as.numeric(lapply(ind, function(XX) sum(!is.na(match(ind0, XX)))))
  nFN <- nz - nTP
  nFP <- mapply(function(xx, yy) length(xx) - yy - 1, ind, nTP)
  nTN <- P-nz - nFP
  predcrite <- lapply(coef, pred, newdat=predat)
  predcrite0 <- do.call(rbind, predcrite)
  brier <- predcrite0[,1]
  accu <- predcrite0[,2]
  return(list(MSE = unlist(mse), TP = nTP, FP = nFP, TN = nTN, FN = nFN, brier=brier, accuacy=accu))
}
