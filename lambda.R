getlambdamax <- function(dat, weights = NULL, threshold = F, refit = F)
{
  nobs <- nrow(dat$y)
  Q <- ncol(dat$y)
  P <- ncol(dat$x)
  
  coef <- matrix(0, nrow = Q, ncol = P+1)
  if (is.null(weights)) weights <- rep(1, nobs)
  ybar <- apply(as.matrix(dat$y), 2, function(u) weighted.mean(u, w = weights))
  if(any(c(ybar, 1-sum(ybar)) == 0)){
    stop("ill-conditioned data: at least one of the response categories was not observed")
  } 
  #warning(paste("ill-conditioned data: at least one of the response categories was not observed"))  
  ## if there are ties in the relative frequencies of the categories, break them
  if(length(unique(ybar)) != length(ybar)){
    ybar <- c(1 - sum(ybar), ybar) + runif(length(ybar) + 1, 0, 1e-2)
    lyb <- sum(ybar)
    ybar <- ybar / lyb
    ybar <- ybar[-1]
  } 
  coef[,1] <- log(ybar/(1-sum(ybar)))
  
  ## handling the intercept
  which.intercept <- which(apply(dat$x, 2, function(u){all(u==1)}) == T)
  if(length(which.intercept) == 0)
    stop("intercept is missing.")
  if((length(which.intercept) == 1) && which.intercept != 1){
    dat$x[,c(1, which.intercept)] <- dat$x[,c(which.intercept, 1)] 
    if(!all(coef == 0))
      coef[,c(1,which.intercept)] <- coef[,c(which.intercept,1)] 
  }
  if(length(which.intercept) > 1)
    stop("more than one intercept column specified")
  
  ## remove all penalized predictors from the model 
  dat.unpen <- dat
  dat.unpen$x <- matrix(dat$x[,1], nrow=nobs)
  coef.unpen <- coef[,1]
  
  tuning <- list(lambda1 = 0, lambda2= 0)
  
  ## now the call to fista where the core fitting takes place.
  updated <- fista(coef.init=coef.unpen, dat=dat.unpen, weights=weights, tuning=tuning)
  mu0 <- updated$mu
  
  ## the gradient based on the model where all the penalized parameters are 0
  grad0 <- gradient(dat = dat, mu = mu0, weights = weights)[,-1]
  
  ## now the norms of the gradient of the penalized groups, based on the null model.
  gradnorms.x <- apply(grad0, 2, function(v) {norm(as.matrix(v), type="f")})
  lambdamax <- max(gradnorms.x)
  
  return(lambdamax)
}



getlambdaseq <- function(lambdamax, ratio=0.95, scl=0.05) { 
  x <- lambdamax
  res <- integer() 
  while(x >= scl * lambdamax) { 
    res <- c(res, x) 
    x <- x * ratio
  } 
  res 
}



