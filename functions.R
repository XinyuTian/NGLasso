## loglik:    a function with arguments "y", "eta" and "weights"
##            implementing the log-likelihood function. weights must be a
##            vector of length nobs.
##            sum(rowSums(y*eta)-log(rowSums(exp(eta))+1))
loglik <- function(y, mu, weights,  ...){
  muc <- cbind(1 - rowSums(mu), mu)
  muc[which(muc <= 1e-6)] <- 1e-8
  muc[which(muc >= 1 - 1e-6)] <- 1 - 1e-8
  yc <- cbind(1 - rowSums(y), y)
  logl <- sum(weights * yc * log(muc))/nrow(y)   
  return(logl)
}
## ploglik:   penalized log likelihood
ploglik = function(y, mu, coef, weights, tuning, Lmatrix){
  lambda2 <- tuning[[2]]
  if(!is.matrix(coef)){
    coef <- matrix(coef, nrow=ncol(y))
  }
  l <- loglik(y, mu, weights) 
  if (lambda2==0) return(l)  else {
    coef <- as.matrix(coef[,-1])
    Lmatrix <- as.matrix(Lmatrix)
    if (ncol(coef)!=ncol(Lmatrix) | ncol(coef)!=nrow(Lmatrix)) {
      warning(paste("error in network-constraint: dimension of Lmatrix doesnt match"))
      return(l)
    }
    l2 <- apply(coef, 1, function(v)  v %*% Lmatrix %*% v)
    return(l - lambda2*sum(l2))
  }
}

## gradient:  a function with arguments "x", "y", "mu" and "weights"
##            implementing the  gradient of the log-likelihood
##            function. returns a q*p matrix that belongs to x. 
##            x must be an nobs*p matrix. y must be a nobs*q matrix.
gradient = function(dat, mu, weights, ...){
  grad <- crossprod((dat$y - mu), weights * dat$x)/nrow(dat$y)
  return(grad)
}
## penalized gradient
pgradient <- function(dat, mu, coef, weights, tuning, Lmatrix){
  lambda2 <- tuning[[2]]
  if(!is.matrix(coef)){
    coef <- matrix(coef, nrow=ncol(dat$y))
  }
  g <- gradient(dat, mu, weights)
  if (lambda2==0) return(g)  else {
    coef <- as.matrix(coef[,-1])
    Lmatrix <- as.matrix(Lmatrix)
    if (ncol(coef)!=ncol(Lmatrix) | ncol(coef)!=nrow(Lmatrix)) {
      warning(paste("error in network-constraint: dimension of Lmatrix doesnt match"))
      return(g)
    }
    g2 <- 2* coef %*% Lmatrix
    g2 <- cbind(0,g2)
    return(g - lambda2*g2)
  }
}

update.eta <- function(dat, coef, weights, iter.count=0)
{
  if(!is.matrix(dat$x) | !is.matrix(coef)){
    #warning(paste("x or coef is not a matrix, iter.count=",iter.count))
    coef <- matrix(coef, nrow=ncol(dat$y))
  } 
  eta <- tcrossprod(dat$x, coef)
  eta
}

## update.mu:   a function with arguments "eta" implementing the inverse link
##            function. eta must be a matrix of size nobs * Q, with Q = K-1
##            for K response categories.
update.mu <- function(eta){
  exp.eta <- exp(eta)
  exp.eta[which(exp.eta==Inf)]=1e6
  mu <- exp.eta/(rowSums(exp.eta)+1)
  return(mu)
}

## l2norm and generic soft thresholding function
tresh <- function(u, lambda, w){
#  w <- sqrt(length(c(u))) * w
  st <- 1 - lambda * w / norm(as.matrix(u), "F")
  u <- max(st, 0) * u
  u
}
## the solution to the optimization problem
fistaProximal <- function(coef, tuning, penweights){
  coef1 <- as.matrix(coef[,-1])
  coef1 <- apply(coef1, 2, tresh, lambda=tuning[[1]], w=penweights)
  return(cbind(coef[,1], coef1))
}

## penalty lambda1*J1
penalty <- function(coef, tuning, Q, penweights){
  if(!is.matrix(coef)){
    coef <- matrix(coef, nrow=Q)
  }
  coef <- as.matrix(coef[,-1])
  coef <- sweep(coef, 2, penweights, `*`)
  J1 <- sum(apply(as.matrix(coef), 2, function(u) norm(as.matrix(u),"F")))
  return(tuning[[1]] * J1)
}


## obj objective function: penalized likelihood, which we want to minimize
obj <- function(coef1, dat, weights, tuning, penweights){
  y <- dat$y
#  if(!is.matrix(coef1)){
#    coef1 <- matrix(coef1, nrow=ncol(dat$y))
#  }
  eta <- update.eta(dat, coef1, weights, iter.count=1)
  mu <- update.mu(eta)
  aa <- ploglik(dat$y, mu, coef=coef1, weights=weights, tuning, dat$Lmatrix)
  bb <- 0
  if (tuning[[1]] != 0) bb <- penalty(coef1, tuning, ncol(y), penweights)
  return(-aa + bb)
} 

## calculate the weights to ADAPTIVE LASSO by regression on each variable
getpenweights <- function(dat) {
  ## a function calls multinom to do one-variable regression
  oneml <- function(u) {
    fit<-multinom(cbind(1-rowSums(dat$y), dat$y) ~ u)
    norm(as.matrix(summary(fit)$coefficients[,2]), type="F")
  }
  weights <- apply(dat$x[,-1], 2, oneml)
  weights <- 1/weights
  return (weights)
}

## generate coviance matrix of X
## symmetric, diagonal is rho, geometric 
getcov <- function(rho, P, ctl=1e-6) {
  v1 <- c((P-1):0, 1:(P-1))
  longv <- rho ^ v1
  longv1 <- sapply(longv, function(x) ifelse(x<ctl,0,x))
  sapply(seq(P), function(i) longv1[(P+1-i):(2*P-i)])
}

## the indices of selected variables, degree of freedom
nvar <- function(coef){
  ind <- apply(as.matrix(coef), 2, function(u) any(u!=0))
  return(which(ind==1))
}


## find the mode of a sequence, that is, the most frequent element
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
