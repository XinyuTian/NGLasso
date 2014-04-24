blockMatrixDiagonal<-function(...){  
  matrixList<-list(...)  
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
} 

repblockMatrixDiagonal <- function(matrix, rep) {
  MatrixList <- list()
  for (i in 1:rep) MatrixList[[i]] <- matrix
  return (blockMatrixDiagonal(MatrixList))
}

#Amatrix <- repblockMatrixDiagonal(matrix(1,nrow=4,ncol=4), rep=3)
#Amatrix <- repblockMatrixDiagonal(matrix(1,nrow=10,ncol=10), rep=7)
#Dmatrix <- diag(rowSums(Amatrix))
#Lmatrix <- Dmatrix -Amatrix
#di <- 1/sqrt(diag(Lmatrix))
#Lmatrix <- t(t(Lmatrix*di)*di)
#Lmatrix <- Lmatrix %*% diag(1/diag(Lmatrix))

multinom.simdata <- function(nobs, P, K, Lmatrix = Lmatrix, 
                             coef=NULL, weights = rep(1, nobs)){
#  source(file="~/My R Files/NGLasso/functions.R")
  if (is.null(coef)) {
    warning(paste("coefficients are missing, random values is used"))
#    coef1 <- rnorm(K-1)
    coef1 <- c(2, -1.5, 1)
    coef <- cbind(coef1,coef1,coef1,coef1,matrix(0,nrow=K-1,ncol=8))
#    coef <- cbind(coef1,coef1,coef1,coef1,coef1,coef1,coef1,coef1,coef1,coef1,matrix(0,nrow=K-1,ncol=60))
    coef <- cbind(rnorm(K-1), coef)
  }
  Q <- nrow(coef)
  X <- matrix(rnorm(nobs*P,mean=0), nrow=nobs, ncol=P)
  X <- cbind(1,X)
#  if(missing(Lmatrix)) Lmatrix <- diag(P); warning(paste("Lmatrix is missing, identity is used"))
  
  eta <- tcrossprod(X, coef)
  mu <- exp(eta)/(rowSums(exp(eta))+1)
  muk <- cbind(mu, 1-rowSums(mu))
  muk[muk < 0] <- 0
  
  sane.permutation <- 0
  while(sane.permutation < 1){
    if (all(weights == 1)) {
      y <- apply(muk, 1, rmultinom, n = 1, size = 1)
      y <- t(y)
    } else {
      y <- matrix(nrow=nobs, ncol=K)
      for(i in seq(nobs)){
        y[i,] <- rmultinom(1, size=weights[i], prob = muk[i,])
      }
      y <- sweep(y, 1, weights, FUN = "/")
    }
    ysum <- colSums(y)
    if(all(ysum != 0)) sane.permutation <- 1    
  }  
  y <- y[,-K]
  ## column centered data
  #X <- apply(X, 2, function(u) scale(u, center=T, scale=F))
  #y <- apply(y, 2, function(u) scale(u, center=T, scale=F))
  dat <- list()
  dat$x <- X
  dat$y <- y
  dat$Lmatrix = Lmatrix
  dat$coef <- coef
  pwt <- getpenweights(dat)
  dat$pwt <- pwt
  return(dat)
}

crtcoef <- function (coef0 = c(0.5, 1.0, -2), coef1 = c(2, -1.5, 1), P, nz) {
  coef=numeric(0)
  sapply(seq(nz), function (x) coef <<- cbind(coef, coef1))
  sapply(seq(P-nz), function (x) coef <<- cbind(coef, 0))
  coef <- cbind(coef0, coef) 
  return(coef)
}

crtLmat <-  function (P, nz) {
  k <- P/nz
  ## !!!!!!!!!!!!! if(!is.integer(k)) stop("P is not a multiple of nz")
  Amatrix <- repblockMatrixDiagonal(matrix(1,nrow=nz,ncol=nz), rep=k)
  #Amatrix <- repblockMatrixDiagonal(matrix(1,nrow=10,ncol=10), rep=7)
  Dmatrix <- diag(rowSums(Amatrix))
  Lmatrix <- Dmatrix -Amatrix
  di <- 1/sqrt(diag(Lmatrix))
  Lmatrix <- t(t(Lmatrix*di)*di)  
  return(Lmatrix)
}

crtcoef1 <- function (K=4, P, nz) {
  coef <- matrix(0, K-1, nz)
  sapply(seq(K-1), function(i) coef[i,] <<- sample(seq(0.4, 2,by=0.2), size=nz, replace=T)*sample(c(-1,1),1))
  sapply(seq(P-nz), function (x) coef <<- cbind(coef, 0))
  coef <- cbind(sample(seq(-2, 2,by=0.2), size=K-1, replace=T), coef) 
  return(coef)
}
