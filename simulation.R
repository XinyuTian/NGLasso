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

multinom.simdata <- function(nobs, P, K, Amatrix = Amatrix, rho = 0.5, 
                             coef=NULL, cov=NULL, weights = rep(1, nobs)){
  if (is.null(coef)) {
    warning(paste("coefficients are missing, random values is used"))
    #    coef1 <- rnorm(K-1)
    coef1 <- c(2, -1.5, 1)
    coef <- cbind(coef1,coef1,coef1,coef1,matrix(0,nrow=K-1,ncol=8))
    #    coef <- cbind(coef1,coef1,coef1,coef1,coef1,coef1,coef1,coef1,coef1,coef1,matrix(0,nrow=K-1,ncol=60))
    coef <- cbind(rnorm(K-1), coef)
  }
  Q <- nrow(coef)
  if (is.null(cov)) {
    cov <- getcov(rho=rho, P=P)
  }
  X <- mvrnorm(n = nobs, mu=rep(0,P), cov)
  Lmatrix <- getLmat(Amatrix, signed=TRUE, X=X)
  X <- cbind(1,X)
  
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
  dat <- list()
  dat$x <- X
  dat$y <- y
  dat$Amatrix = Amatrix
  dat$Lmatrix = Lmatrix
  dat$coef <- coef
  pwt <- getpenweights(dat)
  dat$pwt <- pwt
  dat$signed.Lmat <- FALSE
  return(dat)
}

crtcoef <- function (coef0 = c(0, 0, 0), coef1 = c(0.4, -0.25, 0.1), P, nz) {
  coef=numeric(0)
  sapply(seq(nz), function (x) coef <<- cbind(coef, coef1))
  sapply(seq(P-nz), function (x) coef <<- cbind(coef, 0))
  coef <- cbind(coef0, coef) 
  return(coef)
}

crtcoef1 <- function (K=4, P, nz) {
  coef <- matrix(0, K-1, nz)
  sapply(seq(K-1), function(i) coef[i,] <<- sample(seq(0.05, 0.5,by=0.05), size=nz, replace=T)*sample(c(-1,1),1))
  sapply(seq(P-nz), function (x) coef <<- cbind(coef, 0))
  coef <- cbind(c(0, 0, 0), coef) 
  return(coef)
}

crtcoef2 <- function (K=4, P, nz) {
  coef <- matrix(sample(seq(-0.5, 0.5,by=0.05), size=(K-1)*nz, replace=T), K-1, nz)
  sapply(seq(P-nz), function (x) coef <<- cbind(coef, 0))
  coef <- cbind(c(0, 0, 0), coef) 
  return(coef)
}

## switch the first nz columns with the following nz/2 columns
switchCoef <- function(coef, P, nz) {
  coef[ , (nz+2):(nz*3/2+1)] <- coef[ , 2:(nz/2+1)]
  coef[ , 2:(nz/2+1)] <- 0
  return(coef)
}

crtAmat <- function(upper.lim, P, nz) {
  lower.lim <- 1-upper.lim
  k <- P / nz
  m.block <- matrix(lower.lim, nrow=nz, ncol=nz)
  m.rep <- repblockMatrixDiagonal(m.block, k)
  m.base <- matrix(runif(P^2, min=0, max=upper.lim), nrow=P)
  m.sum <- m.base + m.rep
  Amatrix <- forceSymmetric(m.sum)
  diag(Amatrix) <- 1
  return(Amatrix)
}

crtOverlapAmat <- function(P, nz, overlap=5) {
  k <- P / nz
  Amatrix <- matrix(0, nrow=P, ncol=P)
  for (i in seq(k)) {
    ind <- ((i-1)*nz+1) : (i*nz+overlap)
    ind <- ind[ind<=P]
    Amatrix[ind, ind] <- 1
  }
  return(Amatrix)
}

crtBlockWiseCov <- function(P, nz, rho) {
  k <- P/nz
  ## !!!!!!!!!!!!! if(!is.integer(k)) stop("P is not a multiple of nz")
  m.block <- matrix(rho, nrow=nz, ncol=nz)
  diag(m.block) <- 1
  cov <- repblockMatrixDiagonal(m.block, k)
  return(cov)
}
