fista <- function(dat, weights=rep(1,nrow(dat$y)), tuning, coef.init=NULL, fitype = NULL, 
                  stepsize.init = 4,  max.iter = 500, rel.tol = 1e-6, worsen.max=20)
{
#  sancheck <- model@sancheck
#  environment(sancheck) <- sys.frame(sys.nframe())
#  do.sancheck <- TRUE
  nobs <- nrow(dat$y)
  P <- ncol(dat$x)-1
  Q <- ncol(dat$y)
  K <- Q+1
  iter.count <- 0
  Lmatrix <- dat$Lmatrix
  if(is.null(fitype)) fitype <- "ordinary"
  if(fitype!="ordinary" & fitype!="adapt" & fitype!="refit" & fitype!="adprf") 
  stop("'fitype' should be one of 'ordinary', 'adapt', 'refit', 'adprf'")

  ## initializing
  if (is.null(coef.init)) coef.init <- matrix(1,nrow=K-1,ncol=P+1)
  coef.old1 <- coef.init
  coef <- coef.init
  ## internally, coef must be a list, therefore:
#  if(!is.list(coef.init)){
#    coef.old1 <- list(coef.old1)
#    coef <- list(coef)
#  }
#  eta <- update.eta(dat=dat, coef=coef, weights=weights, iter.count)
  eta <- matrix(1, nrow=nobs, ncol=K-1)
  mu <- update.mu(eta)
  
  if (fitype == "adapt" | fitype == "adprf") {
    penweights <- dat$pwt
  } else  penweights <- rep(1, P)
  
  fista.alpha <- numeric(max.iter + 2)
  fista.alpha[1] <- 0
  fista.alpha[2] <- 1
  stepsize <- stepsize.init
  
  loglik.old <- ploglik(y=dat$y, mu=mu, coef, weights=weights, tuning, Lmatrix) 
  pen.old <- penalty(coef=coef, tuning=tuning, Q, penweights=penweights)
  fn.val.old <- -loglik.old + pen.old
  
  d.fn <- 1
  worsen.count <- 0
  
  best.iter <- 0
  best.coef <- coef
  best.pen <- pen.old
  best.l <- loglik.old
  best.l.approx <- loglik.old
  best.fn.val <- fn.val.old
  best.fn.val.approx <- fn.val.old
  best.eta <- eta
  best.mu <- mu
  best.stepsize <- stepsize
#  best.penweights <- penweights
  best.d.fn <- d.fn
  best.proxim <- 0
  
#  print(paste("this is the ", iter.count, "th iteration."))
  
  ## the loop
  while(d.fn > rel.tol){
    iter.count <- iter.count + 1
    if(iter.count > max.iter){
      warning(paste("Maximum number of iterations reached"))
      break
    } 
    coef.old2 <- coef.old1    ## this is used for computations
    coef.old1 <- coef         ## this is used for (potential) convergence checks
    
    ## perform a sanity check to detect if the model is deteriorating towards a
    ## degenerate solution.
#    if(do.sancheck){
#      sancheck(coef = coef, coef.old2 = coef.old2, mu = mu, eta = eta,
#               weights = weights, Proximal.args = Proximal.args) 
#    }
    
    fista.beta <- (fista.alpha[iter.count] - 1)/fista.alpha[iter.count + 1]
    
    ## compute the search point and its gradient: alpha^(t)
    coefs <- Map('+', coef, lapply(Map('-', coef, coef.old2), function(u){fista.beta*u}))  
    coefs <- matrix(unlist(coefs), nrow=Q)
    etas <- update.eta(dat=dat, coef=coefs, weights=weights, iter.count)
    mus <- update.mu(etas)
    grad <- pgradient(dat=dat, mu=mus, coef=coefs, weights=weights, tuning, Lmatrix)
    
    ## stuff for the stepsize search 
    lS <- ploglik(y=dat$y, mu=mus, coef=coefs, weights=weights, tuning, Lmatrix)
    fn.val <- 1
    fn.val.approx <- 0
    scalefac <- 1
    
    ## stepsize search and main computations
    ## a note for the third line: the loss function is the negative loglik, so 
    ## subtracting the gradient of the loss means adding the gradient of loglik.
    while(fn.val > fn.val.approx & stepsize > 1e-30){
      stepsize <- stepsize / scalefac
      coeft <- Map('+', coefs, lapply(grad, function(u){u * stepsize}))
      coeft <- matrix(unlist(coeft), nrow=Q)
      tuning.scaled <- lapply(tuning, function(u){u * stepsize})
      coefprox <- fistaProximal(coef=coeft, tuning=tuning.scaled, penweights=penweights)
      etaprox <- update.eta(dat=dat, coef=coefprox, weights=weights, iter.count)
      muprox <- update.mu(etaprox)
      lprox <- ploglik(y=dat$y, mu=muprox, coef=coefprox, weights=weights, tuning, Lmatrix)
      coefdiff <- Map('-', coefprox, coefs)
      loss.approx <-  lS + Reduce('+', Map(function(u,v){sum(u*v)}, grad, coefdiff))
      proxim <- Reduce('+', lapply(coefdiff, function(u) sum(u^2)))
      if(proxim <= 1e-20) break 
      fn.val <- -lprox
      fn.val.approx <- -loss.approx + (proxim/ stepsize/2)
      scalefac <- 2
#      print(list(dim1=dim(dat$y), dim2=dim(muprox), w=length(weights)))
#      print(paste("this is the ", iter.count, "th iteration."))
    } ## end while(fn.val > ...)
    
    if(stepsize < 1e-30){warning(paste("FISTA couldnt find a feasible stepsize during iteration",
                                              iter.count))
                                break} 

    coef <- coefprox
    eta <- etaprox
    mu <- muprox
    penprox <- penalty(coef=coefprox, tuning=tuning, Q, penweights=penweights)
    
    fn.val <- fn.val + penprox
    fn.val.approx <- fn.val.approx + penprox  
    
    fista.alpha[iter.count+2] <- (1 + sqrt(1 + 4*(fista.alpha[iter.count+1])^2)) / 2
    
    d.fn <- abs(fn.val - fn.val.old)/(rel.tol/10 + abs(fn.val))
    
    if(fn.val < fn.val.old){
      best.iter <- iter.count
      best.coef <- coef
      best.eta <- eta
      best.mu <- mu
      best.fn.val <- fn.val
      best.fn.val.approx <- fn.val.approx
#      best.penweights <- penweights
      best.pen <- penprox
      best.proxim <- proxim
      best.l <- lprox
      best.l.approx <- loss.approx
      best.d.fn <- d.fn
      best.stepsize <- stepsize
      worsen.count <- max(0, worsen.count - 0.5)
    }else{worsen.count <- worsen.count + 1}
    if(worsen.count >= worsen.max){
      warning(paste("fista terminated after", worsen.max,
                     "consecutive iterations in which the objective function worsened"))
      break
    }
    
    fn.val.old <- fn.val
  } ## end while(d.fn > ...)
  
  ## coef was transformed to be a list. before fista returns its output, coef thus
  ## has to be transformed back into the form that coef.init had. if coef.init 
  ## was a vector or matrix, the artificial listversion of coef has length 1.
  if (fitype == "refit" | fitype == "adprf") {
    ind <- nvar(coefprox) 
    if (length(ind) == 1) out <- list(coef=NA) else {
      if(ind[1] == 1) ind0 <- ind[-1]-1 else ind0 <- ind-1
      newdata <- list(y=dat$y, x=dat$x[, ind], Lmatrix=Lmatrix[ind0,ind0], pwt = dat$pwt[ind])
      out <- fista(newdata, tuning = list(0, tuning[[2]]))
      tempcoef <- matrix(0, Q, P+1)
      tempcoef[, ind] <- out$coef
      out$coef <- tempcoef
    }
  }  else {
    out <- list(coef = best.coef,
                stepsize = best.stepsize,
                #              penweights = best.penweights,
                eta = best.eta,
                mu = best.mu,
                loglikval = best.l,
                loglikapprox = best.l.approx,
                penval = best.pen,
                proximval = best.proxim,
                fn.val = best.fn.val,
                fn.val.approx = best.fn.val.approx,
                iter.count = iter.count,
                best.iter = best.iter,
                d.fn = best.d.fn)    
    
  }
  
    return(out)
  }
  
