

# These functions are benchmarked against the Arma{signal} and impz{signal}
#   functions in package:signal. To validate, form coefficient vectors and
#   create a model:   model <- Arma(b = exogcoefs, a = c(1,- endogcoefs ))
#   The IRF is calculated by impz(model, n=numirfs)



# Extracts vector of coefficients from coef(felm) that correspond to varname as
#   listed in original data (varname)
extractTermCoefficients <- function(varname, coef) {

  positions <- grep(varname, names(coef), fixed=T)
  rval <- coef[positions]
  attr(rval, "varname") <- varname
  
  lags <- as.integer(gsub(varname, "", names(rval), fixed=T))
  if (length(lags)==0 || is.na(lags)==T) lags <- rep(0,1)
  
  attr(rval, "lags") <- lags
  return(rval)
}

# Given a set of coefs with lags attributes, makes the vector of coefs
makeCoefVector <- function(x, maxlag = NULL, minlag=NULL ){

  lags <- attr(x,"lags")
  if (is.null(maxlag)) maxlag <- max(lags)
  if (is.null(minlag)) minlag <- min(lags)    
  len <- maxlag - minlag + 1
  
  rval <- rep(0,len)
  rval[lags - minlag + 1 ] <- x[1:length(x)]
  names(rval)[lags - minlag + 1 ] <- names(x)
  attr(rval, "lags") <- sort(lags)
  
  return(rval)
}


# Form A, B, J, beta matrices as in Lutkepohl (2007) p403
makeIRFXElements <- function(x){

  # rhs and lhs sides of formulas
  Fm       <- pFormula(eval(x$call$formula))
  rhsFm    <- attr(Fm,"rhs")[[1]] 
  lhsFm    <- attr(Fm,"lhs")[[1]]
  Fm2      <- Formula:::paste_formula(lhsFm,rhsFm)
  
  # Terms in the model
  rhs.labels <- attr(terms(Fm2),"term.labels")
  lhs.label  <- capture.output(lhsFm)
  LhsOnRhs   <- grep(lhs.label, rhs.labels, fixed=T)
  if (is.na(LhsOnRhs)) stop("cannot find lhs variable on rhs")
  idx <- 1:length(rhs.labels)
  idx <- idx[idx!=LhsOnRhs]
  exog.labels <- rhs.labels[idx]
  endog.labels <- rhs.labels[LhsOnRhs]
 
  # List where each elt is coefficients for one of variables
  coeflist <- lapply(rhs.labels, extractTermCoefficients, coef=coef(x) ) 
  names(coeflist) <- rhs.labels
    
  # Vector of endogenous coefs w/ lags
  Avec     <- makeCoefVector(coeflist[[LhsOnRhs]])
  coeflist[[LhsOnRhs]] <- NULL
  
  # Dimensions for matrix (see Luktepohl p 403)
  K <- 1
  p <- length(Avec)
  Kp <- K*p
  s  <- max(sapply(coeflist, function(i) max(attr(i, "lags"))))
  M  <- length(exog.labels)
  Ms <- M*s

  # List of coef vectors for exogenous variables
  exogvec     <- lapply(coeflist, makeCoefVector, minlag=0, maxlag=s)  # possibly fix here
  Bvec        <- as.vector(do.call(rbind, exogvec))  
  names(Bvec) <- as.vector(do.call(rbind,lapply(exogvec, attr, "names"))) 
  
  # browser()
  
  # Make companion A matrix
  A <- matrix(0, nrow=(Kp+Ms), ncol=(Kp+Ms))
  A[1, 1:Kp] <- Avec
  A[1, 1:Ms + Kp] <- Bvec[(M+1):(Ms+M)]
  A[2:Kp, (1:Kp-1)] <- diag(Kp-1)
  if(Ms > M) {
    A[Kp + M + 1:(Ms-M), Kp+(1:(Ms-M))] <- diag(Ms-M)
  }
  colnames(A) <- c(names(Avec), names(Bvec)[(M+1):(Ms+M)])
  
  # Make companion B matrix
  B <- matrix(0, nrow=(Kp+Ms), ncol=M )
  B[1, 1:M] <- Bvec[1:M]
  B[Kp+1:M, 1:M] <- diag(M)
  colnames(B) <- names(Bvec)[1:M]

  # Make J selection matrix
  J <- matrix(0, nrow=K, ncol=(Kp+Ms))
  J[1:K, 1:K] <- diag(K)

  # Coefficient vector
  beta <- c(A[1,], B[1,])
  
  rval <- list(beta=beta, A=A, B=B, J=J, 
                lhs.label=lhs.label, LhsOnRhs=LhsOnRhs,
                exog.labels=exog.labels, endog.labels=endog.labels,
                K=K, p=p, Kp=Kp, s=s, M=M, Ms=Ms)
  return(rval)
}

# This takes output from IRFElements
makeIRFX <- function(x, maxIR=24) {
    irf <- as.matrix( sapply(0:maxIR, function(i) x$J %*% (x$A%^%i) %*% x$B ) )
    if(ncol(irf)>1) irf <- t(irf)
    rownames(irf) <- paste("t=", 0:maxIR, sep="")
    colnames(irf) <- x$exog.labels
    return(irf)
}

makeCIRFX <- function(irf) apply(irf, 2, cumsum)

# Makes vcov for use in calculating var(irf)
betavcov <- function(x,vcov){
  oldnames <- colnames(vcov)
  vcov <- cbind(rbind(vcov,0),0)
  colnames(vcov) <- rownames(vcov) <- c(oldnames, "zero")
  newnames <- names(x$beta)
  newnames[is.na(newnames)] <- "zero"
  
  rval <- vcov[newnames,newnames]
  if(!is.na(newnames["zero"])) {
    newnames["zero"] <- NA
    colnames(rval) <- rownames(rval) <- newnames
  }
  
  return(rval)
}


varIRFX <- function(x, vcov, maxIR = 24) {
  
  require(arrayhelpers)
  
  K  <- x$K
  M  <- x$M
  Ms <- x$Ms
  p  <- x$p
  s  <- x$s
  Kp <- x$Kp
  A  <- x$A
  B  <- x$B
  J  <- x$J
  
  dimGi <- c(M*K, (Kp+Ms)*K + M*K )
  
  vcov <- betavcov(x, vcov) # Fixes vcov so that it is conformable w/ beta
  
  # Make set of Gi = delta vec(Di) / delta beta'
  G      <- array(, dim=c(dimGi, maxIR+1))
  G[,,1] <- cbind( matrix(0, nrow=dimGi[1], ncol=dimGi[2]-M*K), diag(M*K) )
  
  for (i in 1:maxIR) {
    Gi <- array(0, c(dimGi, i) ) 
    for (j in 1:i-1) {
      left <- kronecker( t(B) %*% (t(A)%^%(i-1-j)) ,  J %*% (A%^%j) %*% t(J) )
      right <- kronecker( diag(M), J %*% (A%^%i) %*% t(J) )
      Gi[,,j+1] <- cbind( left, right)
    }
    
    G[,,i+1] <- rowSums(Gi, dim=2, drop=F)
  }
  
  
  # Variance of IRFX
  #browser()
  varIRFX <- array(0, c(M*K, M*K, maxIR+1))     
  for (i in 0:maxIR) {
    Gi <- matrix(G[,,i+1], nrow=dim(G)[1], ncol=dim(G)[2])
    varIRFX[,,i+1] <- Gi %*% vcov %*% t(Gi)
  }
  
  # Variance of CIRFX
  varCIRFX <- array(0, c(M*K, M*K, maxIR+1))      
  for (i in 0:maxIR) {
    Gis <- slice(G, k=0:i+1, drop=FALSE)
    cumG            <- matrix( rowSums(Gis, dim=2, drop=F), nrow=dim(G)[1], ncol=dim(G)[2])
    varCIRFX[,,i+1] <- cumG %*% vcov %*% t(cumG)
  }
  
  # Std errors
  seIRFX   <- matrix(0, nrow=dim(varIRFX)[3], ncol=dim(varIRFX)[1])
  seCIRFX  <- matrix(0, nrow=dim(varIRFX)[3], ncol=dim(varCIRFX)[1])
  for(i in 1:dim(varIRFX)[3]) {
    seIRFX[i,]  <- sqrt(diag(as.matrix( varIRFX[,,i])))
    seCIRFX[i,] <- sqrt(diag(as.matrix(varCIRFX[,,i])))
  }
  
  # Name things
  dimnames(varCIRFX) <- dimnames(varIRFX) <- list( x$exog.labels, x$exog.labels, paste("t=",0:maxIR,sep=""))
  dimnames(seCIRFX) <- dimnames(seIRFX)  <- list( paste("t=",0:maxIR,sep=""), x$exog.labels )
  # dimnames(varIRFX) <- list( colnames(B), colnames(B), paste("t=",0:maxIR,sep=""))
  # dimnames(seIRFX)  <- list( paste("t=",0:maxIR,sep=""), colnames(B) )
  
  rval <- list(seIRFX=seIRFX, varIRFX=varIRFX, varCIRFX=varCIRFX, seCIRFX=seCIRFX)
  return(rval)
  
}  



plotIRFX <- function(irf, varirf, mainpre="", pch=20, ... ) {
  lower <- irf - 1.96*abs(varirf$seIRFX)
  upper <- irf + 1.96*abs(varirf$seIRFX)
   
  timeaxis <- 1:nrow(irf)-1
  for(j in 1:ncol(lower)) {
    y <- cbind(irf[,j], lower[,j], upper[,j])
    matplot(timeaxis, cbind(irf[,j], lower[,j], upper[,j]), 
      type="b", lty=c(1,2,2),  col=c(1,2,2),
      pch=pch,  mar=c(5, 4, 4, 2) + 0.1, 
      xlab="", main=paste(mainpre, colnames(irf)[j], sep=""), ylab="", ... )
    abline(h=0)
  }
}
 