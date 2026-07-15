


testLRMs <- function(obj, endogName, exogName, vcov=NULL) {
  
  coef        <- coef(obj)
  if(is.null(vcov)) vcov <- obj$vcv
  
  whichEndog  <- grep(gsub(" ","",endogName), gsub(" ","", names(coef)), fixed=T)  
  whichExog   <- grep(gsub(" ","",exogName),  gsub(" ","", names(coef)), fixed=T)
  
  if(length(whichEndog) == 0 || length(whichExog) == 0) stop("variables not found!")
  
  gsub(" ","", names(coef))
  
  sumEndog    <- paste(paste("b", whichEndog, sep=""), collapse=" - ")
  sumExog     <- paste(paste("b", whichExog , sep=""), collapse=" + ")
    
  useCoef <- min(whichEndog, whichExog):max(whichEndog, whichExog)
  
  fm          <- vector(mode="character", length=2)  
  fm[1]       <- sumExog
  fm[2]       <- paste("(", sumExog, ")/(1 -", sumEndog,")")
  
  p.names     <- paste("b", 1:length(coef), sep="")
  names(coef) <- p.names
  dimnames(vcov)[[1]] <- dimnames(vcov)[[2]] <- p.names
  
  dm           <- matrix(0,nrow=2,ncol=4)
  rownames(dm) <- c("SumBeta","LRM")
  colnames(dm) <- c("Estimate", "SE", "T-stat", "Pr(est > 0)")    
    
  jnk          <- as.matrix(t(sapply(fm, function(i) 
    deltaMethod(coef[useCoef], i, vcov.=vcov[useCoef,useCoef], parameterNames=p.names[useCoef] )
  )))
  mode(jnk)    <- "numeric"
  dm[,1:2]     <- jnk
  dm[,3]       <- dm[,1]/dm[,2]
  dm[,4]       <- pt(dm[,3], df=obj$df, lower.tail=FALSE)
  
  return(dm)
}



addStats <- function(model, endogName, exogName, DK=TRUE) {
  
  model$vcvSCC     <- DriscollKraay(model)
  model$testsOLS   <- testLRMs(model, endogName=endogName, exogName=exogName, vcov=model$vcv)
  model$testsclu   <- testLRMs(model, endogName=endogName, exogName=exogName, vcov=model$clustervcv)
  model$testsSCC   <- testLRMs(model, endogName=endogName, exogName=exogName, vcov=model$vcvSCC)
    
  model$irfElts    <- makeIRFXElements(model)
  model$irf        <- makeIRFX(model$irfElts)
  model$cirf       <- makeCIRFX(model$irf)
  
  model$varirfOLS  <- varIRFX(model$irfElts,model$vcv)
  if (!is.null(model$clustervcv)){
    model$varirfClu  <- varIRFX(model$irfElts,model$clustervcv)
  }
  if (DK) {
    model$varirfSCC  <- varIRFX(model$irfElts,model$vcvSCC)
  }
  
  return(model)
}


addIRF <- function(model, endogName, exogName, DK=TRUE) {
  
  model$vcvSCC     <- DriscollKraay(model)
  model$irfElts    <- makeIRFXElements(model)
  model$irf        <- makeIRFX(model$irfElts)
  model$cirf       <- makeCIRFX(model$irf)
  
  model$varirfOLS  <- varIRFX(model$irfElts,model$vcv)
  if (!is.null(model$clustervcv)){
    model$varirfClu  <- varIRFX(model$irfElts,model$clustervcv)
  }
  if (DK) {
    model$varirfSCC  <- varIRFX(model$irfElts,model$vcvSCC)
  }
  
  return(model)
}



doLinTests <- function(R,model,test="Chisq") {

  vcvs <- list(OLS=model$vcv, Cluster=model$clustervcv, SCC=DriscollKraay(model))
  lHyps <- lapply(vcvs, function(i) linearHypothesis(model, R, vcov.=i, test=test))
  
  out <- list()
  
  out$Rb <- R %*% coef(model)
  out$Df <- lHyps[[1]]$Df[2]
  out$Res.Df <- lHyps[[1]]$Res.Df
  out$Chisq        <- sapply(lHyps, function(i) i$Chisq[2])
  out$`Pr(>Chisq)` <- sapply(lHyps, function(i) i$`Pr(>Chisq)`[2])
  
  return(out)
}




#------------------------------------------------------------------------------




# extension for plm objects (from the plm package)
extract.felm <- function(model, include.bic = TRUE, include.pdim=TRUE, 
                         use.tstats=FALSE, indicate=NULL, indicate.fe=NULL, include.LRM=FALSE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(s$coef)
  coefficients      <- s$coef[, 1]
  standard.errors   <- s$coef[, 2]
  significance      <- s$coef[, 4]
  if(use.tstats)    standard.errors <- coefficients / standard.errors
  
  rs <- s$r.squared[1]
  adj <- s$r.squared[2]
  n <- length(s$resid)
  
  pd <- plm:::pdim.pdata.frame(model$yz)

  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  if(!is.null(indicate.fe)) {
    feInModel   <- match(indicate.fe, names(model$fe) )
    feInModel   <- !is.na(feInModel)
    gof         <- c(gof, feInModel)
    gof.names   <- c(gof.names, names(indicate.fe))
    gof.decimal <- c(gof.decimal, rep(FALSE, length(feInModel)) )
  }
  if(!is.null(indicate)) {
    termsInModel <- sapply(indicate, function(i) any(grepl(i, names(coef(model)), fixed=T)) )
    gof          <- c(gof, termsInModel)
    gof.names    <- c(gof.names, names(indicate))
    gof.decimal  <- c(gof.decimal, rep(FALSE, length(indicate)))
  }
  if (include.pdim == TRUE) {
    gof <- c(gof, pd$nT$n, pd$nT$T, pd$nT$N, pd$balanced)
    gof.names <- c(gof.names, "N", "T", "Obs.", "Balanced Panel")
    gof.decimal <- c(gof.decimal, rep(FALSE, 4))
  }
  if (include.LRM == TRUE) {
    tests       <- model$MultiplierTests[,c(1:2,4)]
    gof         <- c(gof, as.vector(t(tests)))
    gof.names   <- c(gof.names, "SumBeta","SE(SumBeta)","Pr(SumBeta < 0)","LRM","SE(LRM)","PR(LRM < 0)")
    gof.decimal <- c(gof.decimal, rep(TRUE,6))
  }  
  if (include.bic == TRUE) {
    gof <- c(gof, BIC(model))
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  tr <- createTexreg(
    coef.names = coefficient.names, 
    coef = coefficients, 
    se = standard.errors, 
    pvalues = significance, 
    gof.names = gof.names, 
    gof = gof, 
    gof.decimal = gof.decimal
  )
  
  return(tr)
}

setMethod("extract", signature = className("felm", "lfe"), definition = extract.felm)



felmClean <- function(x, se.type="OLS", scaleBy=100) {
  
  kill <- c("robustvcv","clustervcv","cse","ctval","cpval","clustervar","rse","rtval","rpval","vcvSCC","testOLS","testClu","testSCC","varirfOLS","varirfClu","varirfSCC")
  keep <- c("vcv","se","tval","pval")
  scale <- c("beta","se")
  
  if (se.type=="OLS") {
    x$se   <- sqrt(diag(x$vcv))
    x$tval <- x$beta / x$se
    x$pval <- 2 * pt(abs(x$tval), df=x$df, lower.tail = FALSE)
    x$MultiplierTests <- x$testsOLS
    
  } else if (se.type=="cluster") {
    x$vcv  <- x$clustervcv
    x$se   <- x$cse
    x$tval <- x$ctval
    x$pval <- x$cpval
    x$MultiplierTests <- x$testsclu
    
  } else if (se.type=="SCC") {
    x$vcv  <- x$vcvSCC
    x$se   <- sqrt(diag(x$vcv))
    x$tval <- x$beta / x$se
    x$pval <- 2 * pt(abs(x$tval), df=x$df, lower.tail = FALSE)
    x$MultiplierTests <- x$testsSCC
  }
  
  x[kill] <- NULL  
  x$beta  <- x$beta*scaleBy
  x$coef  <- x$coef*scaleBy
  x$se    <- x$se*scaleBy
  x$MultiplierTests[,1:2] <- x$MultiplierTests[,1:2]*scaleBy
  return(x)
}

          
#------------------------------------------------------------------------------


