
# http://adv-r.had.co.nz/Namespaces.html
# http://tolstoy.newcastle.edu.au/R/e6/help/09/03/6725.html
# http://rfunction.com/file/BuildingRPackagesPaper.pdf
# http://cran.r-project.org/doc/Rnews/Rnews_2003-3.pdf#page=29



felm <- function(formula, data, iv=NULL, clustervar=NULL, exactDOF=FALSE, subset, na.action, contrasts=NULL,...) {

  # In a later version we're moving clustervar and iv out of the argument list, into the ... list
  # For now, make sure people specify them by name, so that it will continue to work
  # I.e. if they're non-null, but not in pmatch(names(sys.call())), they're
  # without names.  Warn them.
  knownargs <- c('iv','clustervar','cmethod')
  sc <- names(sys.call())[-1]
  named <- knownargs[pmatch(sc,knownargs)]
  for(arg in c('iv', 'clustervar')) {
    if(!is.null(eval(as.name(arg))) && !(arg %in% named)) {
        warning("Please specify the '",arg,"' argument by name, or use a multi part formula. Its position in the argument list will change in a later version")
      }
  }

  mf <- match.call(expand.dots = FALSE)

  # Currently there shouldn't be any ... arguments
  # check that the list is empty

#  if(length(mf[['...']]) > 0) stop('unknown argument ',mf['...'])
  
  # When moved to the ... list, we use this:
  # we do it right away, iv and clustervar can't possibly end up in ... yet, not with normal users

  cmethod <- 'cgm'

  args <- list(...)
  ka <- knownargs[pmatch(names(args),knownargs, duplicates.ok=FALSE)]
  names(args)[!is.na(ka)] <- ka[!is.na(ka)]
  env <- environment()
  lapply(intersect(knownargs,ka), function(arg) assign(arg,args[[arg]], pos=env))

  if(!(cmethod %in% c('cgm','gaure'))) stop('Unknown cmethod: ',cmethod)

  # also implement a check for unknown arguments
  unk <- setdiff(names(args), knownargs)
  if(length(unk) > 0) stop('unknown arguments ',paste(unk, collapse=' '))


  if(missing(data)) data <- environment(formula)
  pf <- parent.frame()
  pform <- parseformula(formula,data)
  

  if(!is.null(iv) && !is.null(pform[['iv']])) stop("Specify EITHER iv argument(deprecated) OR multipart terms, not both")
  if(!is.null(pform[['cluster']]) && !is.null(clustervar)) stop("Specify EITHER clustervar(deprecated) OR multipart terms, not both")
  if(!is.null(pform[['cluster']])) clustervar <- structure(pform[['cluster']], method=cmethod)

  if(is.null(iv) && is.null(pform[['iv']])) {
    # no iv, just do the thing
    fl <- pform[['fl']]
    formula <- pFormula(pform[['formula']])
    mf[['formula']] <- formula
#    if(!is.null(clustervar)) warning("argument clustervar is deprecated, use multipart formula instead")
    psys <- project(mf,fl,data,contrasts,clustervar,pf)
    gc()

    z <- doprojols(psys,exactDOF=exactDOF)
    z$xz <- psys$yxz$x
    z$yz <- psys$yxz$y
    rm(psys)
    gc()
    z$call <- match.call()
    return(z)
  }

}

reassignInPackage ("felm", pkgName="lfe", felm, keepOld=FALSE)
rm(felm)


# --------------------------------------------------------------------------


bread.felm <- function(x, ...) {
  return(x$inv)
}

estfun.felm <- function(x, ...) {
  idx <- index(x$y)
  h <- x$xz
  hsum <- aggregate(h, by=list(idx[,2]), sum)
  rval <- as.matrix(zoo(hsum[,2:ncol(hsum)], order.by=hsum[,1]))
  return(rval)  
}

sandwich.default <- sandwich:::sandwich

sandwich <- function(x, ...)
{
  UseMethod("sandwich")
}

assignInNamespace("sandwich", sandwich, pos="package:sandwich")
rm(sandwich)

# --------------------------------------------------------------------------

sandwich.felm <- function(x, bread. = bread, meat. = meat, ...)
{
  if(is.function(bread.)) bread. <- bread.(x)
  if(is.function(meat.)) meat. <- meat.(x, ...)
  n <- NROW(estfun(x))
  return(bread. %*% (n* meat.) %*% bread.)
}

# --------------------------------------------------------------------------

DriscollKraay <- function(x, lag = NULL,
                      order.by = NULL, prewhite = FALSE, adjust = FALSE, 
                      diagnostics = FALSE, sandwich = TRUE, ar.method = "ols", data = list(),
                      verbose = TRUE)
{
  idx <- attr(x$yz, "index")    
  pdim <- pdim(idx[[1]], idx[[2]])
  n <- pdim$nT$n
  T. <- pdim$nT$T
  nT <- pdim$nT$N  
  
  if(is.null(lag)) lag<-floor(T.^(1/4))
    
  if(verbose) cat(paste("\nLag truncation parameter chosen:", lag, "\n"))
  
  myweights <- seq(1, 0, by = -(1/(lag + 1)))
  vcv <- vcovHAC(x, order.by = order.by, prewhite = prewhite,
          weights = myweights, adjust = adjust, diagnostics = diagnostics,
          sandwich = sandwich, ar.method = ar.method, data = data)
  colnames(vcv) <- rownames(vcv) <- names(coef(x))
  return(vcv)
}

# --------------------------------------------------------------------------

model.matrix.felm <- function(x) {
  return(x$xz)
}


logLik.felm <- function(x) {
  N <- x$N
  rval <- -(N/2)*( log(2*pi) + log(crossprod(x$residuals)) - log(N) + 1 )
  attr(rval, "nobs") <- N
  attr(rval, "df") <- x$p 
  class(rval) <- "logLik"
  return(rval)
}


fixClusterVCV <- function(x) {
  
  if ( any(eigen(x$clustervcv)$values< 0) ) {
    cat("\nSome eigenvalues of clustervcv are negative. Fixing...")
    x$clustervcvOLD <- x$clustervcv
    e <- eigen(x$clustervcv)
    x$clustervcv <- e$vectors %*% diag(e$values * (e$values > 0)) %*% t(e$vectors)
    dimnames(x$clustervcv) <- dimnames(x$clustervcvOLD)
    x$cse <- sqrt(diag(x$clustervcv))
    x$ctval <- x$beta / x$cse
    x$cpval <- 2 * pt(abs(x$tval), df=x$df, lower.tail = FALSE)
  }
  else {
    cat("\nEigen values all > 0. OK.")
  }
  return(x)
}


