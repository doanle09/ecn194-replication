


stateACFs <- function(model) {
  
  res    <- residuals(model)
  idx    <- index(res, 1)
  lev    <- levels(idx)
  
  acfs <- lapply(lev, function(i){
    ac <- acf(res[idx==i], plot=FALSE)
    x <- as.matrix(ac$acf)  
    rownames(x) <- 0:24
    colnames(x) <- i
    attributes(x)$n.used <- ac$n.used
    return(x)
  }   )
  names(acfs) <- lev
  m <- melt(acfs)
  m$L1 <- NULL
  names(m) <- c("lag","state","value")
  m$min <- m$value * (m$value < 0)
  m$max <- m$value * (m$value > 0)
  m$value <- NULL
  ci <- qnorm((1 + .95)/2)/sqrt(attributes(acfs[[1]])$n.used)
  p <- ggplot(m, aes(x=lag, ymin=min, ymax=max)) + 
    geom_linerange() + 
    facet_wrap(~ state, nrow=10) + 
    geom_hline(yintercept=c(-ci,ci), colour="red", linetype=2) + 
    geom_hline(yintercept=0) +
    scale_x_continuous(breaks=0:10*12) +
    xlab("") + ylab("")
  return(p)
}






plotIRF <- function(irf, varirfs, cum=FALSE, varlabels=NULL) {
  
  if (cum) {
    
    cirf            <- apply(irf, 2, cumsum)
    dimnames(cirf)  <- dimnames(irf)
    irf             <- cirf
    whichse         <-"seCIRFX"
  }
  else{
    whichse <- "seIRFX"
  }
  namesCI <- names(varirfs) 
  times <- 1:length(irf)-1
  
  irfDF <- list(IRF   = irf,
                Upper = lapply(varirfs, function(v)    cbind(irf + 1.96 * v[[whichse]] ) ),
                Lower = lapply(varirfs, function(v)    cbind(irf - 1.96 * v[[whichse]] ) )
  )
  m           <- melt(irfDF) 
  
  m$Lead      <- as.numeric(gsub("t=","",m$Var1))
  m$exogVar   <- m$Var2
  m$CI        <- m$L1
  m$SEType    <- m$L2
  m$Var1  <- m$Var2 <- m$L1 <- m$L2 <- NULL
  
  if (!is.null(varlabels)) levels(m$exogVar) <- varlabels
  
  mysize <- 1
  p <- ggplot(m, aes(x=Lead, y=value, col=SEType, shape=SEType )) +
    geom_line( data=subset(m, CI=="Upper"), size=mysize) + 
    geom_line( data=subset(m, CI=="Lower"), size=mysize) + 
    geom_point(data=subset(m, CI %in% c("Upper","Lower")) , size=3.5   ) +
    xlab("") + ylab("") + 
    geom_hline(yintercept=0) + 
    scale_x_continuous(breaks=0:4*6) + 
    labs(col="95% Conf.\nInterval",shape="95% Conf.\nInterval") +
    # labs(col="",shape="") +
    # scale_colour_brewer(palette="Set1") +
    geom_line( data=subset(m, CI=="IRF")  , size=mysize) +
    guides(col = guide_legend(reverse = TRUE), shape=guide_legend(reverse = TRUE) )
  
  if(length(levels(m$exogVar))>1) p <- p + facet_wrap("exogVar")
  attributes(p)$data <- m
  return(p)
  
} 
  
#------------------------------------------------------------------------------
 
  
mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="fct") { 
    value[value=="aIRF"] <- "IRF"
    value[value=="CIRF"] <- "CIRF"
  } 
  if (var=="pickFormula") {
    value[value=="Lag1"] <- "Truncated lags"    
    value[value=="Lag13"] <- "13 lags"    
    value[value=="Lag24"] <- "24 lags"    
  }
  return(value)
}