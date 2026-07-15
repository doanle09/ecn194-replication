rm    (list=ls())

source("R/load libraries and data.R")

recessionstates <- data.table(data)[as.yearmon(date) >= 2008 & as.yearmon(date) <= 2009.5, .(avguenmp = mean(lau_unemp_sa/lau_pop)), keyby=state][order(-avguenmp)][1:10, state]
data$subsetrecession <- with(data, subset48 == T & !(state %in% recessionstates))

# # want rigs per thousand people
# data$rigs_pop90          <- with(data, rigs_land        / (pop1990                  / 10^3) )
# data$rigs_pop00          <- with(data, rigs_land        / (pop2000                  / 10^3) )
# data$rigs_pop10          <- with(data, rigs_land        / (pop2010                  / 10^3) )
# data$rigs_popinterp      <- with(data, rigs_land        / (pop                      / 10^3) )
# data$rigs_poplau         <- with(data, rigs_land        / (lau_pop                  / 10^3) )
# data$rigs_poplaglau      <- with(data, rigs_land        / (lag(lau_pop)             / 10^3) )

#----------------------- Formulas --------------------------------

# NOTE!!!! The functions will break if more than one (1) variable appears in a given term (e.g., rigs_total/pop1990)!!!
# formula has four parts: (1) the formula, (2) the effects, (3) an IV (none) and (4) the multi-way cluster
# Must access felm in lfe package since that is what was modified

# Time Fixed Effects
# Base model, lag-lengths, without clustering, population normalizations, breaks
fmTPopBase  <-   diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:12) + lag(diff(rigs_pop00/1000),0:10)       | state + date | 0 | state + date
fmTPopLags  <- list(                                                                                               
                 diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:12) + lag(diff(rigs_pop00/1000),0:1 )       | state + date | 0 | state + date, 
                 diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:13) + lag(diff(rigs_pop00/1000),0:13)       | state + date | 0 | state + date,
                 diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:24) + lag(diff(rigs_pop00/1000),0:24)       | state + date | 0 | state + date
                 )                                                                                                 
fmTPopNoClus <-  diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:12) + lag(diff(rigs_pop00/1000),0:10)       | state + date 
fmTPopPNorms <- list(                                                                                              
                 diff(emp_popsa90*1000)     ~ 0 + lag(diff(emp_popsa90*1000)    , 1:12) + lag(diff(rigs_pop90/1000    ), 0:10)  | state + date | 0 | state + date,
                 diff(emp_popsa00*1000)     ~ 0 + lag(diff(emp_popsa00*1000)    , 1:12) + lag(diff(rigs_pop00/1000    ), 0:10)  | state + date | 0 | state + date,
                 diff(emp_popsa10*1000)     ~ 0 + lag(diff(emp_popsa10*1000)    , 1:12) + lag(diff(rigs_pop10/1000    ), 0:10)  | state + date | 0 | state + date,
                 diff(emp_popsainterp*1000) ~ 0 + lag(diff(emp_popsainterp*1000), 1:12) + lag(diff(rigs_popinterp/1000), 0:10)  | state + date | 0 | state + date,
                 diff(emp_popsalau*1000)    ~ 0 + lag(diff(emp_popsalau*1000)   , 1:12) + lag(diff(rigs_poplau/1000   ), 0:10)  | state + date | 0 | state + date,
                 diff(emp_popsalaglau*1000) ~ 0 + lag(diff(emp_popsalaglau*1000), 1:12) + lag(diff(rigs_poplaglau/1000), 0:10)  | state + date | 0 | state + date 
                 )                                                                                                 
fmTPopBreak  <-  diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:12) + lag(diff(PRErigs_pop00/1000),0:10) + lag(diff(POSTrigs_pop00/1000),0:10)  | state + date | 0 | state + date

# Controls
# Base model, lag-lengths, without clustering, breaks
fmCPopBase  <-   diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:12) + lag(diff(rigs_pop00/1000),0:10) + lag(diff(natl_emp_popsa00*1000),0:12+1) +  lag(diff(ipi_sa),0:12)  | state | 0 | state + date
fmCPopLags  <- list(
                 diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:12) + lag(diff(rigs_pop00/1000),0:1 ) + lag(diff(natl_emp_popsa00*1000),0:12+1) +  lag(diff(ipi_sa),0:12)  | state | 0 | state + date, 
                 diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:13) + lag(diff(rigs_pop00/1000),0:13) + lag(diff(natl_emp_popsa00*1000),0:13+1) +  lag(diff(ipi_sa),0:13)  | state | 0 | state + date,
                 diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:24) + lag(diff(rigs_pop00/1000),0:24) + lag(diff(natl_emp_popsa00*1000),0:24+1) +  lag(diff(ipi_sa),0:24)  | state | 0 | state + date
                 )
fmCPopBreak  <-  diff(emp_popsa00*1000) ~ 0 + lag(diff(emp_popsa00*1000), 1:12) + lag(diff(PRErigs_pop00/1000),0:10) + lag(diff(POSTrigs_pop00/1000),0:10) + lag(diff(natl_emp_popsa00),0:12+1) +  lag(diff(ipi_sa),0:12)  | state | 0 | state + date                                         
  
#---------------------- Run regressions ---------------------------------
                                           
ResTPopBase    <- lfe:::felm(fmTPopBase,  data=data, subset=data$subset48)
ResTPopBaseSub <- lfe:::felm(fmTPopBase,  data=data, subset=data$subsetrecession)
ResCPopBase    <- lfe:::felm(fmCPopBase,  data=data, subset=data$subset48)

ResTPopBreak   <- lfe:::felm(fmTPopBreak,  data=data, subset=data$subset48)
ResCPopBreak   <- lfe:::felm(fmCPopBreak,  data=data, subset=data$subset48)

ResTPopPNorms <- ResCPopLags <- ResTPopLags <- list()
for (i in 1:length(fmTPopLags)) ResTPopLags[[i]] <- lfe:::felm(fmTPopLags[[i]],  data=data, subset=data$subset48) 
for (i in 1:length(fmCPopLags)) ResCPopLags[[i]] <- lfe:::felm(fmCPopLags[[i]],  data=data, subset=data$subset48) 
for (i in 1:length(fmTPopPNorms)) ResTPopPNorms[[i]] <- lfe:::felm(fmTPopPNorms[[i]],  data=data, subset=data$subset48) 


# # Robustness check excluding states
# states <- levels(factor(unique(subset(data, subset=subset48, select=state)[[1]])))
# RestTLogStates         <- list()
# length(RestTLogStates) <- length(states)
# names(RestTLogStates)  <- states
# 
# cat("\n")
# for (s in states) {  
#   cat(paste(s, "...", sep=""))
#   x                    <- lfe:::felm(fmTPopNoClus,  data=data, subset=(data$subset48) & (index(data)[[1]]!= s) )
#   x$testsOLS           <- testLRMs(x, endogName="lag(diff(emp_popsa00", exogName="rigs_pop00", vcov=x$vcv)
#   RestTLogStates[[s]]  <- felmClean(x, "OLS", scaleBy=1)
# }
# cat("Done w states. \n")
# 


#----------------------- Fix Cluster VCV matrices --------------------------------

ResTPopBase                                          <- fixClusterVCV(ResTPopBase       )
ResTPopBaseSub                                       <- fixClusterVCV(ResTPopBaseSub    )
ResCPopBase                                          <- fixClusterVCV(ResCPopBase       )
ResTPopBreak                                         <- fixClusterVCV(ResTPopBreak      )
ResCPopBreak                                         <- fixClusterVCV(ResCPopBreak      )
for (i in 1:length(fmTPopLags))   ResTPopLags[[i]]   <- fixClusterVCV(ResTPopLags[[i]]  )
for (i in 1:length(fmCPopLags))   ResCPopLags[[i]]   <- fixClusterVCV(ResCPopLags[[i]]  )
for (i in 1:length(fmTPopPNorms)) ResTPopPNorms[[i]] <- fixClusterVCV(ResTPopPNorms[[i]])

#----------------------- Add multipliers and IRFS --------------------------------

ResTPopBase                                          <- addStats(ResTPopBase       , endogName="lag(diff(emp_popsa", exogName="rigs_pop" )
ResTPopBaseSub                                       <- addStats(ResTPopBaseSub    , endogName="lag(diff(emp_popsa", exogName="rigs_pop" )
ResCPopBase                                          <- addStats(ResCPopBase       , endogName="lag(diff(emp_popsa", exogName="rigs_pop" )
ResTPopBreak                                         <- addStats(ResTPopBreak      , endogName="lag(diff(emp_popsa", exogName="rigs_pop" )
ResCPopBreak                                         <- addStats(ResCPopBreak      , endogName="lag(diff(emp_popsa", exogName="rigs_pop" )
for (i in 1:length(fmTPopLags))   ResTPopLags[[i]]   <- addStats(ResTPopLags[[i]]  , endogName="lag(diff(emp_popsa", exogName="rigs_pop" )
for (i in 1:length(fmCPopLags))   ResCPopLags[[i]]   <- addStats(ResCPopLags[[i]]  , endogName="lag(diff(emp_popsa", exogName="rigs_pop" )
for (i in 1:length(fmTPopPNorms)) ResTPopPNorms[[i]] <- addStats(ResTPopPNorms[[i]], endogName="lag(diff(emp_popsa", exogName="rigs_pop" )

ResTPopBreak$irfElts       <- makeIRFXElements(ResTPopBreak)
ResTPopBreak$irf           <- makeIRFX(ResTPopBreak$irfElts)
ResTPopBreak$varirfOLS     <- varIRFX( ResTPopBreak$irfElts, ResTPopBreak$vcv)
ResTPopBreak$varirfClu     <- varIRFX( ResTPopBreak$irfElts, ResTPopBreak$clustervcv)
ResTPopBreak$varirfSCC     <- varIRFX( ResTPopBreak$irfElts, DriscollKraay(ResTPopBreak))

ResCPopBreak$irfElts       <- makeIRFXElements(ResCPopBreak)
ResCPopBreak$irf           <- makeIRFX(ResCPopBreak$irfElts)
ResCPopBreak$varirfOLS     <- varIRFX( ResCPopBreak$irfElts, ResCPopBreak$vcv)
ResCPopBreak$varirfClu     <- varIRFX( ResCPopBreak$irfElts, ResCPopBreak$clustervcv)
ResCPopBreak$varirfSCC     <- varIRFX( ResCPopBreak$irfElts, DriscollKraay(ResCPopBreak))


#------------------------ Test for structural break -------------------------------

# All betas
RT <- cbind(matrix(0, nrow=11, ncol=12), diag(11), -diag(11) )
# RC <- cbind(RT, matrix(0, nrow=11, ncol=13*2) )
test <- "Chisq"

linearHypothesis(ResTPopBreak, RT, vcov.=ResTPopBreak$vcv           , test=test)
linearHypothesis(ResTPopBreak, RT, vcov.=ResTPopBreak$clustervcv    , test=test)
linearHypothesis(ResTPopBreak, RT, vcov.=DriscollKraay(ResTPopBreak), test=test)

# linearHypothesis(ResCPopBreak, RC, vcov.=ResCPopBreak$vcv           , test=test)
# linearHypothesis(ResCPopBreak, RC, vcov.=ResCPopBreak$clustervcv    , test=test)
# linearHypothesis(ResCPopBreak, RC, vcov.=DriscollKraay(ResCPopBreak), test=test)


# Sum of betas
RT <- c(rep(0,12), rep(-1,11), rep(1,11))
# RC <- c(RT, rep(0,26))

RT %*% coef(ResTPopBreak)
# RC %*% coef(ResCPopBreak)

linearHypothesis(ResTPopBreak, RT, vcov.=ResTPopBreak$vcv           , test=test)
linearHypothesis(ResTPopBreak, RT, vcov.=ResTPopBreak$clustervcv    , test=test)
linearHypothesis(ResTPopBreak, RT, vcov.=DriscollKraay(ResTPopBreak), test=test)

# linearHypothesis(ResCPopBreak, RC, vcov.=ResCPopBreak$vcv           , test=test)
# linearHypothesis(ResCPopBreak, RC, vcov.=ResCPopBreak$clustervcv    , test=test)
# linearHypothesis(ResCPopBreak, RC, vcov.=DriscollKraay(ResCPopBreak), test=test)


c(rep(0,12), rep(1,11), rep(0,11)           ) %*% coef(ResTPopBreak)
c(rep(0,12), rep(0,11), rep(1,11)           ) %*% coef(ResTPopBreak)
# c(rep(0,12), rep(1,11), rep(0,11), rep(0,26)) %*% coef(ResCPopBreak)
# c(rep(0,12), rep(0,11), rep(1,11), rep(0,26)) %*% coef(ResCPopBreak)

#----------------------- Create LRMs for Struct Break cases --------------------------------

#Time FE

endogName <- "lag(diff(emp_popsa"
exogName  <- "PRErigs_pop" 
TpreTestLRM  <- list(
  OLS     = testLRMs(ResTPopBreak, endogName=endogName, exogName=exogName, vcov=ResTPopBreak$vcv),
  Cluster = testLRMs(ResTPopBreak, endogName=endogName, exogName=exogName, vcov=ResTPopBreak$clustervcv),
  SCC     = testLRMs(ResTPopBreak, endogName=endogName, exogName=exogName, vcov=ResTPopBreak$vcvSCC)
)

exogName  <- "POSTrigs_pop" 
TpostTestLRM  <- list(
  OLS     = testLRMs(ResTPopBreak, endogName=endogName, exogName=exogName, vcov=ResTPopBreak$vcv),
  Cluster = testLRMs(ResTPopBreak, endogName=endogName, exogName=exogName, vcov=ResTPopBreak$clustervcv),
  SCC     = testLRMs(ResTPopBreak, endogName=endogName, exogName=exogName, vcov=ResTPopBreak$vcvSCC)
)

# # Controls
# 
# exogName  <- "PRErigs_pop" 
# CpreTestLRM  <- list(
#   OLS     = testLRMs(ResCPopBreak, endogName=endogName, exogName=exogName, vcov=ResCPopBreak$vcv),
#   Cluster = testLRMs(ResCPopBreak, endogName=endogName, exogName=exogName, vcov=ResCPopBreak$clustervcv),
#   SCC     = testLRMs(ResCPopBreak, endogName=endogName, exogName=exogName, vcov=ResCPopBreak$vcvSCC)
# )
# 
# exogName  <- "POSTrigs_pop" 
# CpostTestLRM  <- list(
#   OLS     = testLRMs(ResCPopBreak, endogName=endogName, exogName=exogName, vcov=ResCPopBreak$vcv),
#   Cluster = testLRMs(ResCPopBreak, endogName=endogName, exogName=exogName, vcov=ResCPopBreak$clustervcv),
#   SCC     = testLRMs(ResCPopBreak, endogName=endogName, exogName=exogName, vcov=ResCPopBreak$vcvSCC)
# )
# 

ResTPopBase$testsOLS
ResTPopBase$testsclu
ResTPopBase$testsSCC
ResTPopBase$testsclu[2,2]*1.96

TpreTestLRM
TpostTestLRM

