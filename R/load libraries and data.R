rm(list=ls())

library(data.table)
library(foreign)
library(plm)
library(lfe)
library(R.utils)
library(sandwich)
library(zoo)
library(lmtest)
library(expm)
library(car)
library(plyr)
library(reshape2)
library(ggplot2)
library(texreg)
library(arrayhelpers)
library(scales)

# see http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/
# curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
# sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
# library(devtools)
# devtools::install_url('https://cran.r-project.org/src/contrib/Archive/lfe/lfe_1.8-1441.tar.gz')

# These load in modified felm routines and IRF estimation
source("R/functions - estimation.R")
source("R/functions - post estimation.R")
source("R/functions - irf.R")
source("R/functions - plotting.R")


# read in stata dataset 

# load("data/addedDIData.Rdata")

data <- as.data.table( read.dta("dat/state.dta") )
natl <- as.data.table( read.dta("dat/natl.dta") )

# convert dates from "months since Jan 1960" to R yearmon{zoo} format
data$date           <- as.yearmon(1960 + data$date/12)
natl$date           <- as.yearmon(1960 + natl$date/12)
natl <- natl[date >= 1987,]

statepop            <- data[date %in% (c(1990,2000,2010)+.5), .(date,state,pop)]
natpop              <- natl[date %in% (c(1990,2000,2010)+.5), .(date,natl_population)]
statepop            <- dcast(statepop, state ~ date, value.var="pop")
natpop              <- dcast(natpop  , .     ~ date, value.var="natl_population")
names(statepop)     <- names(natpop) <- c("state","pop1990","pop2000","pop2010")
names(natpop)       <- paste("natl_", names(natpop), sep="")
natpop[[1]]         <- NULL
data[,     pop1990 := NULL]
natl[,natl_pop1990 := NULL]

data                <- merge(merge(data, statepop, by="state"), natl, by="date")
# data$subset         <- ( data$date >= as.yearmon("1992-02") & data$date <= as.yearmon("2015-07") )
data$subset         <- ( data$date >= as.yearmon("1992-02") & data$date <= as.yearmon("2014-02") )
data$subset48       <- data$subset & with(data, state != "Alaska" & state != "Hawaii" )
data$subsetWC       <- with(data, subset48 & date <= "May 2014")
data$datenum        <- as.Date(data$date)
data                <- pdata.frame(data, index=c("state", "date") )

rm(statepop)

# want rigs per million people
data$rigs_pop90          <- with(data, rigs_land        / (pop1990                  / 10^6) )
data$rigs_pop00          <- with(data, rigs_land        / (pop2000                  / 10^6) )
data$rigs_pop10          <- with(data, rigs_land        / (pop2010                  / 10^6) )
data$rigs_popinterp      <- with(data, rigs_land        / (pop                      / 10^6) )
data$rigs_poplau         <- with(data, rigs_land        / (lau_pop                  / 10^6) )
data$rigs_poplaglau      <- with(data, rigs_land        / (lag(lau_pop)             / 10^6) )

# employment per capita
data$emp_pop90           <- with(data, ces_emp_nsa      / (pop1990                  / 10^3) )
data$emp_pop00           <- with(data, ces_emp_nsa      / (pop2000                  / 10^3) )
data$emp_pop10           <- with(data, ces_emp_nsa      / (pop2010                  / 10^3) )
data$emp_popinterp       <- with(data, ces_emp_nsa      / (pop                      / 10^3) )
data$emp_poplau          <- with(data, ces_emp_nsa      / (lau_pop                  / 10^3) )
data$emp_poplaglau       <- with(data, ces_emp_nsa      / (lag(lau_pop)             / 10^3) )

data$emp_popsa90           <- with(data, ces_emp_sa      / (pop1990                  / 10^3) )
data$emp_popsa00           <- with(data, ces_emp_sa      / (pop2000                  / 10^3) )
data$emp_popsa10           <- with(data, ces_emp_sa      / (pop2010                  / 10^3) )
data$emp_popsainterp       <- with(data, ces_emp_sa      / (pop                      / 10^3) )
data$emp_popsalau          <- with(data, ces_emp_sa      / (lau_pop                  / 10^3) )
data$emp_popsalaglau       <- with(data, ces_emp_sa      / (lag(lau_pop)             / 10^3) )

# Employment per capita
data$natl_emp_pop90      <- with(data, natl_ces_emp_nsa / (natpop[1,"natl_pop1990"] / 10^3) )
data$natl_emp_pop00      <- with(data, natl_ces_emp_nsa / (natpop[1,"natl_pop2000"] / 10^3) )
data$natl_emp_pop10      <- with(data, natl_ces_emp_nsa / (natpop[1,"natl_pop2010"] / 10^3) )
data$natl_emp_popinterp  <- with(data, natl_ces_emp_nsa / (natl_population          / 10^3) )
data$natl_emp_poplau     <- with(data, natl_ces_emp_nsa / (natl_cps_pop             / 10^0) )
data$natl_emp_poplaglau  <- with(data, natl_ces_emp_nsa / (lag(natl_cps_pop)        / 10^0) )

data$natl_emp_popsa90      <- with(data, natl_ces_emp_sa / (natpop[1,"natl_pop1990"] / 10^3) )
data$natl_emp_popsa00      <- with(data, natl_ces_emp_sa / (natpop[1,"natl_pop2000"] / 10^3) )
data$natl_emp_popsa10      <- with(data, natl_ces_emp_sa / (natpop[1,"natl_pop2010"] / 10^3) )
data$natl_emp_popsainterp  <- with(data, natl_ces_emp_sa / (natl_population          / 10^3) )
data$natl_emp_popsalau     <- with(data, natl_ces_emp_sa / (natl_cps_pop             / 10^0) )
data$natl_emp_popsalaglau  <- with(data, natl_ces_emp_sa / (lag(natl_cps_pop)        / 10^0) )


data$post2008       <- as.yearmon(data$date) >= as.yearmon("2008-01")
data$PRErigs_pop00  <- as.numeric(!data$post2008) * data$rigs_pop00
data$POSTrigs_pop00 <- as.numeric( data$post2008) * data$rigs_pop00

rm(natpop)

# ## Check on magnitudes
# x <- melt(data, id.vars=c("state","date", "statemonth"), na.rm=T)
# x$value <- as.numeric(x$value)
# dcast(x, state ~ variable, mean, na.rm=T,
#       subset=.(variable %in% c("rigs_pop00","rigs_poplau","emp_pop90","emp_poplau","natl_emp_pop90","natl_emp_poplau","pop2000","lau_pop")))
