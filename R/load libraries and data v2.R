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

load("dat/natl-data-2016-08-11.Rdata")
load("dat/state-data-2016-08-11.Rdata")

natl <- natl.dat
data <- state.dat
rm(natl.dat, state.dat)

data$subset         <- ( data$date >= as.yearmon("1992-02") & data$date <= as.yearmon("2014-06") )
data$subset48       <- data$subset & with(data, state != "Alaska" & state != "Hawaii" )
data$subset44       <- data$subset & with(data, !(state %in% c("Alaska", "Hawaii", "Florida","Nevada","California","Arizona" )) )
data$datenum        <- as.Date(data$date)


# want rigs per million people
data[, rigs_pop       := rigs_land      / (pop2000      / 10^6) ]
data[, rigs_pop90     := rigs_land      / (pop1990      / 10^6) ]
data[, rigs_pop00     := rigs_land      / (pop2000      / 10^6) ]
data[, rigs_pop10     := rigs_land      / (pop2010      / 10^6) ]
data[, rigs_popinterp := rigs_land      / (pop_interp   / 10^6) ]
data[, rigs_poplau    := rigs_land      / (lau_pop      / 10^6) ]
data[, rigs_poplaglau := rigs_land      / (lag(lau_pop) / 10^6) ]

# splits in rigs
data[, post2008     := as.yearmon(date) >= as.yearmon("2008-01")]
data[, PRErigs_pop  := as.numeric(!post2008) * rigs_pop]
data[, POSTrigs_pop := as.numeric( post2008) * rigs_pop]

data[natl, OILrigs_pop   := rigs_pop * (natl_rigs_oil/natl_rigs_total)    , on="date"]
data[natl, GASrigs_pop   := rigs_pop * (natl_rigs_gas/natl_rigs_total)    , on="date"]
data[natl, HORZrigs_pop  := rigs_pop * (  natl_rigs_horiz/natl_rigs_total), on="date"]
data[natl, VERTrigs_pop  := rigs_pop * (1-natl_rigs_horiz/natl_rigs_total), on="date"]

data[natl, PREOILrigs_pop    := as.numeric(!post2008) * rigs_pop * (    natl_rigs_oil / natl_rigs_total), on="date"]
data[natl, PREGASrigs_pop    := as.numeric(!post2008) * rigs_pop * (    natl_rigs_gas / natl_rigs_total), on="date"]
data[natl, PREHORZrigs_pop   := as.numeric(!post2008) * rigs_pop * (  natl_rigs_horiz / natl_rigs_total), on="date"]
data[natl, PREVERTrigs_pop   := as.numeric(!post2008) * rigs_pop * (1-natl_rigs_horiz / natl_rigs_total), on="date"]
data[natl, POSTOILrigs_pop   := as.numeric(post2008)  * rigs_pop * (    natl_rigs_oil / natl_rigs_total), on="date"]
data[natl, POSTGASrigs_pop   := as.numeric(post2008)  * rigs_pop * (    natl_rigs_gas / natl_rigs_total), on="date"]
data[natl, POSTHORZrigs_pop  := as.numeric(post2008)  * rigs_pop * (  natl_rigs_horiz / natl_rigs_total), on="date"]
data[natl, POSTVERTrigs_pop  := as.numeric(post2008)  * rigs_pop * (1-natl_rigs_horiz / natl_rigs_total), on="date"]



# jobs per 1000 ppl
data[, emp_pop_sa       := ces_emp_sa / (pop2000      / 10^3) ]
data[, emp_pop_sa90     := ces_emp_sa / (pop1990      / 10^3) ]
data[, emp_pop_sa00     := ces_emp_sa / (pop2000      / 10^3) ]
data[, emp_pop_sa10     := ces_emp_sa / (pop2010      / 10^3) ]
data[, emp_pop_sainterp := ces_emp_sa / (pop_interp   / 10^3) ]
data[, emp_pop_salau    := ces_emp_sa / (lau_pop      / 10^3) ]
data[, emp_pop_salaglau := ces_emp_sa / (lag(lau_pop) / 10^3) ]

data[, emp_pop_nsa       := ces_emp_sa / (pop2000      / 10^3) ]
data[, emp_pop_nsa90     := ces_emp_sa / (pop1990      / 10^3) ]
data[, emp_pop_nsa00     := ces_emp_sa / (pop2000      / 10^3) ]
data[, emp_pop_nsa10     := ces_emp_sa / (pop2010      / 10^3) ]
data[, emp_pop_nsainterp := ces_emp_sa / (pop_interp   / 10^3) ]
data[, emp_pop_nsalau    := ces_emp_sa / (lau_pop      / 10^3) ]
data[, emp_pop_nsalaglau := ces_emp_sa / (lag(lau_pop) / 10^3) ]

# interpolate population
natl$natl_pop_interp <- with(natl, approx(x=date, y=natl_population, xout=date))$y

natl[, natl_emp_pop_sa       := natl_ces_emp_sa/(natl_pop2000         / 10^3)]
natl[, natl_emp_pop_sa90     := natl_ces_emp_sa/(natl_pop1990         / 10^3)]
natl[, natl_emp_pop_sa00     := natl_ces_emp_sa/(natl_pop2000         / 10^3)]
natl[, natl_emp_pop_sa10     := natl_ces_emp_sa/(natl_pop2010         / 10^3)]
natl[, natl_emp_pop_sainterp := natl_ces_emp_sa/(natl_pop_interp      / 10^3)]
natl[, natl_emp_pop_salau    := natl_ces_emp_sa/(natl_cps_pop_sa      / 10^3)]
natl[, natl_emp_pop_salaglau := natl_ces_emp_sa/(lag(natl_cps_pop_sa) / 10^3)]

natl[, natl_emp_pop_nsa       := natl_ces_emp_nsa/(natl_pop2000          / 10^3)]
natl[, natl_emp_pop_nsa90     := natl_ces_emp_nsa/(natl_pop1990          / 10^3)]
natl[, natl_emp_pop_nsa00     := natl_ces_emp_nsa/(natl_pop2000          / 10^3)]
natl[, natl_emp_pop_nsa10     := natl_ces_emp_nsa/(natl_pop2010          / 10^3)]
natl[, natl_emp_pop_nsainterp := natl_ces_emp_nsa/(natl_pop_interp       / 10^3)]
natl[, natl_emp_pop_nsalau    := natl_ces_emp_nsa/(natl_cps_pop_nsa      / 10^3)]
natl[, natl_emp_pop_nsalaglau := natl_ces_emp_nsa/(lag(natl_cps_pop_nsa) / 10^3)]

data <- merge(data, natl, by="date", all=T)[date >= 1987]

data <- pdata.frame(data, index=c("state", "date") )