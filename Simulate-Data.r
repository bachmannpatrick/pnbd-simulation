

#load required packages
library(CLVTools)
library(lubridate)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(stats)

# Source the simulation function
source("functions/pnbd_Simulation_Dynamic.r")

# Set a seed for reproducability
set.seed(123456)

# Load the exogenous covariates (they are not generated!)
load("data/covariates.Rdata")

# Number of customers to be simulated
n = 2000 # max: 11000 customers due to exogenous contextual factors

# Duration of estimation period in weeks
estimation.duration <- 416

## Parameters ----
# Base parameters
r = 0.8 # Homogeneity purchase process
alpha = 4.5 # Shape parameter purchase process
s = 0.5 # Homogeneity attrition process
beta = 40 # Shape parameter attrition process

# Lifetime contextual factor parameters
gamma4=0.6
gamma5=1
gamma6=0.0

# Transaction contextual factor parameters
gamma1=0.2
gamma2=0.5
gamma3=0.0

#Three covariates are available
covariate.names.trans <- c("high.season", "gender")
covariate.names.life <- c("high.season", "gender")
#covariate.names.trans <- c("direct.marketing")
#covariate.names.life <- c("direct.marketing")

#combine parameters
params = c(r, alpha, s, beta)
#trans.gammas = c(gamma1=gamma1)
#life.gammas = c(gamma4=gamma4)
#trans.gammas = c(gamma1=gamma1, gamma2=gamma2)
#life.gammas = c(gamma4=gamma4, gamma5=gamma5)
trans.gammas = c(gamma1=gamma1, gamma2=gamma2)
life.gammas = c(gamma4=gamma4, gamma5=gamma5)

# Get intervals of important dates
cal.start.date      <- min(covariates.dynamic$Cov.Date)
holdout.end.date    <- max(covariates.dynamic$Cov.Date)

## Preparations ----
# Covariate table
covariates.dynamic[ , ( names(trans.gammas) ):= data.table( t(trans.gammas) ) ]
covariates.dynamic[ , ( names(life.gammas) ):= data.table( t(life.gammas) ) ]
covariates.dynamic[, exp.gX.P:= exp(rowSums(covariates.dynamic[,.SD,.SDcols=covariate.names.trans]*covariates.dynamic[,.SD,.SDcols=names(trans.gammas)]))]
covariates.dynamic[, exp.gX.L:= exp(rowSums(covariates.dynamic[,.SD,.SDcols=covariate.names.life]*covariates.dynamic[,.SD,.SDcols=names(life.gammas)]))]

## Simulation ----
# Generate the transactions for every customer
datai = list()

#Make sure we use all cores (minus 1)
no.cores <- max(detectCores()-1,1)
cl <- makeCluster(no.cores)
registerDoParallel(cl)

# Generate data for every customer i
datai <- foreach(i=covariates.dynamic[,unique(Id)][1:n], .packages = c("lubridate", "data.table"))%dopar%{
  
  # Call the simulation function
  output<- pnbd_Simulation_Dynamic(dt.life.i=covariates.dynamic[Id==i,c("Id", "Cov.Date","exp.gX.L")], 
                                   dt.trans.i=covariates.dynamic[Id==i,c("Id", "Cov.Date","exp.gX.P")], 
                                   params=params , date.estimation.start=cal.start.date, 
                                   date.holdout.end = (cal.start.date+weeks(estimation.duration)))
                                   #date.holdout.end = holdout.end.date)
  
  dates <- output[[1]] # Save dates of simulated transactions
  omg <- output[[2]] # save omega (length of lifetime), just for debugging...
  
  # Use Id's from the dataset directly (important to afterwards filter out the correct covariates)
  # Note: We do not simulate prices. All transactions prices are set to 1.
  return(data.table(Id = rep(i, length(dates)), Date=dates, Price=1) )
}

simdata=rbindlist(datai)
stopCluster(cl)

#Make sure that there are no NaNs in omega
simdata[is.nan(simdata$omegas)]

#check the number of customer
length(unique(simdata$Id))

# Save data
#workaround preventing errors for timestamp
# Save data
covdata <- covariates.dynamic[Id %in% simdata$Id]
varibles <- unique(c("Id", "Cov.Date", covariate.names.trans, covariate.names.trans))
covdata <- covariates.dynamic[, ..varibles]
save("simdata", "covdata", "covariate.names.trans", "covariate.names.life", file="SimDataCovariates.RData") # With covariates
