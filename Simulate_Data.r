
library(CLVTools)
library(ggplot2)
library(lubridate)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)

source("functions/pnbd_Simulation_Dynamic.R")

#Set a seed
set.seed(12345)

# Load the exogenous covariates (they are not generated!)
# Currently there are covariates for 11000 customers
load("data/covariates.Rdata")

# Number of customers to be simulated
n = 2000 # max: 4206 customers due to exogenous contextual factors

# Duration of estimation period in weeks
estimation.duration <- 104

## Parameters ----
##Specify parameters
# Base parameters
r = 1.5 # Homogeneity purchase process
alpha = 170 # Shape parameter purchase process
s = 0.5 # Homogeneity attrition process
beta = 8 # Shape parameter attrition process

# Transaction contextual factor parameters
gamma1=1 
#gamma2=0.5
#gamma3=0

# Lifetime contextual factor parameters
gamma4=1
#gamma5=0.5
#gamma6=0

covariate.names <- c("direct.marketing")
#covariate.names <- c("direct.marketing", "high.season")

#combine parameters
params = c(r, alpha, s, beta)
trans.gammas = c(gamma1)
life.gammas = c(gamma4)
#trans.gammas = c(gamma1, gamma2)
#life.gammas = c(gamma4, gamma5)
#trans.gammas = c(gamma1, gamma2, gamma3)
#life.gammas = c(gamma4, gamma5, gamma6)



# Get intervals of important dates
cal.start.date      <- min(covariates.dynamic$Cov.Date)
hold.end.date       <- max(covariates.dynamic$Cov.Date)


## Preparations ----
# Covariate table

# Transaction process
transaction.adj.m.all <-  .pnbd_Sim_DynCov_gen.calc_adj_m(covariates.dynamic=covariates.dynamic, cov.gammas = trans.gammas, 
                                                          covariate.names = covariate.names, 
                                                          lower.date = cal.start.date,
                                                          upper.date = hold.end.date )
# Lifetime process
lifetime.adj.m.all <-     .pnbd_Sim_DynCov_gen.calc_adj_m(covariates.dynamic=covariates.dynamic, cov.gammas = life.gammas, 
                                                                     covariate.names = covariate.names,
                                                                     lower.date = cal.start.date,
                                                                     upper.date = hold.end.date )


## Simulation ----
# Generate the transactions for every customer

datai=list()

#Make sure we use all cores
no.cores<-max(detectCores()-1,1)
cl <- makeCluster(no.cores)
registerDoParallel(cl)

# Generate data for every customer i
datai <- foreach(i=1:n, .packages = c("lubridate", "data.table"))%dopar%{

    # Get individual covariates
     lifetime.adj.m.tmp <- lifetime.adj.m.all[i, ]
     transaction.adj.m.tmp <- transaction.adj.m.all[i, ]
    
    
    # Call the simulation function
    output<- pnbd_Simulation_Dynamic(lifetime.adj.m.tmp,transaction.adj.m.tmp, params , trans.gammas, life.gammas, covariates.dynamic.names, cal.start.date, (cal.start.date+weeks(estimation.duration)))
    
    dates <- output[[1]] # Save dates of simulated transactions
    omg <- output[[2]] # save omega (length of lifetime), just for debugging...
    

    # Use Id's from the dataset directly (important to afterwards filter out the correct covariates)
    return(data.table(Id = rep(rownames(lifetime.adj.m.all)[i], length(dates)), Date=dates, omegas=omg, Price=1) ) 
    
  }

simdata=rbindlist(datai)

stopCluster(cl)

#Make sure that there are no NaNs in omega
simdata[is.nan(simdata$omegas)]

#check the number of customer
length(unique(simdata$Id))

# Save data
covdata <- covariates.dynamic[Id %in% simdata$Id]
save("simdata", "covdata",  file="SimDataCovariates.RData") # With covariates
#save("simdata",  file="SimData.RData") # base model
