


.pnbd_Sim_DynCov_gen.calc_adj_m <- function(covariates.dynamic, cov.gammas, covariate.names, upper.date, lower.date){
  
  all.cov.mats <- lapply(X=covariate.names, FUN=function(name){
      acast(data = covariates.dynamic[Cov.Date >= lower.date & Cov.Date <= upper.date], 
            formula = Id~Cov.Date, value.var = name)} )
  adj.m <- exp(Reduce("+", Map("*", all.cov.mats, cov.gammas)))
  
  return(adj.m)
}


## Lifetime process ----
## model the probability that a customer with a given transaction history is still "alive" 
.lifetimecdf <- function(omega, lifetime.adj.m, startdate, date.holdout.end, date.estimation.start){
  # Calculate the cdf of omega from startdate to enddate with the covariates in lifetime.adj.m 
  
  # omega: length of lifetime in weeks, can be a vector!
  
  # lifetime.adj.m = lifetime.adj.m * lambda0 !!

  # Make sure we deal with a matrix even if we just look at one customer
  if (is.null(dim(lifetime.adj.m))){                ## Converts lifetime.adj.m to matrix form if it's not one
    names<-rownames(lifetime.adj.m)
    lifetime.adj.m<-matrix(lifetime.adj.m, nrow = 1)
    colnames(lifetime.adj.m)<-names
  }
  
  
  ## omega - time from the onset of customer's transactions (Tr.#1) until the end of observation period

  # Calculate domega!
  d.time <- interval(startdate, ceiling_date(startdate + seconds(1), "week"))
  d.omega <- as.numeric(d.time, unit="week") # length of the interval (measurement unit)
  
  deathdate <- startdate + days(round(omega * 7))
  

  # We use floordate and the ceiling data for cal.start -> this is to make sure in every last special case the number of periods is taken out correctly
  time.to.startdate <- interval(start = ceiling_date(date.estimation.start+ seconds(1), "week"), 
                               end = floor_date(startdate, "week"))
  time.startdate.death <- interval(start = ceiling_date(startdate+ seconds(1), "week"), 
                                end = floor_date(deathdate, "week"))
  
  # Number of periods (starting from the beginning) until the first purchase per customer
  positionstartdate <- as.numeric(time.to.startdate, units="week") + 2
  
  komega <- as.numeric(time.startdate.death, units="week") + 2 ## because we count both the very first and the very last intervals
  
  #Number of intervals from our start to end of estimation period (induced by omega). What position in the lifetime.adj.m matrix is this?
  positionenddate = positionstartdate + komega - 1 
  
  
  
  
  
  #last part
  last = lifetime.adj.m[ , positionenddate]*(omega - d.omega - as.numeric(komega >= 2)*(komega - 2))

  middle = 0
  # Middle part (2: komega-1):
  if (komega >= 3){
    middle = sum(lifetime.adj.m[,(positionstartdate + 1): (positionenddate - 1) ])
  }
    
  # First part
  first = lifetime.adj.m[ , positionstartdate]* d.omega
  
  cdfval = unlist(unname(1-exp(-(first + middle + last))))
  
  return(cdfval)
}




pnbd_Simulation_Dynamic<- function(lifetime.adj.m.tmp, transaction.adj.m.tmp, params, trans.gammas, life.gammas, covariates.dynamic.names, date.estimation.start, date.holdout.end){
 # We need a special form of walk now (same as in expectation functions)
  
  # startdate = beginning of the universe  
  # enddate = end of the universe  
  
  r = params[1]
  alpha = params[2]
  s = params[3]
  beta = params[4]
  
  
  # Step 1: Simulate from Gamma (heterogenity) distributions
  
  #
  #lambda0=rgamma(1,shape=r, rate=alpha)
  #mu0=rgamma(1,shape=s, rate=beta)
  
  lambda0 = rgamma(1,shape = r,scale = 1/alpha)
  mu0 = rgamma(1,shape = s, scale = 1/beta)
  
  
  # Now add lambda0 and mu0 to get lambda(t) and mu(t)!!
  lifetime.adj.m <- lifetime.adj.m.tmp * mu0
  transaction.adj.m <- transaction.adj.m.tmp * lambda0
  
  
  
  # Step 2: During the whole period simulate from a nonhomogeneous Poisson process in [startuniverse, enduniverse]=[firstcovariate, lastcovariate] -------

  
  # Find lambda that bounds lambda(t) for all t:
  lambdabound = max(transaction.adj.m)
  lambdabound = lambdabound +0.5
  
  # t is always the number of weeks we need to add to startdate!

  t = 0 
  k = 1
  
  T= as.integer(difftime(date.holdout.end, date.estimation.start, units = "weeks"))
  S <- c() # S are event times, when a transaction occurs!
  S[1] = 0 # First element of S is always 0! (which is just the first transaction), not anymore!
  
  
  while (t < T){
    rand = runif(1)
    t = t - log(rand)/lambdabound
    
    if (t >= T){  
      break
    }
    
    rand2 = runif(1)
    
    # Get lambda(t)=lambdat
    # - - - - - - - - - - - - - - - - - - - - - 
    datet <- date.estimation.start + days(round(t * 7))
    
    intervalt <- interval(start = ceiling_date(date.estimation.start, "week"), end = floor_date(datet+seconds(1), "week"))
    
    # Number of periods (starting from the beginning) until the first purchase per customer
    kt <- as.numeric(intervalt, units = "week") + 2
    
    lambdat <- transaction.adj.m[kt] # Awesome :)
    # - - - - - - - - - - - - - - - - - - - - - 
    if (rand2 <= lambdat/lambdabound){
      
      k = k + 1
      S[k] = t #S(j)=t_j, time since beginning of the universe of transaction j 
    }
  }
  
  
  if(length(S) > 0){  # If there is at least one transaction
  
  datestmp<-date.estimation.start + days(round(S * 7)) # Start at the beginning of the universe
  # This gives us transaction dates for the whole period! (the whole universe so to speak)
  
  
  # Step 3: Simulate from lifetime distribution, given the first purchase (time 0) -------
  
  # Simulate from uniform distribution on [0,1]
  
  # Use floordate and the ceiling data for cal.start -> this is to make sure in every last special case the number of periods is taken out correctly
  time.startdate.enddate <- interval(start = datestmp[1], end = date.holdout.end) #S[1] is the time of the first transaction, 
  
  # Number of periods (starting from the beginning) until the first purchase per customer
  maxval <- as.numeric(time.startdate.enddate, units = "week")
  
  U = runif(1)
  # Find omega such that U-F(omega)=0
  fun <- function (omega) (.lifetimecdf(omega, lifetime.adj.m, datestmp[1], date.holdout.end, date.estimation.start) - U)^2
  # Idea: omega can be arbitrarily large in theory, but of course in practice it only can be as long as we have covariates
  # So the idea to search for the solution of the equation below between 0 to the farthest point at which we still have covariates.
  opt <- optim(par = 0, fun, lower = 0.1, upper = maxval, method="L-BFGS-B", hessian = F)
  omega <- opt$par
  

  # So far we have defined the very first transaction of the customer, which is transaction 0!
  death <-  S[1] + omega
  deathdate <- datestmp[1] + days(round(omega * 7))
  

  # Step 4: Cut all transactions that happen above the lifetime of a given customer -------
  
  # Cut S:
  Scut<-S[S < death]
  dates<-datestmp[datestmp < deathdate]
  
  print(U)
  print(omega)
  print(death)
  print(S)
  print(Scut)
  
  return(list(dates, omega))
  
  
}else{
  
  
  return(list(NaN, NaN))  
  
  
  }
  
  
}



