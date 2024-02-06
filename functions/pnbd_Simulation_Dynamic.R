

## Function for lifetime process ----
## model the probability that a customer with a given transaction history is still "alive" 
.lifetimecdf <- function(omega, dt.life.i, startdate){
  # Calculate the cdf of omega from startdate to enddate with the covariates in lifetime.adj.m 
  
  # omega: length of lifetime in weeks!
  ## omega - time from the onset of customer's transactions (Tr.#1) until the end of observation period
  # Calculate domega!
  endofintervalone<-dt.life.i[Cov.Date > startdate][Cov.Date==min(Cov.Date)]$Cov.Date
  d.time <- interval(startdate, endofintervalone-days(1))
  d.omega <- as.numeric(d.time, unit="week") # length of the interval (measurement unit)
  
  deathdate <- startdate + days(round(omega * 7))
  intervaldeath<-dt.life.i[Cov.Date <= deathdate][Cov.Date==max(Cov.Date)]$Cov.Date
  
  # |-x-|--|--|--|-x-| => Nr of intervals between the two events (x) + the ones they are included in 
  komega <- as.numeric(nrow(dt.life.i[Cov.Date < intervaldeath & Cov.Date >= endofintervalone]))+2
  
  #last part
  last = dt.life.i[Cov.Date <= deathdate][Cov.Date==max(Cov.Date)]$mu_i*
    (omega - d.omega - as.numeric(komega >= 2)*(komega - 2))
  
  # Middle part (2: komega-1):
  middle = 0
  if (komega >= 3){
    middle = sum(dt.life.i[Cov.Date < intervaldeath & Cov.Date >= endofintervalone]$mu_i)
  }
  
  # First part
  first = dt.life.i[Cov.Date <= startdate][Cov.Date==max(Cov.Date)]$mu_i* d.omega
  
  cdfval = unlist(unname(1-exp(-(first + middle + last))))
  
  return(cdfval)
}




pnbd_Simulation_Dynamic<- function(dt.life.i, dt.trans.i, params, date.estimation.start, date.holdout.end){
  # We need a special form of walk (same as in expectation functions)
  
  # startdate = beginning of the universe  
  # enddate = end of the universe  
  
  r = params[1]
  alpha = params[2]
  s = params[3]
  beta = params[4]
  
  
  ## Step 1: Simulate from Gamma (heterogenity) distributions -----
  
  #purchase rate
  lambda0 = stats::rgamma(1, shape = r,rate = alpha)
  #attrition rate
  mu0 = stats::rgamma(1, shape = s, rate = beta)
  
  # Now add lambda0 and mu0 to get lambda(t) and mu(t)!!
  dt.life.i[, mu_i:=exp.gX.L*mu0]
  dt.trans.i[, lambda_i:=exp.gX.P*lambda0]
  

  ## Step 2: During the whole period simulate from a nonhomogeneous Poisson process in [startuniverse, enduniverse]=[firstcovariate, lastcovariate] -------
  
  # Find lambda that bounds lambda(t) for all t:
  lambdabound = max(dt.trans.i$lambda_i)
  lambdabound = lambdabound +0.5
  
  # t is always the number of weeks we need to add to startdate!
  t = 0 
  k = 1
  
  T.cal= as.integer(difftime(date.holdout.end, date.estimation.start, units = "weeks"))
  
  S <- c() # S are event times, when a transaction occurs!
  S[1] = 0 # First element of S is always 0! (which is just the first transaction), not anymore!
  
  
  while (t < T.cal){
    rand = stats::runif(1)
    t = t - log(rand)/lambdabound
    
    if (t >= T.cal){  
      break
    }
    
    rand2 = stats::runif(1)
    
    # Get lambda(t)=lambdat
    # - - - - - - - - - - - - - - - - - - - - - 
    datet <- date.estimation.start + days(round(t * 7))
  
    lambdat<-dt.trans.i[Cov.Date <= datet][Cov.Date==max(Cov.Date)]$lambda_i
    
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
    
    U = stats::runif(1)

    # Find omega such that U-F(omega)=0
    fun <- function (omega) {((.lifetimecdf(omega=omega,dt.life.i=dt.life.i,startdate=datestmp[1]) - U))^2*1000}

    # Idea: omega can be arbitrarily large in theory, but of course in practice it only can be as long as we have covariates
    # So the idea to search for the solution of the equation below between 0 to the farthest point at which we still have covariates.
    opt <- optim(par = 0, fun, lower = 0, upper = maxval, method="L-BFGS-B", hessian = F)
    omega <- opt$par
    
    # So far we have defined the very first transaction of the customer, which is transaction 0!
    death <-  S[1] + omega
    deathdate <- datestmp[1] + days(round(omega * 7))
    
    
    # Step 4: Cut all transactions that happen above the lifetime of a given customer -------
    # Cut S:
    Scut<-S[S <= death]
    dates<-datestmp[datestmp <= deathdate]
    
    #print(U)
    #print(omega)
    #print(death)
    #print(S)
    #print(Scut)
    
    return(list(dates, omega))
    
  }else{
    return(list(NaN, NaN))  
  }
}



