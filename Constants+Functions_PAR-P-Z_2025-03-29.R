# Constants & Functions for PAR-P-Z model 2025-01-16

# CONSTANTS ********************************************************************

PAR0 = 600 # surface irradiance
eps0 <- 0.0213  # Light extinction by pure water, m-1
epsDOC = 0.0514  # DOC light extinction coef, m2 g-1
epsP <- 0.0177   # Phytoplankton, m2 (mg chl)-1
CChl = 60 # C:Chl mass ratio
CC.chl = 60 # Chl carrying capacity = rmax * shading coef / Competition coef
CC.C = CC.chl*CChl # C carrying capacity
ghalf = 0.005*CC.C # nominal is 0.005*CC.C
q = 1 # exponent in grazing term
f = 0.001  # refuge exchange parameter (consider flow from upstream lake)

# Fitted constants to 2024 data all 3 lakes
epsDOCw = 0.1513  # DOC + water light extinction, mean of 3 lakes 2024
rmax = 0.456942350422245  # ML estimate
hmaxML = 0.155761628853678  # ML estimate
hmax = 1.557616   # eyeball value that makes G steep enough (could adjust H instead)
fitsigma = 6.356  # ML sigma in Chl units

# FUNCTIONS **********************************************************************

# Function returns one-step negative log likelihood for one-step prediction of chl
# Input: parameters to be estimated, A as chl, Zt m, C is DOC mg/L, N is vector length
OneStep = function(pars,P0vec,A0vec,A1vec,Ztvec,kNC0,ZBvec,N) {   
  # unpack parameters
  epar = exp(pars)
  rmax = epar[1]
  hmax = epar[2]
  sigma = epar[3]
  
  # simulate
  Ahat = rep(0,N)
  for(i in 1:N) {   # loop over A0 vector
    Chl = A0vec[i]
    A = Chl*CChl  # phyto biomass as C
    ZT = Ztvec[i]
    H = ZBvec[i]
    # break the step into 10 steps
    dt = 0.1
    x0 = A  # C units in loop, except for shading
    Peffect = P0vec[i]/5  # scale P effect to the loading at max chl in prior expts
    CC.CP = CC.C*Peffect  # adjust carrying capacity for P load
    for(j in 1:10) {
      # light extinction epilimnion
      x.chl = x0/CChl # convert C to Chl for extinction calculation
      eps.pool = kNC0[i] + epsP*x.chl 
      PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
      # phytoplankton dynamics
      cc = rmax*(PARmean/PAR0)/CC.CP # competition parameter adjusted for P load
      #G = ( hmax*H*x0/(ghalf + x0) ) # Michaelis-Menten grazing with q = 1
      G = ( hmax*H*x0^q/(ghalf^q + x0^q) ) # type 3 grazing for handling & attack setup
      refuge = f*(CC.CP - x0)
      dx = rmax*(PARmean/PAR0)*x0 - cc*x0*x0 - G + refuge
      x1 = x0 + (dx*dt)
      x0 = x1
    }
    Ahat[i] = x0/CChl  # save Chl for comparison to data
  }
  #devmodel = Ahat-A0vec
  #devobs = A1vec - A0vec
  dev = Ahat - A1vec # compare 1 step prediction to observation
  #dev = devmodel - devobs # compare rates
  #NLL = sum(dev*dev) # least squares, sigma estimate will be bad
  NLL = 0.5*N*log(2*pi) + 0.5*N*log(sigma*sigma) + 
    (sum(dev*dev)/(2*sigma*sigma) )
  return(NLL)
}  # end loss function

# ================================================================================

# Function returns one-step prediction of chlorophyll, A
# Input: parameters to be estimated, A as chl, Zt m, C is DOC mg/L, N is vector length
Yhat = function(pars,P0vec,A0vec,A1vec,Ztvec,kNC0,ZBvec,N) {   
  # unpack parameters
  epar = exp(pars)
  rmax = epar[1]
  hmax = epar[2]
  sigma = epar[3]
  
  # simulate
  Ahat = rep(0,N)
  for(i in 1:N) {   # loop over A0 vector
    Chl = A0vec[i]
    A = Chl*CChl  # phyto biomass as C
    ZT = Ztvec[i]
    H = ZBvec[i]
    # break the step into 10 steps
    dt = 0.1
    x0 = A  # C units in loop, except for shading
    Peffect = P0vec[i]/5    # scale P effect to the loading at max chl in prior expts
    CC.CP = CC.C*Peffect  # adjust carrying capacity for P load
    for(j in 1:10) {
      # light extinction epilimnion
      x.chl = x0/CChl # convert C to Chl for extinction calculation
      kNC = kNC0[i]
      eps.pool = kNC + epsP*x.chl 
      PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
      # phytoplankton dynamics
      cc = rmax*(PARmean/PAR0)/CC.CP # competition parameter adjusted for P load
      #G = ( hmax*H*x0/(ghalf + x0) ) # Michaelis-Menten grazing with q = 1
      G = ( hmax*H*x0^q/(ghalf^q + x0^q) ) # type 3 grazing for handling & attack setup
      refuge = f*(CC.CP - x0)
      dx = rmax*(PARmean/PAR0)*x0 - cc*x0*x0 - G + refuge
      x1 = x0 + (dx*dt)
      x0 = x1
    }
    Ahat[i] = x0/CChl  # save Chl for comparison to data
  }
  return(Ahat)
}  # end one-step projection

# =================================================================================

# general phytoplankton growth & grazing function
# Phytoplankton growth function
# HERE THE DRIVERS ARE SINGLE VALUES TO COMPARE RATES VERSUS A
grograze = function(kNC,p0,ZT,Zb,rmax,hmax) { 
  peffect = p0/5  # p effect on CC for test case; load / saturating load
  gradlen = 1000 # 200 is nominal
  Agrad = seq(0,1.1*CC.C*peffect,length.out=gradlen)  # C units
  Agro = rep(0,gradlen)
  for(i in 1:gradlen) {   # loop for bits that PARmean updates
    Chl = Agrad[i]/CChl  # convert to chlorophyll to compute shading
    eps.pool = kNC + epsP*Chl 
    PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
    CC.CP = CC.C*peffect  # carrying capacity adjusted for P load
    cc = rmax*(PARmean/PAR0)/(CC.CP) # competition parameter
    refuge = f*(CC.CP - Agrad[i])
    Agro[i] = rmax*(PARmean/PAR0)*Agrad[i] - cc*Agrad[i]^2 + refuge 
  }
  # Grazing loss function - OK as vector calculation
  G = ( hmax*Zb*Agrad/(ghalf + Agrad) )  # equation for M-M case
  #G = (hmax*Zb*Agrad^q/(ghalf^q + Agrad^q) )  # equation for sigmoid case
  # Net growth
  Net = Agro - G
  # Equilibria for plotting
  srate = sign(Agro-G)
  dsrate = diff(srate)
  ieq = which(dsrate!=0)
  Aeq.C = Agrad[ieq]
  # List to return
  outlist = list(Agrad,Agro,G,Net,Aeq.C)
  return(outlist)
}

# ===============================================================================

# Function for dA/dt deterministic
# inputs are scalar values
dAdet = function(Chl,kNC,p0,ZT,Zb,rmax,hmax) { 
  peffect = p0/5  # p effect on CC for test case; load / saturating load
  AC = Chl*CChl  # update in carbon units
  eps.pool = kNC + epsP*Chl 
  PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
  CC.CP = CC.C*peffect  # carrying capacity adjusted for P load
  cc = rmax*(PARmean/PAR0)/(CC.CP) # competition parameter
  refuge = f*(CC.CP - AC)
  # deterministic growth rate
  Agro = rmax*(PARmean/PAR0)*AC - cc*AC^2 + refuge 
  # Grazing loss 
  G = ( hmax*Zb*AC/(ghalf + AC) )  # equation for M-M case
  #G = (hmax*Zb*Agrad^q/(ghalf^q + Agrad^q) )  # equation for sigmoid case
  # Net growth with additive noise
  NetChl = (Agro - G)/CChl  # increment converted to Chl units
  return(NetChl)
}

# ===============================================================================

# Function for dA/dt with additive noise
# inputs are scalar values
dAdt = function(Chl,kNC,p0,ZT,Zb,rmax,hmax,shock,sigma,dt,dtn) { 
  peffect = p0/5  # p effect on CC for test case; load / saturating load
  AC = Chl*CChl  # update in carbon units
  eps.pool = kNC + epsP*Chl 
  PARmean = (PAR0/ZT)*(1/eps.pool)*(1 - exp(-eps.pool*ZT))
  CC.CP = CC.C*peffect  # carrying capacity adjusted for P load
  cc = rmax*(PARmean/PAR0)/(CC.CP) # competition parameter
  refuge = f*(CC.CP - AC)
  # deterministic growth rate
  Agro = rmax*(PARmean/PAR0)*AC - cc*AC^2 + refuge 
  # Grazing loss 
  G = ( hmax*Zb*AC/(ghalf + AC) )  # equation for M-M case
  #G = (hmax*Zb*Agrad^q/(ghalf^q + Agrad^q) )  # equation for sigmoid case
  # Net growth with additive noise
  NetC = (Agro - G)*dt + sigma*shock*dtn
  Chlincr = Chl + (NetC/CChl)  # convert increment to Chl units
  Chl1 = max(Chlincr,0.1)  # 0.1 is minimum Chl1 to prevent negative values due to shocks
  return(Chl1)
}

# ==============================================================================

# Functions for critical criterion
#dAdet = function(Chl,kNC,p0,ZT,Zb,rmax,hmax)
crit.pload = function(par,cvec) {  # par is log units
  # parameters
  epar = exp(par)
  Chl = epar[1]
  p0 = epar[2]
  # constants
  kNC = cvec[1]
  ZT = cvec[2]
  ZB = cvec[3]
  rmax = cvec[4]
  hmax = cvec[5]
  rate = dAdet(Chl,kNC,p0,ZT,ZB,rmax,hmax) # deterministic net change in Chl units
  Aup = Chl*1.001
  Adn = Chl*0.999
  rate.up = dAdet(Aup,kNC,p0,RZT,RZbiom,rmax,hmax)
  rate.dn = dAdet(Adn,kNC,p0,RZT,RZbiom,rmax,hmax)
  drate = ( (rate.up-rate)+(rate-rate.dn) )/2
  dA = (Aup-Adn)/2
  dratedA = drate/dA
  loss = rate^2 + dratedA^2 # loss function to be minimized
  return(loss)
}
# ===============================================================================
crit.kNC = function(par,cvec) {
  # parameters
  Chl = par[1]
  kNC = par[2]
  # constants
  p0 = cvec[1]
  ZT = cvec[2]
  ZB = cvec[3]
  rmax = cvec[4]
  hmax = cvec[5]
  rate = dAdet(Chl,kNC,p0,ZT,ZB,rmax,hmax) # deterministic net change in Chl units
  Aup = Chl*1.001
  Adn = Chl*0.999
  rate.up = dAdet(Aup,kNC,p0,RZT,RZbiom,rmax,hmax)
  rate.dn = dAdet(Adn,kNC,p0,RZT,RZbiom,rmax,hmax)
  drate = ( (rate.up-rate)+(rate-rate.dn) )/2
  dA = (Aup-Adn)/2
  dratedA = drate/dA
  loss = rate^2 + dratedA^2 # loss function to be minimized
  return(loss)
}

# =================================================================================

relax = function(acf.x) {   # relaxation time
  xx = acf.x$lag
  yy = acf.x$acf
  y=log(yy)
  lmy = lm(y ~ xx)
  c = unname(lmy$coefficients[2])
  TR = 1/abs(c)
  return(TR)
}

# ==============================================================================

DLM <- function(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat) {
  
  # Online algorithm for Dynamic linear regression 
  # Copyright 2016 by Stephen R. Carpenter
  
  # Description and definitions:
  
  # Observation equation is
  # Y_t = F_t'*theta_t + eta_t where
  # Y_t is the prediction
  # F_t is a vector of predictors at the beginning of the time step
  # theta_t is the parameter vector
  # eta_t is an individual prediction error
  
  # System equation is:
  # theta_t = theta_t-1 + omega_t
  # where theta is defined above and omega_t is an individual process error
  
  # Inputs to the function are:
  # delta, the discount factor
  # n.gamma, the initial number of observations (usually 1)
  # d.gamma, the initial shape parameter for prediction errors
  #  (prior estimate of prediction variance = d.gamma / n.gamma)
  # mvec, the initial guess of regression coefficients
  # Cpar, the initial guess of the covariance matrix of regression coefficients
  # Yvec, the vector of the observed response variate
  # Fmat, the matrix of predictors
  
  # Outputs are:
  # predix, the one-step-ahead predictions of the response variate
  # varpredix, the prediction variance at start of time step before error is measured
  # pars, the updated parameter estimates using the most recent prediction error
  # parvar, the variances of the parameters
  # Svec, the update (after error is measured within a time step) of varpredix
  
  # Updating follows the equations on p. 176-179 of Carpenter 2003,
  # Regime Shifts in Lake Ecosystems: Pattern and Variation
  
  # Determine constants
  npar <- length(mvec)
  Nobs <- length(Yvec)
  S0 <- d.gamma/n.gamma
  
  # Set up vectors to hold results
  predix <- rep(0,Nobs)
  varpredix <- rep(0,Nobs)
  Svec = rep(0,Nobs)
  pars <- matrix(0,nrow=Nobs,ncol=npar)
  parvar = matrix(0,nrow=Nobs,ncol=npar)
  
  for(i in 1:Nobs)  {  #Start DLM loop
    # Generate predictions
    Fvec <- Fmat[i,] # vector of predictors
    predix[i] <- sum(Fvec*mvec)
    # Compute error and update estimates
    error <- Yvec[i]-predix[i]
    Rmat <- Cpar/delta
    varpredix[i] <- (t(Fvec) %*% Rmat %*% Fvec) + S0
    n.gamma <- (delta*n.gamma)+1
    d.gamma <- (delta*d.gamma)+(S0*error*error/varpredix[i])
    S1 <- d.gamma/n.gamma
    Svec[i] = S1  # save updated variance
    Avec <- (Rmat %*% Fvec)/varpredix[i]
    mvec <- mvec + (Avec*error)
    pars[i,] <- mvec
    Cpar <- (S1/S0)*(Rmat - (Avec %*% t(Avec))*varpredix[i])
    # Disallow negative variances on the diagonal
    for(idiag in 1:npar) {
      Cpar[idiag,idiag] <- max(0,Cpar[idiag,idiag])
    }
    parvar[i,] = diag(Cpar)
    S0 <- S1 # roll over S
  } # End DLM loop
  
  DLM.out <- list(predix,varpredix,pars,parvar,Svec)
  return(DLM.out)
} # END DLM FUNCTION

# ================================================================================

# ==============================================================================

# DLM with forecast
# additional input is kh, integer forecast horizon

DLMfore <- function(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat,kh) {
  
  # Online algorithm for Dynamic linear regression 
  # Copyright 2016 by Stephen R. Carpenter, forecast added 2025
  
  # Description and definitions:
  
  # Observation equation is
  # Y_t = F_t'*theta_t + eta_t where
  # Y_t is the prediction
  # F_t is a vector of predictors at the beginning of the time step
  # theta_t is the parameter vector
  # eta_t is an individual prediction error
  
  # System equation is:
  # theta_t = theta_t-1 + omega_t
  # where theta is defined above and omega_t is an individual process error
  
  # Inputs to the function are:
  # delta, the discount factor
  # n.gamma, the initial number of observations (usually 1)
  # d.gamma, the initial shape parameter for prediction errors
  #  (prior estimate of prediction variance = d.gamma / n.gamma)
  # mvec, the initial guess of regression coefficients
  # Cpar, the initial guess of the covariance matrix of regression coefficients
  # Yvec, the vector of the observed response variate
  # Fmat, the matrix of predictors
  
  # Outputs are:
  # predix, the one-step-ahead predictions of the response variate
  # varpredix, the prediction variance at start of time step before error is measured
  # pars, the updated parameter estimates using the most recent prediction error
  # parvar, the variances of the parameters
  # Svec, the update (after error is measured within a time step) of varpredix
  
  # Updating follows the equations on p. 176-179 of Carpenter 2003,
  # Regime Shifts in Lake Ecosystems: Pattern and Variation
  
  # Determine constants
  npar <- length(mvec)
  Nobs <- length(Yvec)
  S0 <- d.gamma/n.gamma
  
  # Set up vectors to hold results
  predix <- rep(0,Nobs)
  varpredix <- rep(0,Nobs)
  Svec = rep(0,Nobs)
  pars <- matrix(0,nrow=Nobs,ncol=npar)
  parvar = matrix(0,nrow=Nobs,ncol=npar)
  fcast = matrix(0,nrow=Nobs,ncol=4)  # columns: kh-ahead forecast, obs, err, Q
  
  for(i in 1:Nobs)  {  #Start DLM loop
    # Generate predictions
    Fvec <- Fmat[i,] # vector of predictors
    predix[i] <- sum(Fvec*mvec)
    # Compute error and update estimates
    error <- Yvec[i]-predix[i]
    Rmat <- Cpar/delta
    varpredix[i] <- (t(Fvec) %*% Rmat %*% Fvec) + S0
    n.gamma <- (delta*n.gamma)+1
    d.gamma <- (delta*d.gamma)+(S0*error*error/varpredix[i])
    S1 <- d.gamma/n.gamma
    Svec[i] = S1  # save updated variance
    Avec <- (Rmat %*% Fvec)/varpredix[i]
    mvec <- mvec + (Avec*error)
    pars[i,] <- mvec
    Cpar <- (S1/S0)*(Rmat - (Avec %*% t(Avec))*varpredix[i])
    # Disallow negative variances on the diagonal
    for(idiag in 1:npar) {
      Cpar[idiag,idiag] <- max(0,Cpar[idiag,idiag])
    }
    parvar[i,] = diag(Cpar)
    if(i < (Nobs-kh)) {   # kh-ahead forecast
      ikh = i+kh
      forekh = sum(Fmat[ikh,]*mvec)
      efore = forekh - Yvec[ikh]
      Qfore = t(Fmat[ikh,])%*%Rmat%*%Fmat[ikh,]+(kh*S1)
      fcast[i,] = c(forekh,Yvec[ikh],efore,Qfore)
    } # end kh-ahead forecast from timestep i

    S0 <- S1 # roll over S
  } # End DLM loop
  
  DLM.out <- list(predix,varpredix,pars,parvar,Svec,fcast)
  return(DLM.out)
} # END DLM FUNCTION

# ================================================================================

# ONLINE DYNAMIC LINEAR MODEL (DLM) ESTIMATION FOR MULTIVARIATE CASE
# from MAR1_No-intercept-onlineDLM+bootstraps_Peter2024_2024-12-22.R

MDLM <- function(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat) {
  
  # Online algorithm for Dynamic linear regression 
  # Copyright 2016 by Stephen R. Carpenter
  
  # Description and definitions:
  
  # Observation equation is
  # Y_t = F_t'*theta_t + eta_t where
  # Y_t is the prediction
  # F_t is a vector of predictors at the beginning of the time step
  # theta_t is the parameter vector
  # eta_t is an individual observation error
  
  # System equation is:
  # theta_t = theta_t-1 + omega_t
  # where theta is defined above and omega_t is an individual process error
  
  # Inputs to the function are:
  # delta, the discount factor
  # n.gamma, the initial number of observations (usually 1)
  # d.gamma, the initial shape parameter for prediction errors
  #  (prior estimate of prediction variance = d.gamma / n.gamma)
  # mvec, the initial guess of regression coefficients
  # Cpar, the initial guess of the covariance matrix of regression coefficients
  # Yvec, the vector of the observed response variate
  # Fmat, the matrix of predictors
  
  # Outputs are:
  # predix, the one-step-ahead predictions of the response variate
  # varpredix, the prediction variance at start of time step before error is measured
  # pars, the updated parameter estimates using the most recent prediction error
  # parvar, the variances of the parameters
  # Svec, the update (after error is measured within a time step) of varpredix
  
  # Updating follows the equations on p. 176-179 of Carpenter 2003,
  # Regime Shifts in Lake Ecosystems: Pattern and Variation
  
  # Determine constants
  npar <- length(mvec)
  Nobs <- length(Yvec)
  S0 <- d.gamma/n.gamma
  
  # Set up vectors to hold results
  predix <- rep(0,Nobs)
  varpredix <- rep(0,Nobs)
  Svec = rep(0,Nobs)
  pars <- matrix(0,nrow=Nobs,ncol=npar)
  parvar = matrix(0,nrow=Nobs,ncol=npar)
  
  for(i in 1:Nobs)  {  #Start DLM loop
    # Generate predictions
    Fvec <- Fmat[i,] # vector of predictors
    predix[i] <- sum(Fvec*mvec)
    # Compute error and update estimates
    error <- Yvec[i]-predix[i]
    Rmat <- Cpar/delta
    varpredix[i] <- (t(Fvec) %*% Rmat %*% Fvec) + S0
    n.gamma <- (delta*n.gamma)+1
    d.gamma <- (delta*d.gamma)+(S0*error*error/varpredix[i])
    S1 <- d.gamma/n.gamma
    Svec[i] = S1  # save updated variance
    Avec <- (Rmat %*% Fvec)/varpredix[i]
    mvec <- mvec + (Avec*error)
    pars[i,] <- mvec
    Cpar <- (S1/S0)*(Rmat - (Avec %*% t(Avec))*varpredix[i])
    # Disallow negative variances on the diagonal
    for(idiag in 1:npar) {
      Cpar[idiag,idiag] <- max(0,Cpar[idiag,idiag])
    }
    parvar[i,] = diag(Cpar)
    S0 <- S1 # roll over S
  } # End DLM loop
  
  DLM.out <- list(predix,varpredix,pars,parvar,Svec)
  return(DLM.out)
} # END DLM FUNCTION

