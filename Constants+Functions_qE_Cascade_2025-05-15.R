# Constants & Functions for qE rate experiments Cascade

# CONSTANTS ********************************************************************

#
#  FOOD WEB PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Phytoplankton ************

I0 <- 300  # Surface irradiance, microEinsteins m-1 s-1
# P-I curve parameters, Follows et al.
k_sat <- 0.012 # per microEinsteins m-1 s-1
k_inh <- 0.004 # per microEinsteins m-1 s-1 (nominal 0.004, range 0.001-0.007)
# Light extinction parameters, Carpenter et al. L&O 1998
DOC <- 5  # Assumed DOC conc, g m-3
eps0 <- 0.0213 + 0.0514*DOC  # Baseline, per m (nominal DOC par 0.0514)
epsP <- 0.0177   # Phytoplankton, per m2 (mg phosphorus)-1

# Derived parameter from Follows et al.
Fmax <- ((k_sat + k_inh)/k_sat)*exp(-(k_inh/k_sat)*log(k_inh/(k_sat+k_inh)))

# Information for depth integration
Zmix <- 4  # Mixed layer depth, m
nZ <- 10  # Steps for vertical integration
dZ <- Zmix/nZ
Zvec <- seq(0,Zmix,by=dZ)

rP <- 3  # Phytoplankton growth parameter per unit phosphorus load
Load <- 0.6  # Daily phosphorus load
mP <- 0.1  # Phytoplankton daily mortality 
#sigmaP = 1.08  # SD of additive noise for phytoplankton

# Zooplankton ***********************
Ho <- 4    # Refuge biomass    # 1 in OLD
DH <- 0.5  # Diffusion parameter  # 0.5 in OLD
cHF <- 0.1  # Consumption rate by planktivore
alf <- 0.3  # Conversion efficiency of consumed phytoplankton to zooplankton
cPH <- 0.25  # Consumption rate of phytoplankton by zooplankton
#sigmaH = 0.02  # SD of additive noise for herbivores

# Fish *****************************
qELO <- 1  # First Catchability x Effort
qEHI <- 4  # Second Catchability x Effort
fA <- 2  # Fecundity of adult piscivore (2 in OLD)
cJA <- 0.1  # Density dependent mortality rate of juveniles
cJF <- 0.5  # Consumption of juveniles by planktivores
cFA <- 0.3  # Consumption of planktivores by adult piscivores
vuln <- 80  # Vulnerability coefficient
hide <- 80  # Hiding coefficient
surv <- 0.6  # Overwinter survivorship of adults
Fo <- 200  # Refuge density of planktivores  # 100 in OLD
DF <- 0.09  # Diffusion parameter for planktivores
# additive noise sigmas from
# /2025_RateDependence_kNC_P_qE/RelaxationTimes_dailies_Squeal1_2008-2011.R
sigmaF = 2.1  # planktivores 
sigmaH = 0.02 # zooplankton
sigmaP = 1.08  # chl

A2biom <- 0.2  # Convert A to kg / ha
J2biom <- 0.05  # Convert J to kg / ha
F2biom <- 1  # Convert F to kg / ha

#
# INITIAL CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~
#
Ainit <- 25  # Critical A is 88.3
Finit <- 4
Hinit <- 6
Pinit <- 7

# FUNCTIONS ======================================================================

#
#Phytoplankton Growth ********************************************************
GAMMA <- function(z,Pbar) {
  Iz <- I0*exp(-z*(eps0+epsP*Pbar))
  rate <- (1/Fmax)*(1 - exp(-k_sat*Iz))*exp(-k_inh*Iz)
}  # End Phyto Growth function ***********************************************

# Simulation of one time step *************************************************
FWsim.step <- function(qE,A0,F0,J0,H0,P0,dt,dtZ,nvec) { # Start simulation function
  
  NoiseF <- nvec[1]
  NoiseH <- nvec[2]
  NoiseP <- nvec[3]
  
  # Fish dynamics
  #Arate <- (surv/nint)*J0 - qE*A0 - ((1-surv)/nint)*A0
  Arate <- (surv)*J0 - qE*A0 - ((1-surv))*A0
  Frate <- DF*(Fo-F0) - cFA*F0*A0
  Jpredloss <- (-cJA*J0*A0)-(cJF*vuln*J0*F0/(hide + vuln + cJF*F0) ) # Note this is negative
  #Jrate <- (fA/nint)*A0 + Jpredloss - (surv/nint)*J0 
  Jrate <- (fA)*A0 + Jpredloss - (surv)*J0
  A1 <- A0 + (Arate*dt)   # Update A numerically
  F1 <- F0 + (Frate*dt) + (sigmaF*NoiseF*dtZ)
  J1 <- J0 + (Jrate*dt)
  A1 <- max(A1,0.1)  # Force A1 greater than 0.1
  F1 <- max(F1,0.1)  # Force F greater than 0.1
  J1 <- max(J1,0.1)  # Force J greater than 0.1
  # Zooplankton dynamics
  Hrate <- DH*(Ho-H0) + alf*cPH*H0*P0 - cHF*H0*F0
  H1 <- H0 + (Hrate*dt) + (sigmaH*NoiseH*dtZ)
  H1 <- max(H1,0.1)  # Force H greater than 0.01
  # Phytoplankton dynamics
  Pbar <- P0   # Set P value for vertical integration
  gamvec <- GAMMA(Zvec,Pbar) # vertical sequence of light effect on growth
  gamI <- dZ*sum(gamvec)  # vertically integrated light effect on growth
  Prate <- (rP*Load*gamI*P0) - (mP*P0) - (cPH*H0*P0)  
  P1 <- P0 + (Prate*dt) + (sigmaP*NoiseP*dtZ)
  P1 <- max(P1,0.1)  # Force P greater than 0.1
  
  # Construct list for output
  SimList <- list(A1,F1,J1,H1,P1)
  return(SimList)
  
}  # End one-step simulation function  ************************************

# Simulation by iteration *************************************************
FWsim.iter <- function(qE,A0,F0,J0,H0,P0,dt,dtZ,nint) { # Start simulation function
  # iterate nint times with a time step
  
  noise = rnorm(3*nint)
  noisemat = matrix(noise,nrow=nint,ncol=3)
  
  for (j in 1:nint) {
    nvec = noisemat[j,]
    NoiseF <- nvec[1]
    NoiseH <- nvec[2]
    NoiseP <- nvec[3]
    # Fish dynamics
    #Arate <- (surv/nint)*J0 - qE*A0 - ((1-surv)/nint)*A0
    Arate <- (surv)*J0 - qE*A0 - ((1-surv))*A0
    Frate <- DF*(Fo-F0) - cFA*F0*A0
    Jpredloss <- (-cJA*J0*A0)-(cJF*vuln*J0*F0/(hide + vuln + cJF*F0) ) # Note this is negative
    #Jrate <- (fA/nint)*A0 + Jpredloss - (surv/nint)*J0 
    Jrate <- (fA)*A0 + Jpredloss - (surv)*J0
    A1 <- A0 + (Arate*dt)   # Update A numerically
    F1 <- F0 + (Frate*dt) + (sigmaF*NoiseF*dtZ)
    J1 <- J0 + (Jrate*dt)
    A1 <- max(A1,0.1)  # Force A1 greater than 0.1
    F1 <- max(F1,0.1)  # Force F greater than 0.1
    J1 <- max(J1,0.1)  # Force J greater than 0.1
    # Zooplankton dynamics
    Hrate <- DH*(Ho-H0) + alf*cPH*H0*P0 - cHF*H0*F0
    H1 <- H0 + (Hrate*dt) + (sigmaH*NoiseH*dtZ)
    H1 <- max(H1,0.1)  # Force H greater than 0.01
    # Phytoplankton dynamics
    Pbar <- P0   # Set P value for vertical integration
    gamvec <- GAMMA(Zvec,Pbar) # vertical sequence of light effect on growth
    gamI <- dZ*sum(gamvec)  # vertically integrated light effect on growth
    Prate <- (rP*Load*gamI*P0) - (mP*P0) - (cPH*H0*P0)  
    P1 <- P0 + (Prate*dt) + (sigmaP*NoiseP*dtZ)
    P1 <- max(P1,0.1)  # Force P greater than 0.1
    # reset initials for next iteration
    A0 = A1
    F0 = F1
    J0 = J1
    H0 = H1
    P0 = P1
  } # end iteration over nint
  
  # Construct list for output
  SimList <- list(A1,F1,J1,H1,P1)
  return(SimList)
  
}  # End iterative simulation function  ************************************

# Simulation by iteration with low sigma *************************************************
FWsim.iter.lowsig <- function(qE,A0,F0,J0,H0,P0,dt,dtZ,nint) { # Start simulation function
  # iterate nint times with a time step
  
  noise = 0.1*rnorm(3*nint)  # 0.1 decreases the noise
  noisemat = matrix(noise,nrow=nint,ncol=3)
  
  for (j in 1:nint) {
    nvec = noisemat[j,]
    NoiseF <- nvec[1]
    NoiseH <- nvec[2]
    NoiseP <- nvec[3]
    # Fish dynamics
    #Arate <- (surv/nint)*J0 - qE*A0 - ((1-surv)/nint)*A0
    Arate <- (surv)*J0 - qE*A0 - ((1-surv))*A0
    Frate <- DF*(Fo-F0) - cFA*F0*A0
    Jpredloss <- (-cJA*J0*A0)-(cJF*vuln*J0*F0/(hide + vuln + cJF*F0) ) # Note this is negative
    #Jrate <- (fA/nint)*A0 + Jpredloss - (surv/nint)*J0 
    Jrate <- (fA)*A0 + Jpredloss - (surv)*J0
    A1 <- A0 + (Arate*dt)   # Update A numerically
    F1 <- F0 + (Frate*dt) + (sigmaF*NoiseF*dtZ)
    J1 <- J0 + (Jrate*dt)
    A1 <- max(A1,0.1)  # Force A1 greater than 0.1
    F1 <- max(F1,0.1)  # Force F greater than 0.1
    J1 <- max(J1,0.1)  # Force J greater than 0.1
    # Zooplankton dynamics
    Hrate <- DH*(Ho-H0) + alf*cPH*H0*P0 - cHF*H0*F0
    H1 <- H0 + (Hrate*dt) + (sigmaH*NoiseH*dtZ)
    H1 <- max(H1,0.1)  # Force H greater than 0.01
    # Phytoplankton dynamics
    Pbar <- P0   # Set P value for vertical integration
    gamvec <- GAMMA(Zvec,Pbar) # vertical sequence of light effect on growth
    gamI <- dZ*sum(gamvec)  # vertically integrated light effect on growth
    Prate <- (rP*Load*gamI*P0) - (mP*P0) - (cPH*H0*P0)  
    P1 <- P0 + (Prate*dt) + (sigmaP*NoiseP*dtZ)
    P1 <- max(P1,0.1)  # Force P greater than 0.1
    # reset initials for next iteration
    A0 = A1
    F0 = F1
    J0 = J1
    H0 = H1
    P0 = P1
  } # end iteration over nint
  
  # Construct list for output
  SimList <- list(A1,F1,J1,H1,P1)
  return(SimList)
  
}  # End iterative low-sigma simulation function  ************************************

# Compute rolling window autocorrelation time
ACtime = function(x) {
  zzz=acf(x,lag=1,plot=FALSE)
  zz=zzz$acf[2]
  ACt=-1/log(zz)
  return(ACt)
}

# Compute rolling window autocorrelation
AClag1 = function(x) {
  zzz=acf(x,lag=1,plot=FALSE)
  zz=zzz$acf[2]
  #ACt=-1/log(zz)
  return(zz)
}

## END MODEL FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

