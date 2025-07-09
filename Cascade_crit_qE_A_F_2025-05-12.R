# Program to find critical qE (harvest coef), A (adult bass), F (planktivores) 
#  for Cascade model
# Based on FishThresh2.R
# SRC 2025-05-12

rm(list = ls())
graphics.off()

# constants ====================================================================

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

# Zooplankton ***********************
Ho <- 1    # Refuge biomass
DH <- 0.5  # Diffusion parameter
cHF <- 0.1  # Consumption rate by planktivore
alf <- 0.3  # Conversion efficiency of consumed phytoplankton to zooplankton
cPH <- 0.25  # Consumption rate of phytoplankton by zooplankton


Fo <- 100  # Refuge density of planktivores
DF <- 0.1  # Diffusion parameter for planktivores
Finit <- Fo
Hinit <- 1
Pinit <- 6

fA <- 2  # Fecundity of adult piscivore
cJA <- 0.001  # Density dependent mortality rate of juveniles
cJF <- 0.5  # Consumption of juveniles by planktivores
cFA <- 0.3  # Consumption of planktivores by adult piscivores
vuln <- 1  # Vulnerability coefficient
hide <- 8  # Hiding coefficient
surv <- 0.5  # Overwinter survivorship of adults
Fo <- 100  # Refuge density of planktivores
DF <- 0.1  # Diffusion parameter for planktivores
sigma <- 0.1  # SD of additive noise for planktivores (0.1 in May 07)
A2biom <- 0.2  # Convert A to kg / ha
J2biom <- 0.05  # Convert J to kg / ha
F2biom <- 1  # Convert F to kg / ha

# functions =======================================================================

#Phytoplankton Growth 
GAMMA <- function(z,Pbar) {
  Iz <- I0*exp(-z*(eps0+epsP*Pbar))
  rate <- (1/Fmax)*(1 - exp(-k_sat*Iz))*exp(-k_inh*Iz)
}

# Loss function for qE and A
Alossfct <- function(parvec) {
  # Basic idea: if f := A(t+1) - A(t) then the critical point
  # occurs where f=0 and df/dA = 0
  qE <- parvec[1]
  Ainit <- parvec[2]
  funA <- FWdet(qE,Ainit)
  Aplus <- Ainit+0.1
  funAplus <- FWdet(qE,Aplus)
  dfundA <- (funA - funAplus)/0.1
  Aloss <- (funA*funA) + (dfundA*dfundA)
}

# Simulation -- Deterministic model
FWdet <- function(qE,Ainit) { # Start simulation function
  
  A0 <- Ainit  # Set up initial conditions
  F0 <- Finit
  J0 <- fA*Ainit
  H0 <- Hinit
  P0 <- Pinit
  
  AY <- Ainit      # Vectors to hold "annual" statistics
  FY <- Finit
  JY <- fA*Ainit
  HY <- Hinit
  PY <- Pinit
  
  Ait <- rep(0,nint*nY)    # Vectors to hold detailed time series
  Fit <- rep(0,nint*nY)
  Jit <- rep(0,nint*nY) 
  Hit <- rep(0,nint*nY) 
  Pit <- rep(0,nint*nY) 
  
  # Compute the simulation
  ico <- 0
  for(iY in 1:nY)  {   # start loop over years
    #NoiseF <- rnorm(nint)*0 # Zero out the noise terms
    #NoiseH <- rnorm(nint)*0
    #NoiseP <- rnorm(nint)*0
    J0 <- fA*A0  # Juveniles at start of year
    for(it in 1:nint) {  # start integration within a year
      ico <- ico+1
      # Store intermediate results
      Ait[ico] <- A0
      Fit[ico] <- F0
      Jit[ico] <- J0
      Hit[ico] <- H0
      Pit[ico] <- P0
      # Fish dynamics
      Arate <- -qE*A0     # Update A numerically
      Frate <- DF*(Fo-F0) - cFA*F0*A0
      Jrate <- (-cJA*J0*A0)-(cJF*vuln*J0*F0/(hide + vuln + cJF*F0) )
      A1 <- A0 + (Arate*dt)   # Update A numerically
      F1 <- F0 + (Frate*dt) #+ (sigma*NoiseF[it]*dtZ)
      J1 <- J0 + (Jrate*dt)
      A1 <- max(A1,0.1)  # Force A1 greater than 0.1
      F1 <- max(F1,0.1)  # Force F greater than 0.1
      J1 <- max(J1,0.1)  # Force J greater than 0.1
      # Zooplankton dynamics
      Hrate <- DH*(Ho-H0) + alf*cPH*H0*P0 - cHF*H0*F0
      H1 <- H0 + (Hrate*dt) #+ (sigma*NoiseH[it]*dtZ)
      H1 <- max(H1,0.01)  # Force H greater than 0.01
      # Phytoplankton dynamics
      Pbar <- P0   # Set P value for vertical integration
      gamvec <- GAMMA(Zvec,Pbar) # vertical sequence of light effect on growth
      gamI <- dZ*sum(gamvec)  # vertically integrated light effect on growth
      Prate <- (rP*Load*gamI*P0) - (mP*P0) - (cPH*H0*P0) 
      P1 <- P0 + (Prate*dt) #+ (sigma*NoiseP[it]*dtZ) 
      P1 <- max(P1,0.1)  # Force P greater than 0.1
      # Update state variables
      A0 <- A1
      F0 <- F1
      J0 <- J1
      H0 <- H1
      P0 <- P1
    }  # End loop for integration within a year
    
    # Save annual end-of-year statistics
    #Ystart <- ico-nint+1
    #AY[iY] <- mean(Ait[Ystart:ico])
    #FY[iY] <- mean(Fit[Ystart:ico])
    #JY[iY] <- mean(Jit[Ystart:ico])
    #HY[iY] <- mean(Hit[Ystart:ico])
    #PY[iY] <- mean(Pit[Ystart:ico])
    
    # Update A0 for overwinter survivorship to next year
    A0 <- surv*(J1+A1)
    # Zero out the J0's
    J0 <- 0
    
  }  # End loop over years
  
  # Sample from stationary distribution
  #Sstart <- (nY-nYexpt)*nint + 1
  #Fst <- Fit[Sstart:ico]
  #Jst <- Jit[Sstart:ico]
  #Hst <- Hit[Sstart:ico]
  #Pst <- Pit[Sstart:ico]
  
  # Construct list for output
  #SimList <- list(year,AY,FY,JY,HY,PY,Fst,Jst,Hst,Pst)
  
  # Return rate
  Arate <- A0-Ainit
  return(Arate)
}  # End simulation function

# Try to find critical qE and A values using Alossfct <- function(parvec) ============
# run parameters
print('Bounded L-BFGS-B search for optimal qE and Adult piscivore',quote=F)
nint <- 30  # Time steps per year 
print(c('time steps / year = ',nint),quote=F)
dt <- 1/nint
dtZ = sqrt(dt)
nY=100
#
qEguess <- 1.8
Aguess <- 100   
parvec <- c(qEguess,Aguess)
lobound <- c(1,2)
upbound <- c(5,300)
CritPt <- optim(parvec,Alossfct,method="L-BFGS-B",lower=lobound,upper=upbound)
print(CritPt)

# solve for F at critical A
# Frate <- DF*(Fo-F0) - cFA*F0*A0
# CritPt[2] = A0 at critical point which is also an equilibrium
# set Frate to 0 and solve for F:
# F at eq = -1*DF/(DF - cFA*CritPt$par[2])
Feq = -1*DF/(DF - cFA*CritPt$par[2])
print(c('Planktivore at eq = ',Feq),quote=F)

save(CritPt,Feq,file='Crit_qE_AdultPisc.Rdata')

ntest = 300
Agrad = seq(1,300,length.out=ntest)
qEcrit = CritPt$par[1]
rateq = rep(1,ntest)
for(i in 1:ntest) {
  rateq[i] = FWdet(qEcrit,Agrad[i])
}

windows()
plot(Agrad,rateq,type='l',lwd=2,col='blue')
grid()
abline(h=0,lwd=2)


