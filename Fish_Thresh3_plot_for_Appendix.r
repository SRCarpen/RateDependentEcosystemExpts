# Simulation of the full food web for investigating the squeal
# This program is like Fish_Thresh3.r but plots the A rate curve vs A
# SRC 23 June 2007

# To run, type source("Fish_Thresh3_plot.R")
# To clear all variables, type rm(list = ls())
rm(list = ls())

# This line closes all open graphics windows
graphics.off()

#
#  PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# Zooplankton ***********************
Ho <- 1    # Refuge biomass
DH <- 0.5  # Diffusion parameter
cHF <- 0.1  # Consumption rate by planktivore
alf <- 0.3  # Conversion efficiency of consumed phytoplankton to zooplankton
cPH <- 0.25  # Consumption rate of phytoplankton by zooplankton

# Fish *****************************
qELO <- 1  # First Catchability x Effort
qEHI <- 4  # Second Catchability x Effort
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

#
# INITIAL CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~
#
#Ainit <- 120  # 20 for simulations in May 07
#Finit <- Fo
Hinit <- Ho
Pinit <- 6

#
# FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#Phytoplankton Growth 
GAMMA <- function(z,Pbar) {
Iz <- I0*exp(-z*(eps0+epsP*Pbar))
rate <- (1/Fmax)*(1 - exp(-k_sat*Iz))*exp(-k_inh*Iz)
}

# Simulation -- Deterministic model
FWdet <- function(qE,Ainit,Finit) { # Start simulation function

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
NoiseF <- rnorm(nint)*0 # Zero out the noise terms
NoiseH <- rnorm(nint)*0
NoiseP <- rnorm(nint)*0
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
F1 <- F0 + (Frate*dt) + (sigma*NoiseF[it]*dtZ)
J1 <- J0 + (Jrate*dt)
A1 <- max(A1,0.1)  # Force A1 greater than 0.1
F1 <- max(F1,0.1)  # Force F greater than 0.1
J1 <- max(J1,0.1)  # Force J greater than 0.1
# Zooplankton dynamics
Hrate <- DH*(Ho-H0) + alf*cPH*H0*P0 - cHF*H0*F0
H1 <- H0 + (Hrate*dt) + (sigma*NoiseH[it]*dtZ)
H1 <- max(H1,0.01)  # Force H greater than 0.01
# Phytoplankton dynamics
Pbar <- P0   # Set P value for vertical integration
gamvec <- GAMMA(Zvec,Pbar) # vertical sequence of light effect on growth
gamI <- dZ*sum(gamvec)  # vertically integrated light effect on growth
Prate <- (rP*Load*gamI*P0) - (mP*P0) - (cPH*H0*P0) 
P1 <- P0 + (Prate*dt) + (sigma*NoiseP[it]*dtZ) 
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
AFrate <- c(A0-Ainit,F0-Finit)

}  # End simulation function

# Loss function for minimization
Alossfct <- function(parvec) {
# Basic idea: if f := A(t+1) - A(t) then the critical point
# occurs where f=0 and df/dA = 0
qE <- parvec[1]
Ainit <- parvec[2]
Finit <- parvec[3]
funAF <- FWdet(qE,Ainit,Finit)
epsilon <- 1
Aplus <- Ainit+epsilon
funAFplus <- FWdet(qE,Aplus,Finit)
dfundA <- (funAF[1] - funAFplus[1])/epsilon
Aloss <- sum(funAF*funAF) + (dfundA*dfundA)
}

# Time series plots
TimePlot <- function(year,AB,JB,FB,HY,PY) { # Start time series plots
# Plot results
windows()
par(mfrow=c(4,1))  # window with 4 rows, 1 column
Fuplim <- max(c(AB,JB))
Flimits <- c(0,Fuplim)
plot(year,JB,xlab=" ",ylab=" ",ylim=Flimits,type='l',col='purple',lwd=2)
title(main="Adult Blue, Juvenile Purple",xlab="year",ylab="Piscivore")
points(year,AB,type='l',col='blue',lwd=2)
grid()
plot(year,FB,xlab="year",ylab="Planktivore",type='l',col='red',lwd=2)
grid()
plot(year,HY,xlab="year",ylab="Herbivore",type='l',col='blueviolet',lwd=2)
grid()
plot(year,PY,xlab="year",ylab="Phytoplankton",type='l',col='darkgreen',lwd=2)
title(main="Plankton",xlab="year",ylab="Phytoplankton")
grid()
} # End time series plots

#
# RUN SIMULATION  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Setting up for the plots

# First make the rate plots

#nYbase <- 100  # Number of years for burn-in to initial equilibrium
#nYexpt <- 1  # Number of years for statistical study
nY <- 1 #nY <- nYbase + nYexpt
year <- seq(1,nY, by=1)
nint <- 32  # Time steps per year 
dt <- 1/nint
dtZ <- sqrt(dt)

# Load critical qE
#save(CritPt,Feq,file='Crit_qE_AdultPisc.Rdata')
load(file='Crit_qE_AdultPisc.Rdata')
qEcrit = CritPt$par[1]
Acrit = CritPt$par[2]

# Rate plots
qEvec <- c(1.0,1.2,1.4,1.6,qEcrit,1.8,1.9)  # shorter list                         
nQ <- length(qEvec)                                                                                                                               
Amin <- 0.01                                                                    
Amax <- 250                                                                     
nA <- 50                                                                       
Agrad <- seq(Amin,Amax,length.out=nA)
Finit <- 1.7  # 0.5*Fo # value used in supplement V3                                           
                                                                                
# Vectors for rate output                                                       
Arate0 <- rep(0, times=nA*nQ)                                                   
Amat0 <- Arate0                                                                 
count <- 0                                                                      
for(iQ in 1:nQ) {  # Start loop over qE                                         
qE <- qEvec[iQ]                                                                 
for (iA in 1:nA) { # Start loop over nA                                         
count <- count+1                                                                
Aval <- Agrad[iA] 
RateVec <- FWdet(qE,Aval,Finit)                                                               
Arate0[count] <- RateVec[1]                                                     
Amat0[count] <- Aval                                                            
} # End loop over A                                                             
} # End loop over qE                                                            
                                                                                
Arates <- matrix(Arate0, nrow=nA, ncol=nQ)                                      
Amat <- matrix(Amat0, nrow=nA, ncol=nQ)                                         
qElabel <- round(qEvec,digits=4)           
colvec = rainbow(nQ,start=0,end=5/6)

# Plot rate curves versus A for each Q                                          
windows()  # Overlay plots    
par(mar=c(4,4.5,2,1)+0.1,cex.axis=1.4,cex.lab=1.6,cex.main=1.5,font.axis=2,
    font.lab=2,font.main=2)
matplot(Amat,Arates,xlab=' ',ylab=' ',type='l',lty=1,lwd=2,                     
 col=colvec)                                               
 title(main=paste(c("Harvest Rates, red to purple"), sep = "=", 
 collapse= ",  "), 
 xlab='Piscivore density',ylab='Net Population Growth Rate',cex.lab=1.5)                                      
grid()              
points(Acrit,0,type='p',pch=21,cex=3,lwd=2.5,col='black')
legend('bottomleft',legend=qElabel,col=colvec,lwd=rep(2,nQ),
       bg='white',bty='o',box.col='black',box.lwd=1.5,cex=1.4,
       title='Harvest Rate')
                                                                                