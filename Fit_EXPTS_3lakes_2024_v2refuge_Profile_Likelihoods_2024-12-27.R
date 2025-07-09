# Fit DOC-Grazer-Chl model to 2024 data and find profile likelihoods, December 2024
# Modified Fit2015_Zt-by-reg+PubFigs_2+3_drivers_2022-01-28:
#   * remove exp(-b*A) term
#   * add a refuge term f(K - A) as a low-A stabilizer
#   * use q=1, ordinary Michaelis-Menten grazing (see notes ~ how this affects eq)

rm(list = ls())
graphics.off()

# CONSTANTS

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

# FUNCTIONS 

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

# END FUNCTIONS =======================================================

# Input data compiled by Inputs_3lakes_PAR-P-Phyto_2024_2024-12-23.R
# save vectors for fitting the 3 lakes together
# each vector has, in order, Paul, Peter, Tuesday
# DOY0 has first digit 1 for Paul, 2 for Peter, 3 for Tuesday
# fetch is medsquarea of each lake
#ndoy = length(DOY0)
#save(ndoy,Chl0,Chl1,Zt0,ZB0,DOY0,kNC0,kPAR0,Pload0,fetch,
#     file='chl+drivers_3lakes.Rdata')

load(file='chl+drivers_3lakes.Rdata')

# write variates as time 0 and time 1 and drivers as time 0
#Chl0 = Chl0  # NOTE SOME NAMES ARE UNCHANGED
#Chl1 = Chl1
Ztvec = Zt0
#DOCvec = DOC[1:(ndoy-1)]  # not needed; use non-chl PAR ext coef instead
ZBvec = ZB0
#DOY0 = DOY0
#kNC0 = kNC0
#kPAR0 = kPAR0
P0vec = Pload0

windows(height=10,width=5)
par(mfrow=c(5,1),mar=c(2.5,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(DOY0,P0vec,type='l',lwd=2,col='limegreen')
abline(v=2,lty=3,lwd=2,col='black')
abline(v=3,lty=3,lwd=2,col='black')
plot(DOY0,Chl0,type='l',lwd=2,col='forestgreen')
abline(v=2,lty=3,lwd=2,col='black')
abline(v=3,lty=3,lwd=2,col='black')
plot(DOY0,kNC0,type='l',lwd=2,col='sienna')
abline(v=2,lty=3,lwd=2,col='black')
abline(v=3,lty=3,lwd=2,col='black')
plot(DOY0,ZB0,type='l',lwd=2,col='red')
abline(v=2,lty=3,lwd=2,col='black')
abline(v=3,lty=3,lwd=2,col='black')
plot(DOY0,Zt0,type='l',lwd=2,col='blue')
abline(v=2,lty=3,lwd=2,col='black')
abline(v=3,lty=3,lwd=2,col='black')

# component extinction coefs entered under 'constants'
#eps0 <- 0.0213  # Light extinction by pure water, m-1
#epsDOC = 0.0514  # DOC light extinction coef, m2 g-1
#epsP <- 0.0177   # Phytoplankton, m2 (mg chl)-1

# Retrieve multi-lake regression model to predict Zt 
#  by Zmix=f(DOC_etc)_2021-12-31.R
#load(file='ZmixModel.Rdata')  # 7 lakes
#load(file='ZmixModel_w_Hbird.Rdata')  # 8 lakes
load(file='ZmixModel_LakeYears_Hbird.Rdata')  # lake years
print(summary(reg1)) 
pars = reg1$coefficients
#
# Use median squarea sqrt(lake area in ha) values in fetch

# guess parameters rmax, hmax, sigma
sdhat = sd( (Chl1-Chl0),na.rm=T)
#guess = c(2,0.2,sdhat)
guess = log( c(1,1,sdhat) )  # parameters are exponentiated in the model functions
# fit model
fit1 = optim(guess,OneStep,gr=NULL,P0vec,Chl0,Chl1,Ztvec,kNC0,ZBvec,(ndoy),
             method='Nelder-Mead')
#method='SANN')
print('Optim results',quote=F)
#print(c('parameters',fit1$par),quote=F)
print(c('parameters',exp(fit1$par)),quote=F)
print(c('final NLL',fit1$value),quote=F)
print(c('iterations',fit1$counts),quote=F)
print(c('convergence, 0 is good',fit1$convergence),quote=F)

Ahat = Yhat(fit1$par,P0vec,Chl0,Chl1,Ztvec,kNC0,ZBvec,(ndoy))

# name the fitted parameters
rmax = exp(fit1$par[1])
hmax = exp(fit1$par[2])
print(c('rmax = ',rmax,', hmax = ',hmax),quote=F)

windows()
par(mfrow=c(1,1),mar=c(4,4.4,2,2)+0.1,cex.axis=1.8,cex.lab=1.8)
plot(DOY0,Ahat,type='l',lwd=2,col='forestgreen',xlab='DOY',ylab='Chl')
points(DOY0,Chl1,type='p',pch=20,col='lightseagreen')

# Profile likelihood for hmax
rpar1 = seq(-1,0,length.out = 11)
rpar2 = -1*rev(rpar1)
logrepar = c(rpar1[1:10],rpar2)
repar=10^logrepar
# parameter gradient
hgrad = repar*hmax
nh = length(hgrad)
NLLprof = rep(0,nh)
parvec = fit1$par
for(i in 1:nh) {
  parvec[2] = log(hgrad[i])
  NLLprof[i] = OneStep(parvec,P0vec,Chl0,Chl1,Ztvec,kNC0,ZBvec,(ndoy))
}

windows()
par(mfrow=c(1,1),mar=c(4,4.4,2,2)+0.1,cex.axis=1.8,cex.lab=1.8)
plot(hgrad,NLLprof,type='b',lwd=1,pch=20,cex=1.5,col='black',
     xlab='hmax gradient',ylab='Profile NLL')
abline(v=hmax,lty=2,lwd=1,col='red')

windows()
par(mfrow=c(1,1),mar=c(4,4.4,2,2)+0.1,cex.axis=1.8,cex.lab=1.8)
plot(repar,NLLprof,type='b',lwd=1,pch=20,cex=1.5,col='black',
     xlab='multiple of nominal hmax',ylab='Profile NLL')
abline(v=1,lty=2,lwd=1,col='red')
