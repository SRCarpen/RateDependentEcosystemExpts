# Fit DOC-Grazer-Chl model to 2024 data, December 2024
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
abline(v=2,lty=3,lwd=2,col='black')
abline(v=3,lty=3,lwd=2,col='black')

# Test growth and grazing curves

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

# Plot kNC0 versus DOC equivalents
windows()
par(mfrow=c(1,1),mar=c(4.5,4.5,2,1)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(kNC0,(kNC0/epsDOC),type='p',pch=20,col='sienna',
     xlab='Non-Chl PAR Extinction',
     ylab='DOC with same PAR Extinction, g/m^3',
     main='DOC equivalent to Non-Chl PAR Extinction')
grid(lwd=2,col='slategray')

# Calculate growth curves and eq for Peter Lake
p0 = 1 # choose P load
# use mean drivers
kNC.test = quantile(kNC0,probs=0.1)
DOC.test = kNC.test/epsDOC
RZT = median(Ztvec)
RZbiom = quantile(ZBvec,probs=0.75)
# Increase hmax for a steeper initial slope =====================================
hmaxML = hmax # save ML fit hmax
hmax = 10*hmaxML   # factor should be roughly proportional to p0 value
# ===============================================================================
Rout = grograze(kNC.test,p0,RZT,RZbiom,rmax,hmax)
# grograze output:
RAgrad = Rout[[1]]
Rgro = Rout[[2]]
Rgrz = Rout[[3]]
Rnet = Rout[[4]]
Req.C = Rout[[5]]
Req.chl = Req.C/CChl

print(c('equilibria as C = ',Req.C),quote=F)
print(c('equilibria as Chl = ',Req.chl),quote=F)

# Plots for Peter Lake
windows(width=5,height=10)
par(mfrow=c(2,1),mar=c(2.5,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
yrange = range(Rgro,Rgrz)
plot(RAgrad,Rgrz,type='l',lwd=2,col='blue',ylim=yrange,
     xlab='Phytoplankton (C)',ylab='Growth & Grazing',
     main=paste('kNC = ',round(kNC.test,3),'(DOC ~ ',round(DOC.test,1),'), P load = ',p0) )
points(RAgrad,Rgro,type='l',lwd=2,col='forestgreen')
abline(h=0,lty=3,lwd=2)
#
par(mar=c(4,4.5,1,2)+0.1)
plot(RAgrad,Rnet,type='l',lwd=2,col='blue',
     xlab='Phytoplankton (C)',ylab='Growth - Grazing') 
#main=paste('DOC = ',round(DOCf,1),', P load = ',peffect*5) )
abline(h=0,lty=3,lwd=2)
#points(Req.C,rep(0,length(Req.C)),type='p',pch=c(19,21,19),lwd=2,cex=2.5,col='black')

# Equilibria versus kNC over observed gradient
q.kNC = unname(quantile(kNC0,probs = c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)) )
klo = q.kNC[1]
khi = max(kNC0) #q.kNC[9]
nk = 20
kvec = seq(klo,khi,length.out=nk)
neq = rep(0,nk)
eqmat = matrix(0,nr=nk,nc=3)
p0=1
for(ik in 1:nk) {
  ktest = kvec[ik]
  grout = grograze(ktest,p0,RZT,RZbiom,rmax,hmax)
  eq = grout[[5]]
  #print(c(ik,'eq found: ',eq),quote=F)
  neq[ik] = length(eq)
  if(neq[ik] == 3) {eqmat[ik,] = eq}
}

windows()
par(mfrow=c(1,1),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(kvec,neq,type='p',pch=20,col='blue',
     xlab='Non-Chl PAR extinction',ylab='number of equilibria') 

# convert equilibria to chl
eqmat = eqmat/CChl
windows()
par(mfrow=c(1,1),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
yrange = range(eqmat,na.rm=T)
plot(kvec,eqmat[,1],ylim=yrange,type='b',lwd=2,pch=19,col='blue',
     xlab='Non-Chl PAR extinction',ylab='equilibria, mg Chl / m^3',
     main=paste('equilibria, P load = ',p0)) 
points(kvec,eqmat[,2],type='b',lwd=2,pch=21,cex=1.2,col='red')
points(kvec,eqmat[,3],type='b',lwd=2,pch=19,col='forestgreen')
grid(lwd=2,col='slategray')