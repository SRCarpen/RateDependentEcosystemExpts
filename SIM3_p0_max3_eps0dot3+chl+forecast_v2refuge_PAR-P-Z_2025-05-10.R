# Fit DOC-Grazer-Chl model to 2024 data, December 2024
# Modified Fit2015_Zt-by-reg+PubFigs_2+3_drivers_2022-01-28:
#   * remove exp(-b*A) term
#   * add a refuge term f(K - A) as a low-A stabilizer
#   * use q=1, ordinary Michaelis-Menten grazing (see notes ~ how this affects eq)

# built on SIM2
# 2025_RateDependence_kNC_P_qE/SIM2_vinf-Crit_kNC_P0_2024_v2refuge_PAR-P-Z_2025-02-07.R
# with some features including forecast DLM from
# 2025_RateDependence_kNC_P_qE/SIM3_kNCspeeds+chl+forecast_V2refuge_PAR-P-Z_2025-04-12.R

rm(list = ls())
graphics.off()

library(moments)
library(runner)

source('Constants+Functions_PAR-P-Z_2025-03-29.R')

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

# fitted parameters from SIM0
#rmax =            0.456942350422245 
#hmax =          0.155761628853678

# test value of kNC for P0 effects
# Peter was Aquashaded in 2024 so use Paul Lake kNC0 for 'unmanipulated Peter'
kNC.test = quantile(kNC0[1:77],probs=0.5)  # kNC0[1:77] is Paul Lake data
# other average values for 2024
RZT = median(Ztvec)
RZbiom = quantile(ZBvec,probs=0.75)

# Equilibria versus P0 at given kNC =================================================
#p0grad = c(0.05,0.1,0.3,0.5,1,2,3,5)
np0 = 25
p0grad = seq(0.05,2,length.out=np0)
pneq = rep(0,np0)
peqmat = matrix(0,nr=np0,nc=3)
for(ip in 1:np0) {
  p0test = p0grad[ip]
  grout = grograze(kNC.test,p0test,RZT,RZbiom,rmax,hmax)
  eq = grout[[5]]
  #print(c(ip,'eq found: ',eq),quote=F)
  pneq[ip] = length(eq)
  if(pneq[ip] == 3) {peqmat[ip,] = eq}
}

# convert equilibria to chl
peqmat = peqmat/CChl
windows()
par(mfrow=c(1,1),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
yrange = range(peqmat,na.rm=T)
plot(p0grad,peqmat[,1],ylim=yrange,type='b',lwd=2,pch=19,col='blue',
     xlab='P load, mg/(m^2 d)',ylab='equilibria, mg Chl / m^3',
     main=paste('equilibria, kNC = ',round(kNC.test,3))) 
points(p0grad,peqmat[,2],type='b',lwd=2,pch=21,cex=1.2,col='red')
points(p0grad,peqmat[,3],type='b',lwd=2,pch=19,col='forestgreen')
grid(lwd=2,col='slategray')

# TEST -- FIND CRITICAL p0 =================================================

# find critical p0
cvec = as.vector(c(kNC.test,RZT,RZbiom,rmax,hmax)) # vector of constants
print('vector of constants',quote=F)
print(cvec)
guess = log(c(2,0.6))  # guesses of Chl & p0 critical value
critfit = optim(guess,crit.pload,gr=NULL,cvec,method='Nelder-Mead')
#critfit = optim(guess,crit.pload,gr=NULL,cvec,method='BFGS') # experiment with method
parfit = exp(critfit$par)
print(c('critical AChl and p0 ',parfit),quote=F)
critAchl = parfit[1]
p0.crit = parfit[2]
print(c('convergence',critfit$convergence),quote=F)

# END TEST TO FIND CRITICAL values ===============================================

# simulate and calculate EWS over a gradient of p0 ===================================

# simulate on Pload gradient
# noise and driver speed
sigma = fitsigma  # fitsigma is error from model fit to data; try also fitsigma/10
epsilon = 0.3  # speed of p0 relative to Chl; try values from 0.1 to 4
print('',quote=F)
print(c('***********************************************************'),quote=F)
print(c('sigma for chl simulation on p0 gradient',sigma),quote=F)
print(c('speed of p0 ',epsilon),quote=F)
print(c('***********************************************************'),quote=F)
print('',quote=F)

#  RUN CONTROL PARAMETERS
dt = 0.02  # days  0.04 in SIMs before 14 April 2025
dt.day = 1/dt
dtn = sqrt(dt)
nday = 120  # days to simulate  200 in SIMs before 14 April 2025
nstep = nday/dt   # time steps
Achl = rep(0,nstep)  # vector to hold results
#Achl[1] = median(Chl0[78:172])  # Peter Lake is 78 to 172
#Achl[1] = 10  # start near upper stable limb
Achl[1] = 0.1 # start near lower stable limb


# fixed values, in cases where we need them
kNCfix = median(kNC0[1:77])  # Paul Lake is 1 to 77
#kNCfix = median(kNC0[78:172])  # Peter Lake is 78 to 172
kNCmed = kNCfix # set the 'average' kNC for this run
# ecosystem state to simulate
#Chl,kNC,p0,ZT,Zb,rmax,hmax,sigma,dt,dtn
ZTmed = median(Zt0[78:172])  # Peter Lake is 78 to 172
ZBmed = median(ZB0[78:172])
#
# Compute p0 gradient using epsilon:
# this function simulates Ritchie-Sieber lamda given epsilon and nday
# Lamda ranges 0 to 2; we rescale lamda between 0.1 and 1.1
# We set lamda = p0sim
lamsim = function(epsilon,dt,nday,p0.crit) {
  nhalf = nday
  # integer days
  halfstep = c(0:nhalf)*dt
  tlam = c( -1*rev(halfstep),halfstep[2:nhalf] )
  nlam = length(tlam)
  #noise = rnorm(nstep)  # N(0,1) noise
  # calculate and plot lamda versus tstep
  lamdamax=2  # lamda ranges 0 to lamdamax
  lamda = (lamdamax/2)*( tanh(lamdamax*epsilon*tlam/2) +1)
  # plot lamda vs tstep
  windows(width=8,height=4)
  par(mfrow=c(1,2),mar=c(4.5,4.5,3,2)+0.1,cex.axis=1.6,cex.lab=1.6)
  plot(tlam,lamda,type='l',lwd=2,col='sienna',
       xlab = 'time step',ylab='lamda',main='Ritchie-Sieber simulation of p0')
  abline(h=p0.crit,lty=1,lwd=1,col='darkblue')
  grid()
  # rescale lamda from 0 - 2 to range between 0.1 and 1, for example
  # if p0 = a + b*lamda then 
  # 0.1 = a + b*0 = a
  # 1.1 = 0.1 + b*2 ==> b = 1/2
  # therefore p0sim = 0.1 + 0.5*lamda
  # or if Phi = 3
  # 3 = 0.1 + b*2  ==> b = 1.45
  Plo = 0.1  # lowest P load; 0.1 in 12 April expt
  Phi =  3 # 1.1    # highest P load; 1.1 in 12 April expt
  # solve for ahat and bhat
  bhat = (Phi - Plo)/(1.6 - 0.4)  # for epsilon 0.3 lamda range is 0.4 to 1.6
  ahat = Plo - 0.4*bhat  # 0.4 is low lamda for epsilon 0.3
  p0lam = ahat + bhat*lamda
  # rescale time axis to match the chlorophyll simulation
  lamfil = spline(x=tlam,y=p0lam,method='fmm',n=(1/dt)*(length(halfstep)-1))
  nlam = length(lamfil$x)
  print(c('n lamda simulated',nlam),quote=F)
  lamstep = c(1:nlam)
  lamda = lamfil$y
  lamday = lamstep*dt
  plot(lamday,lamfil$y,type='l',lwd=1,col='brown',
       xlab='time rescaled to days',ylab='lamda rescaled to p0',
       main='Ritchie-Sieber rescaled')
  abline(h=p0.crit,lty=1,lwd=1,col='darkblue')
  grid()
  # find time where p0 = p0.crit
  delt = lamda - p0.crit
  sdelt = sign(delt)
  dsdelt = c(0,-diff(sdelt))
  scross = lamstep[which(!dsdelt == 0)]
  dcross = lamday[which(!dsdelt == 0)]
  print(c('time step and day when lamda = p0.crit',scross,dcross),quote=F)
  outlist = list(lamstep,lamday,lamda,scross,dcross)
  return(outlist)
}

# run p0 simulation
# outlist = list(lamstep,lamday,lamda,scross,dcross)
print('',quote=F)
print('Report of p0 simulation',quote=F)
print(c('epsilon = ',epsilon),quote=F)
p0vals = lamsim(epsilon,dt,nday,p0.crit)
p0step = p0vals[[1]]
p0day = p0vals[[2]]
p0sim = p0vals[[3]]
tpoint = p0vals[[4]]
dpoint = p0vals[[5]]
print('End of p0 simulation',quote=F)
print('',quote=F)
#
# simulate
tstep = c(1:nstep)*dt
noise = rnorm(nstep)  # N(0,1) noise
for(i in 2:nstep) {
  # dAdt = function(Chl,kNC,p0,ZT,Zb,rmax,hmax,shock,sigma,dt,dtn)
  Achl[i] = dAdt(Achl[i-1],kNCfix,p0sim[i-1],ZTmed,ZBmed,rmax,hmax,noise[i-1],sigma,dt,dtn)
}

winlen = 25
v.Achl = runner(Achl,k=winlen,f=var)

windows(height=9,width=9)
par(mfrow=c(3,2),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.8)
plot(tstep,p0sim,type='l',lwd=1,col='sienna')
abline(h=p0.crit,lty=3,lwd=3,col='orangered')
grid()
plot(tstep,Achl,type='l',lwd=1,col='forestgreen',ylab = 'Chlorophyll')
grid()
plot(tstep[winlen:nstep],log10(v.Achl[winlen:nstep]),type='l',lwd=1,col='red',
     xlab='tstep',ylab='log10(var(Chl))',
     main=paste('window = ', winlen))
grid()

# run 1-dimensional DLM on Achl
# DLM <- function(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat) 
Acen = Achl - mean(Achl)  # center the data
npar = 1   # if matrix is centered we don't need the intercept
nX = length(Acen)
unit = rep(1,nX)

#
# Quick parameter estimate using LS 
print('Multivariate least squares regression to check the data',quote=F)
Xmat1 = as.matrix(Acen[2:nX])
Xmat0 = as.matrix(Acen[1:(nX-1)])
Xinv = solve(t(Xmat0)%*%Xmat0)
bmat = Xinv%*%t(Xmat0)%*%Xmat1
print('',quote=F)
#print('b matrix',quote=F)
#print(bmat)
#print('eigenvalues of b matrix',quote=F)
#beig = eigen(bmat,only.values=T)
#print(beig$values)
#print('prediction error variances',quote=F)
yhat = Xmat0%*%bmat
err = Xmat1 - yhat
verr = apply(err,2,var)
parcov = Xinv*verr
print('----------------------------------------------',quote=F)
print('LS regression result',quote=F)
print('bmat, verr, parcov',quote=F)
print(bmat)
print(verr)
print(parcov)
print('----------------------------------------------',quote=F)
#
# set up DLM
# delta:  # first round of experiments used 0.99; 0.996 is 1-day memory for this dt;
# 0.9998 is memory as long as the data, close to a ML estimate
delta = 0.996  
n.gamma = 1
d.gamma = verr
kh.day = 7 # forecast horizon, days
kh = dt.day*kh.day
print('***********************************************************')
print(c('DLM delta',delta),quote=F)
print(c('forecast horizon ',kh),quote=F)
print('***********************************************************')

dlm.est = DLMfore(delta,n.gamma,d.gamma,bmat,parcov,Xmat1,Xmat0,kh)
# DLM output list <- list(predix,varpredix,pars,parvar,Svec,fcast)
yhat = dlm.est[[1]]
modvar = dlm.est[[2]]  # process or 'model' variance
ar1 = dlm.est[[3]]  # parameter
sd1 = dlm.est[[4]]  # sd; these are small, ~ 1/10^8
vupdate = dlm.est[[5]]
fcast = dlm.est[[6]]

# add ar1 values to plot
plot(tstep[2:nX],ar1,type='l',lwd=1,col='blue',xlab='tstep',ylab='AR(1) coef')
grid()
abline(h=1)
# add prediction variance to plot
plot(tstep[2:nX],log10(modvar),type='l',lwd=1,col='tomato',xlab='tstep',ylab='log10(process var.)')
grid()
# add stationary variance to plot if ar1 is within the unit circle
# if ar1 is outside unit circle then substitute the model variance
vstat = ifelse(abs(ar1^2) <= 1,modvar/(1 - ar1^2),NA)  # stationary variance
plot(tstep[2:nX],log10(vstat),type='l',lwd=1,col='black',xlab='tstep',ylab='log10(stationary var.)',
     main='NA when ar1^2 > 1')
grid()

# reorganize key DLM results ===============================================================

# save the simulation ----------------------------------------------------------------------------------
#save(dt,dt.day,sigma,nstep,ikNC,speedvec,p0fix,tstep,Achl,Acen,nX,tpoint,
#     file='sim_for_DLM_obssig_eps=0.15.Rdata')
# ------------------------------------------------------------------------------------------------------

windows(height=12,width=12)
par(mfrow=c(3,2),mar=c(4,4.5,1,2)+0.1,cex.axis=1.8,cex.lab=1.8)

plot(tstep,p0sim,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='p0 Simulation')
abline(h=p0.crit,lty=3,lwd=3,col='orangered')
text(x=20,y=0.95*max(p0sim),paste('epsilon = ',epsilon),cex=1.8)
grid()

plot(tstep,Acen,type='l',lwd=4,col='skyblue',xlab='Day',
     ylab='Chlorophyll, centered')
points(tstep[2:nX],yhat,type='l',lwd=1,col='black')
grid()
abline(v=dpoint,lty=2,lwd=2,col='red')
text(x=100,y=0.65,paste('sigma = ',round(sigma,1)),cex=1.8)
legend('topleft',legend=c('simulation','DLM'),col=c('skyblue','black'),
       lwd=c(3,1),cex=1.8,bty='n')

#windows(height=6,width=6)
#par(mfrow=c(2,1),mar=c(2,4.5,1,2)+0.1,cex.axis=1.5,cex.lab=1.5)

plot(tstep[2:nX],ar1,type='l',lwd=1,col='darkblue',xlab='Day',
     ylab='AR(1) coef')
grid()
abline(h=1,lty=3,lwd=3,col='black')
abline(v=dpoint,lty=2,lwd=2,col='red')

plot(tstep[2:nX],sd1,log='y',type='l',lwd=1,col='blue',xlab='Day',
     ylab='s.e. of AR(1) coef)')
grid()
abline(v=dpoint,lty=2,lwd=2,col='red')

plot(tstep[2:nX],modvar,log='y',type='l',lwd=1,col='magenta',xlab='Day',
     ylab='Process Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='red')

# add kh-ahead forecast
# columns of fcast: kh-ahead forecast, obs, err, Q
fest = fcast[(kh+1):(nX-2),1]
fobs = fcast[(kh+1):(nX-2),2]
tfore = tstep[(401+kh):(nX-2)]
Qmat = fcast[(401+kh):(nX-2),4]
# Plot forecast vs obs
#yrange = range(c(fest,fobs),na.rm=T)
#plot(fest,fobs,type='p',pch=20,cex=0.5,col='black',xlab='mean of 250-step forecast',
#     ylab='Observation')
#grid()
plot(tfore,Qmat,log='y',type='l',lwd=1,col='purple',xlab='Day',
     ylab='Forecast Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='red')
text(x=100,y=0.99*max(Qmat),paste('days ahead = ',kh.day),cex=1.8)
#text(x=100,y=0.85*max(Qmat),paste(c(kh/dt.day,' days')),cex=1.8)

# save the simulation ----------------------------------------------------------------------------------
save(epsilon,sigma,dpoint,Pload0,tstep,p0sim,p0.crit,nX,Achl,kh,kh.day,Acen,yhat,
     ar1,sd1,modvar,tfore,Qmat,file='p0sim_EWS_eps0dot3.Rdata')
# ------------------------------------------------------------------------------------------------------

