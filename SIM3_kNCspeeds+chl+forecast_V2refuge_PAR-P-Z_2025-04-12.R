#SIM3_speeds_forecast_v2refuge_PAR-P-Z_2025-04-12

# Built on SIM2
# This version collects control parameters, simulates kNC based on epsilon, 
# then simulates chl and calculates EWS including from forecast variance
# Original code:
#    /SIM2_vinf-testcrit+forecast_v2refuge_PAR-P-Z_2025-03-29.R
# Notes on original code:
# Modified Fit2015_Zt-by-reg+PubFigs_2+3_drivers_2022-01-28:
#   * remove exp(-b*A) term
#   * add a refuge term f(K - A) as a low-A stabilizer
#   * use q=1, ordinary Michaelis-Menten grazing (see notes ~ how this affects eq)

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

# fitted parameters from SIM0 under 'constants'
#rmax =            0.456942350422245 
#hmax =          0.155761628853678

# simulate on kNC gradient ==========================================================

# SET SIMULATION EPSILON AND SIGMA BETWEEN LINES 85-90

# LAKE CONSTANTS FOR THE RUN
p0 = 1 # choose P load
# use fixed drivers
# Peter was Aquashaded in 2024 so use Paul Lake kNC0 for 'unmanipulated Peter'
kNC.test = quantile(kNC0[1:77],probs=0.5)  # kNC0[1:77] is Paul Lake data
DOC.test = kNC.test/epsDOC
RZT = median(Ztvec)
RZbiom = quantile(ZBvec,probs=0.75)
# ecosystem state needed to run a simulation
#Chl,kNC,p0,ZT,Zb,rmax,hmax,sigma,dt,dtn
ZTmed = median(Zt0[78:172])  # Peter Lake is 78 to 172
ZBmed = median(ZB0[78:172])
# noise and driver speed
sigma = fitsigma  # fitsigma is error from model fit to data; try also fitsigma/10
epsilon = 1  # speed of kNC relative to Chl; try values from 0.1 to 4
print('',quote=F)
print(c('***********************************************************'),quote=F)
print(c('sigma for chl simulation on kNC gradient',sigma),quote=F)
print(c('speed of kNC ',epsilon),quote=F)
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
Achl[1] = 10  # start near upper stable limb
#Achl[1] = 2 # start near lower stable limb

# find critical kNC
cvec = as.vector(c(p0,RZT,RZbiom,rmax,hmax)) # vector of constants
print('vector of constants',quote=F)
print(cvec)
guess = c(5.6,1.507)  # guesses of Chl & kNC critical value
critfit = optim(guess,crit.kNC,gr=NULL,cvec,method='Nelder-Mead')
parfit = critfit$par
print('',quote=F)
print(c('critical AChl and kNC ',parfit),quote=F)
kNC.crit = parfit[2]
print(c('convergence',critfit$convergence),quote=F)

# simulate kNC over a gradient that includes the critical value
# this function simulates Ritchie-Sieber lamda given epsilon and nday
# We set lamda = kNCsim
lamsim = function(epsilon,dt,nday,kNC.crit) {
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
       xlab = 'time step',ylab='lamda',main='Ritchie-Sieber simulation of kNC')
  abline(h=kNC.crit,lty=1,lwd=1,col='darkblue')
  grid()
  # rescale time axis to match the chlorophyll simulation
  lamfil = spline(x=tlam,y=lamda,method='fmm',n=(1/dt)*(length(halfstep)-1))
  nlam = length(lamfil$x)
  print(c('n lamda simulated',nlam),quote=F)
  lamstep = c(1:nlam)
  lamda = lamfil$y
  lamday = lamstep*dt
  plot(lamday,lamfil$y,type='l',lwd=1,col='brown',
       xlab='time rescaled to days',ylab='lamda',
       main='Ritchie-Sieber with time rescaled')
  abline(h=kNC.crit,lty=1,lwd=1,col='darkblue')
  grid()
  # find time where kNC = kNC.crit
  delt = lamda - kNC.crit
  sdelt = sign(delt)
  dsdelt = c(0,-diff(sdelt))
  scross = lamstep[which(!dsdelt == 0)]
  dcross = lamday[which(!dsdelt == 0)]
  print(c('time step and day when lamda = kNC.crit',scross,dcross),quote=F)
  outlist = list(lamstep,lamday,lamda,scross,dcross)
  return(outlist)
}

# run kNC simulation
# outlist = list(lamstep,lamday,lamda,scross,dcross)
print('',quote=F)
print('Report of kNC simulation',quote=F)
print(c('epsilon = ',epsilon),quote=F)
kNCvals = lamsim(epsilon,dt,nday,kNC.crit)
kNCstep = kNCvals[[1]]
kNCday = kNCvals[[2]]
kNCsim = kNCvals[[3]]
tpoint = kNCvals[[4]]
dpoint = kNCvals[[5]]
print('End of kNC simulation',quote=F)
print('',quote=F)
#
# simulate
tstep = c(1:nstep)*dt
noise = rnorm(nstep)  # N(0,1) noise
for(i in 2:nstep) {
  Achl[i] = dAdt(Achl[i-1],kNCsim[i-1],p0,ZTmed,ZBmed,rmax,hmax,noise[i-1],sigma,dt,dtn)
}

winlen = 25
v.Achl = runner(Achl,k=winlen,f=var)

windows(height=9,width=9)
par(mfrow=c(3,2),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.8)
plot(tstep,kNCsim,type='l',lwd=1,col='sienna')
abline(h=kNC.crit,lty=3,lwd=3,col='orangered')
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

windows(height=12,width=12)
par(mfrow=c(3,2),mar=c(4,4.5,1,2)+0.1,cex.axis=1.8,cex.lab=1.8)

plot(tstep,kNCsim,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='kNC Simulation')
abline(h=kNC.crit,lty=3,lwd=3,col='orangered')
text(x=20,y=1.8,paste('epsilon = ',epsilon),cex=1.8)
grid()

plot(tstep,Acen,type='l',lwd=4,col='skyblue',xlab='Day',
     ylab='Chlorophyll, centered')
points(tstep[2:nX],yhat,type='l',lwd=1,col='black')
grid()
abline(v=dpoint,lty=2,lwd=2,col='red')
text(x=100,y=0.65,paste('sigma = ',round(sigma,1)),cex=1.8)
legend('bottomleft',legend=c('simulation','DLM'),col=c('skyblue','black'),lwd=c(3,1),cex=1.8)

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
text(x=33,y=0.98*max(Qmat),paste('days ahead = ',kh.day),cex=1.8)
#text(x=100,y=0.85*max(Qmat),paste(c(kh/dt.day,' days')),cex=1.8)

# save the simulation ----------------------------------------------------------------------------------
save(epsilon,sigma,dpoint,Pload0,tstep,kNCsim,kNC.crit,nX,Achl,kh,kh.day,Acen,yhat,
     ar1,sd1,modvar,tfore,Qmat,file='kNCsim_EWS_eps1.Rdata')
# ------------------------------------------------------------------------------------------------------


# SUMMARY:  early warnings for synthesis table ============================================
print('',quote=F)
print('SUMMARY FOR SYNTHESIS TABLE',quote=F)
print(c('epsilon: ',epsilon,', crit time: ',dpoint,', sigma: ',sigma),quote=F)
print('',quote=F)
# find nearest tstep of ar1 crossing
ews1 = as.data.frame(cbind(tstep[2:nX],ar1))
colnames(ews1) = c('t','ar1')
ews1a = subset(ews1,subset=(t > 60 & t < 80))
ews1a$ar12 = (ews1a$ar1-1)^2
tews = which.min(ews1a$ar12)
ewar1 = ews1a$t[tews]
print('Differences from critical times: ',quote=F)
print('positive is early warning, negative is late warning',quote=F)
print(c('time of ar1 = 1 ',ewar1,', difference from crit time = ',dpoint-ewar1),quote=F)
#
# find nearest maxima of model variances
ews2 = as.data.frame(cbind(tstep[2:nX],sd1,modvar))
colnames(ews2) = c('t','sd1','modvar')
ews2a = subset(ews2,subset=(t > 62 & t < 80))
tsd1 = which.max(ews2a$sd1)
ewarsd1 = ews2a$t[tsd1]
print(c('time of max(ar1 sd) ',ewarsd1,', difference from crit time = ',dpoint-ewarsd1),quote=F)
#
tmv = which.max(ews2a$modvar)
ewarmv = ews2a$t[tmv]
print(c('time of max(modvar) ',ewarmv,', difference from crit time = ',dpoint-ewarmv),quote=F)
#
ews3 = as.data.frame(cbind(tfore,Qmat))
colnames(ews3) = c('t','Q')
ews3a = subset(ews3,subset=(t > 62 & t < 80))
tQ = which.max(ews3a$Q)
ewarQ = ews3a$t[tQ]
print(c('forecast lead time, days ',kh.day),quote=F)
print(c('time of max(forecast var) ',ewarQ,', difference from crit time = ',dpoint-ewarQ),quote=F)

#save(epsilon,sigma,dpoint,ewar1,ewarsd1,ewarmv,ewarQ,file='kNCsim_EWS_eps1.Rdata')
