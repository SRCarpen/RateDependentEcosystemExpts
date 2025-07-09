# Fit DOC-Grazer-Chl model to 2024 data, December 2024
# Modified Fit2015_Zt-by-reg+PubFigs_2+3_drivers_2022-01-28:
#   * remove exp(-b*A) term
#   * add a refuge term f(K - A) as a low-A stabilizer
#   * use q=1, ordinary Michaelis-Menten grazing (see notes ~ how this affects eq)

rm(list = ls())
graphics.off()

library(moments)
library(runner)

#source('Constants+Functions_PAR-P-Z_2025-01-16.R')
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

# Retrieve multi-lake regression model to predict Zt 
#  by Zmix=f(DOC_etc)_2021-12-31.R
#load(file='ZmixModel.Rdata')  # 7 lakes
#load(file='ZmixModel_w_Hbird.Rdata')  # 8 lakes
load(file='ZmixModel_LakeYears_Hbird.Rdata')  # lake years
print(summary(reg1)) 
#pars = reg1$coefficients
#
# Use median squarea sqrt(lake area in ha) values in fetch

# Test growth and grazing curves

# Calculate growth curves and eq for Peter Lake
p0 = 1 # choose P load
# use mean drivers
# Peter was Aquashaded in 2024 so use Paul Lake kNC0 for 'unmanipulated Peter'
kNC.test = quantile(kNC0[1:77],probs=0.5)  # kNC0[1:77] is Paul Lake data
DOC.test = kNC.test/epsDOC
RZT = median(Ztvec)
RZbiom = quantile(ZBvec,probs=0.75)

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

# Equilibria versus kNC over observed gradient =========================================
q.kNC = unname(quantile(kNC0,probs = c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)) )
klo = q.kNC[1]
khi = max(kNC0) #q.kNC[9]
nk = 25
kvec = seq(klo,khi,length.out=nk)  # observed range of kNC
#kvec = seq(1.5,1.51,length.out=nk)  # narrow range to define the critical point
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

# convert equilibria to chl
eqmat = eqmat/CChl
windows()
par(mfrow=c(1,1),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
yrange = range(eqmat,na.rm=T)
plot(kvec,eqmat[,1],ylim=yrange,type='p',pch=19,col='blue',
     xlab='Non-Chl PAR extinction',ylab='equilibria, mg Chl / m^3',
     main=paste('equilibria, P load = ',p0)) 
points(kvec,eqmat[,2],type='p',lwd=2,pch=21,cex=1.2,col='red')
points(kvec,eqmat[,3],type='p',pch=19,col='forestgreen')
grid(lwd=2,col='slategray')

# TEST -- FIND CRITICAL p0 and kNC =================================================

# find critical p0
cvec = as.vector(c(kNC.test,RZT,RZbiom,rmax,hmax)) # vector of constants
print('vector of constants',quote=F)
print(cvec)
guess = log(c(2,0.6))  # guesses of Chl & p0 critical value
critfit = optim(guess,crit.pload,gr=NULL,cvec,method='Nelder-Mead')
#critfit = optim(guess,crit.pload,gr=NULL,cvec,method='BFGS') # experiment with method
parfit = exp(critfit$par)
print(c('critical AChl and p0 ',parfit),quote=F)
p0.crit = parfit[2]
print(c('convergence',critfit$convergence),quote=F)

# find critical kNC
cvec = as.vector(c(p0,RZT,RZbiom,rmax,hmax)) # vector of constants
print('vector of constants',quote=F)
print(cvec)
guess = c(5.6,1.507)  # guesses of Chl & kNC critical value
critfit = optim(guess,crit.kNC,gr=NULL,cvec,method='Nelder-Mead')
parfit = critfit$par
print(c('critical AChl and kNC ',parfit),quote=F)
kNC.crit = parfit[2]
print(c('convergence',critfit$convergence),quote=F)

# END TEST TO FIND CRITICAL values ===============================================

# PLOT THE 2 SETS OF STABLE POINTS SIDE BY SIDE
windows(width=10,height=5)
par(mfrow=c(1,2),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
yrange = range(peqmat,na.rm=T)
plot(p0grad,peqmat[,1],ylim=yrange,type='b',lwd=2,pch=19,col='blue',
     xlab='P load, mg/(m^2 d)',ylab='equilibria, mg Chl / m^3',
     main=paste('equilibria, kNC = ',round(kNC.test,3))) 
points(p0grad,peqmat[,2],type='b',lwd=2,pch=21,cex=1.2,col='red')
points(p0grad,peqmat[,3],type='b',lwd=2,pch=19,col='forestgreen')
grid(lwd=1,col='slategray')
abline(v=p0.crit,lty=3,lwd=3,col='black')

yrange = range(eqmat,na.rm=T)
plot(kvec,eqmat[,1],ylim=yrange,type='p',pch=19,col='blue',
     xlab='Non-Chl PAR extinction',ylab='equilibria, mg Chl / m^3',
     main=paste('equilibria, P load = ',p0)) 
points(kvec,eqmat[,2],type='p',lwd=2,pch=21,cex=1.2,col='red')
points(kvec,eqmat[,3],type='p',pch=19,col='forestgreen')
grid(lwd=1,col='slategray')
abline(v=kNC.crit,lty=3,lwd=3,col='black')

# CALCULATE CRITICAL VALUES OVER GRADIENTS ========================================
# Calculate critical p0 for a range of kNC
kNC.v = c(0.1, 0.25, 0.5, 1, 1.5, 1.75, 2, 2.1)
nknc = length(kNC.v)
critp0.v = rep(0,nknc)
for(i in 1:nknc) {
  knci = kNC.v[i]
  cvec = as.vector(c(knci,RZT,RZbiom,rmax,hmax)) # vector of constants
  #print('vector of constants',quote=F)
  #print(cvec)
  guess = log(c(3,0.53))  # guesses of Chl & p0 critical value
  critfit = optim(guess,crit.pload,gr=NULL,cvec,method='Nelder-Mead')
  parfit = exp(critfit$par)
  #print(c('critical AChl and p0 ',parfit),quote=F)
  critp0.v[i] = parfit[2]
}

# calculate critical kNC for a range of p0
p0.v = c(0.3,0.6,1,1.5,2,2.5,3)
np0 = length(p0.v)
critknc.v = rep(0,np0)
for(i in 1:np0) {
  p0i = p0.v[i]
  cvec = as.vector(c(p0i,RZT,RZbiom,rmax,hmax)) # vector of constants
  #print('vector of constants',quote=F)
  #print(cvec)
  guess = c(5.6,1.507)  # guesses of Chl & kNC critical value
  critfit = optim(guess,crit.kNC,gr=NULL,cvec,method='Nelder-Mead')
  parfit = critfit$par
  #print(c('critical AChl and kNC ',parfit),quote=F)
  critknc.v[i] = parfit[2]
  print(c('convergence',critfit$convergence),quote=F)
}

# plot critical values
windows(height=5,width=10)
par(mfrow=c(1,2),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(kNC.v,critp0.v,type='b',lwd=2,pch=20,cex=1,col='blue',xlab='kNC',
     ylab='critical Pload',main='P load critical values')
grid()
text(x=0.8,y=1.2,'2 stable eq',cex=1.6,font=2,col='mediumseagreen')
text(x=1.7,y=0.6,'1 stable eq',cex=1.6,font=2,col='dodgerblue')
plot(p0.v,critknc.v,type='b',lwd=2,pch=20,cex=1,col='sienna',xlab='P Load',
     ylab='critical kNC',main='kNC critical values')
grid()
text(x=1.7,y=1,'2 stable eq',cex=1.6,font=2,col='mediumseagreen')
text(x=1,y=3,'1 stable eq',cex=1.6,font=2,col='dodgerblue')
#abline(h=2.1,lty=3,lwd=3,col='red')

# plot the curves on the same axes
windows()
par(mfrow=c(1,1),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(p0.v,critknc.v,type='p',lwd=2,pch=19,cex=1.2,col='sienna',xlab='P Load',
     ylab='kNC',main='critical values')
points(critp0.v,kNC.v,type='p',lwd=2,pch=17,cex=1.2,col='blue')
grid()
abline(h=2.1,lty=3,lwd=3,col='red')
text(x=2,y=1.5,'2 stable eq',cex=1.6,font=2,col='mediumseagreen')
#text(x=0.15,y=0.15,'1 eq',cex=1.6,font=2,col='forestgreen')
#text(x=2,y=2,'1 stable eq',cex=1.6,font=2,col='blue')
text(x=1,y=3,'1 stable eq',cex=1.6,font=2,col='brown')


# simulate and calculate EWS over a gradient of p0 ===================================

# simulate on Pload gradient
dt = 0.04  # days
dtn = sqrt(dt)
sigma = 1 
nstep = 200/dt  # days to simulate
Achl = rep(0,nstep)  # vector to hold results
#Achl[1] = median(Chl0[78:172])  # Peter Lake is 78 to 172
Achl[1] = 2  # start above unstable eq to grow
# Driver gradient to simulate
p0sim = seq(0.3,1,length.out=nstep)
kNCsim = seq(1,2,length.out=nstep)
# fixed values, in cases where we need them
kNCfix = median(kNC0[1:77])  # Paul Lake is 1 to 77
#kNCfix = median(kNC0[78:172])  # Peter Lake is 78 to 172
kNCmed = kNCfix # set the 'average' kNC for this run
p0fix = 0.3
# ecosystem state to simulate
#Chl,kNC,p0,ZT,Zb,rmax,hmax,sigma,dt,dtn
ZTmed = median(Zt0[78:172])  # Peter Lake is 78 to 172
ZBmed = median(ZB0[78:172])
#
# simulate
tstep = c(1:nstep)*dt
noise = rnorm(nstep)  # N(0,1) noise
for(i in 2:nstep) {
  Achl[i] = dAdt(Achl[i-1],kNCfix,p0sim[i-1],ZTmed,ZBmed,rmax,hmax,noise[i-1],sigma,dt,dtn)
}

winlen = 25
v.Achl = runner(Achl,k=winlen,f=var)

windows(height=12,width=6)
par(mfrow=c(4,1),mar=c(4,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(tstep,p0sim,type='l',lwd=1,col='darkred')
abline(h=0.6,lty=3,lwd=3,col='orangered')
grid()
plot(tstep,Achl,type='l',lwd=1,col='forestgreen')
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
print('b matrix',quote=F)
print(bmat)
#print('eigenvalues of b matrix',quote=F)
#beig = eigen(bmat,only.values=T)
#print(beig$values)
#print('prediction error variances',quote=F)
yhat = Xmat0%*%bmat
err = Xmat1 - yhat
verr = apply(err,2,var)
#print(verr)
parcov = Xinv*verr
#
# set up DLM
delta = 0.99
n.gamma = 1
d.gamma = verr
dlm.est = DLM(delta,n.gamma,d.gamma,bmat,parcov,Xmat1,Xmat0)
ar1 = dlm.est[[3]]  # parameter
sd1 = dlm.est[[4]]  # sd; these are small, ~ 1/10^8
# add ar1 values to plot
plot(tstep[2:nX],ar1,type='l',lwd=1,col='blue',xlab='tstep',ylab='AR(1) coef')
grid()
abline(h=1)

# simulate on kNC gradient ==========================================================
dt = 0.04  # days
dtn = sqrt(dt)
sigma = 1
nstep = 200/dt  # days to simulate
Achl = rep(0,nstep)  # vector to hold results
#Achl[1] = median(Chl0[78:172])  # Peter Lake is 78 to 172
Achl[1] = 10  # start near upper stable limb
#Achl[1] = 2 # start near lower stable limb
# Driver gradient to simulate
p0sim = seq(0.3,2,length.out=nstep)
#kNCsim = seq(1,2,length.out=nstep)  # linear kNC gradient
# simulated kNC using Ritchie-Sieber variable speed model
#    /2024_PAR-DOC-Zoop-P-Chl_model/Test2_RS_speed_control.R
#save(speedvec,kNC.half,kNC.full,file='simkNC.Rdata')
# speedvec is a vector of speeds as multiples of epsilon
# kNC.half simulates kNC from 1 to 2 using the second half 
#  of R-S solution
# kNC.full simulates the full R-S model for kNC between 0 and 2
# in both matrices each column is a speed
load(file='simkNC.Rdata')
kNCsim = kNC.full[,3]   # select column for speed of change in kNC
# fixed values, in cases where we need them
kNCmed = kNCfix
p0fix = 1
# ecosystem state to simulate
#Chl,kNC,p0,ZT,Zb,rmax,hmax,sigma,dt,dtn
ZTmed = median(Zt0[78:172])  # Peter Lake is 78 to 172
ZBmed = median(ZB0[78:172])

#
# simulate
tstep = c(1:nstep)*dt
noise = rnorm(nstep)  # N(0,1) noise
for(i in 2:nstep) {
  Achl[i] = dAdt(Achl[i-1],kNCsim[i-1],p0fix,ZTmed,ZBmed,rmax,hmax,noise[i-1],sigma,dt,dtn)
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
print('b matrix',quote=F)
print(bmat)
#print('eigenvalues of b matrix',quote=F)
#beig = eigen(bmat,only.values=T)
#print(beig$values)
#print('prediction error variances',quote=F)
yhat = Xmat0%*%bmat
err = Xmat1 - yhat
verr = apply(err,2,var)
#print(verr)
parcov = Xinv*verr
#
# set up DLM
delta = 0.99
n.gamma = 1
d.gamma = verr
dlm.est = DLM(delta,n.gamma,d.gamma,bmat,parcov,Xmat1,Xmat0)
#DLM.out <- list(predix,varpredix,pars,parvar,Svec) # output line of DLM
yhat = dlm.est[[1]]
modvar = dlm.est[[2]]  # process or 'model' variance
ar1 = dlm.est[[3]]  # parameter
sd1 = dlm.est[[4]]  # sd; these atre small, ~ 1/10^8
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
plot(tstep[2:nX],log10(vstat),type='l',lwd=1,col='purple',xlab='tstep',ylab='log10(stationary var.)',
     main='NA when ar1^2 > 1')
grid()

