# qE rate-dependence expts for Cascade model squeal 1
# 
# This version derived from
# /2025_RateDependence_kNC_P_qE/Cascade_EcoLetts2008_model_forOikosQD2014_2013-01-04.R
#
# Original model was published in
# Carpenter, S.R., W.A. Brock, J.J. Cole, J.F. Kitchell and M.L. Pace. 2008. 
# Leading indicators of trophic cascades. Ecology Letters 11: 128-138.
# T
# Here we modify the version of the model published in
# Brock, W.A. and S.R. Carpenter. 2012. Early warnings of regime shift when the ecosystem structure is unknown. 
# PLoS ONE 7(9): e45586. doi:10.1371/journal.pone.0045586
#
# Original program by SRC 26 Nov 2012
# Annotated for Oikos by SRC on 4 January 2013
# Modified for rate dependence SRC 15 May 2025

rm(list = ls())
graphics.off()

#library('zoo')

# Constants & functions are collected in 
#  /2025_RateDependence_kNC_P_qE/Constants+Functions_qE_Cascade_2025-05-15.R
source('Constants+Functions_qE_Cascade_2025-05-15.R')

# load critical values from 
# /2025_RateDependence_kNC_P_qE/Cascade_crit_qE_A_F_2025-05-12.R
# save(CritPt,Feq,file='Crit_qE_AdultPisc.Rdata')
load(file='Crit_qE_AdultPisc.Rdata')
qE.crit = CritPt$par[1]
Acrit = CritPt$par[2]
print(c('qE.crit = ',qE.crit,', Acrit = ',Acrit,', F at Acrit = ',Feq),quote=F)

nint <- 30  # Time steps per 'year' 
dt <- 1/nint
dtZ <- sqrt(dt)

nsim = 120
tstep = (1:nsim)


# Set up qEvec for rate dependence
# original from 2008 paper:
#qEvec <- c(rep(0.01,nburn),seq(0.01, 0.05, length.out=(nsim)))
# constant qE
#qEvec = rep(0.001,nburn)

# simulate qE over a gradient that includes the critical value
# this function simulates Ritchie-Sieber lamda given epsilon and nday
# We set lamda = qEsim
lamsim = function(epsilon,dt,nsim,qE.crit) {
  nhalf = nsim
  halfstep = c(0:nhalf)*dt
  tlam = c( -1*rev(halfstep),halfstep[2:nhalf] )
  nlam = length(tlam)
  # calculate and plot lamda versus tstep
  lamdamax=2  # lamda ranges 0 to lamdamax
  lamda = (lamdamax/2)*( tanh(lamdamax*epsilon*tlam/2) +1)
  # plot lamda vs tstep
  windows(width=8,height=4)
  par(mfrow=c(1,2),mar=c(4.5,4.5,3,2)+0.1,cex.axis=1.6,cex.lab=1.6)
  plot(tlam,lamda,type='l',lwd=2,col='sienna',
       xlab = 'time step',ylab='lamda',main='Ritchie-Sieber simulation of lamda')
  abline(h=qE.crit,lty=1,lwd=1,col='darkblue')
  grid()
  # rescale time axis to match the chlorophyll simulation
  lamfil = spline(x=tlam,y=lamda,method='fmm',n=(1/dt)*(length(halfstep)-1))
  nlam = length(lamfil$x)
  print(c('n lamda simulated',nlam),quote=F)
  lamstep = c(1:nlam)
  lamda = lamfil$y
  lamday = lamstep*dt
  plot(lamday,lamfil$y,type='l',lwd=1,col='brown',
       xlab='time rescaled',ylab='lamda',
       main='Ritchie-Sieber with time rescaled')
  abline(h=qE.crit,lty=1,lwd=1,col='darkblue')
  grid()
  # find time where qE = qE.crit
  delt = lamda - qE.crit
  sdelt = sign(delt)
  dsdelt = c(0,-diff(sdelt))
  scross = lamstep[which(!dsdelt == 0)]
  dcross = lamday[which(!dsdelt == 0)]
  print(c('time step and day when lamda = qE.crit',scross,dcross),quote=F)
  outlist = list(lamstep,lamday,lamda,scross,dcross)
  return(outlist)
}

# run qE simulation
# outlist = list(lamstep,lamday,lamda,scross,dcross)
print('',quote=F)
print('Report of qE simulation',quote=F)
epsilon = 1  # if epsilon < ~0.25 then the critical qE is not reached
print(c('epsilon for rate experiment',epsilon),quote=F)
dt.qE = dt  # time step for qE simulation
print(c('time step for qE',dt.qE),quote=F)
print(c('time step for food web',dt),quote=F)
qEvals = lamsim(epsilon,dt.qE,nsim,qE.crit)
qEstep = qEvals[[1]]
#qEday = qEvals[[2]]
qEsim.all = qEvals[[3]]
# subsample qEsim.all to the length of nsim
n.all = length(qEsim.all)
ystep = n.all/nsim
ikeep = seq(1,n.all,by=ystep)
qEsim = qEsim.all[ikeep]  # qE series matching nsim time steps
tpoint = qEvals[[4]]
dpoint = qEvals[[5]]
print('End of qE simulation',quote=F)
print('',quote=F)

# Set up vectors to hold simulation results

# Food web:
At = rep(0,nsim)
Ft = At
Jt = At
Ht = At
Pt = At
At[1] = Ainit
Ft[1] = Finit
Jt[1] = fA*Ainit
Ht[1] = Hinit
Pt[1] = Pinit

# food web simulation

# noise used by FWsim.step; FWsim.iter has endogenous noise
#noise.vec = rnorm(3*nsim)
#noise.mat = matrix(noise.vec,nrow=nsim,ncol=3)

for(i in 2:nsim)  {
  # daily time step
  #FWnext = FWsim.step(qEsim[i-1],At[i-1],Ft[i-1],Jt[i-1],Ht[i-1],Pt[i-1],dt,dtZ,noise.mat[(i-1),])
  # nint time steps per day
  FWnext = FWsim.iter(qEsim[i-1],At[i-1],Ft[i-1],Jt[i-1],Ht[i-1],Pt[i-1],dt,dtZ,nint) 
  At[i] = FWnext[[1]]
  Ft[i] = FWnext[[2]]
  Jt[i] = FWnext[[3]]
  Ht[i] = FWnext[[4]]
  Pt[i] = FWnext[[5]]
}

windows()
par(mfrow=c(3,1),mar=c(2.5,4.5,1,2)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(tstep,qEsim,type='l',lwd=2,col='darkgreen',xlab='time step',ylab='qE')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
plot(tstep,At,type='l',lwd=2,col='blue',xlab='time step',ylab='Adults')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
plot(tstep,Jt,type='l',lwd=2,col='cyan',xlab='time step',ylab='Juveniles')
abline(v=dpoint,lty=2,lwd=2,col='darkred')

windows()
par(mfrow=c(3,1),mar=c(2.5,4.5,1,2)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(tstep,Ft,type='l',lwd=2,col='darkgreen',xlab='time step',ylab='Planktivores')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
plot(tstep,Ht,type='l',lwd=2,col='blue',xlab='time step',ylab='Herbivores')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
plot(tstep,Pt,type='l',lwd=2,col='forestgreen',xlab='time step',ylab='Phytoplankton')
abline(v=dpoint,lty=2,lwd=2,col='darkred')

# Save data for further analysis and plots
save(tstep,nsim,epsilon,qEsim,qE.crit,dpoint,At,Ft,Ht,Pt,Jt,file='Sim_qE_CascadeSqueal_eps1.Rdata')
