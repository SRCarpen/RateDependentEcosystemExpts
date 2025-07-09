# Run a univariate forecast DLM for qE cascade model 2025-05-17
# /2025_RateDependence_kNC_P_qE/qE_rate_expts_Cascade_EcoLetts2008_V2_2025-05-15.R

rm(list = ls())
graphics.off()

# Constants & functions are collected in 
#  /2025_RateDependence_kNC_P_qE/Constants+Functions_qE_Cascade_2025-05-15.R
source('Constants+Functions_qE_Cascade_2025-05-15.R')

# Save data for further analysis and plots
#save(tstep,nsim,epsilon,qEsim,qE.crit,dpoint,At,Ft,Ht,Pt,Jt,file='Sim_qE_CascadeSqueal.Rdata')

load(file='Sim_qE_CascadeSqueal_eps0dot3_lowsigma.Rdata')
Fname = c('qEsim_EWS_lowsig_eps0.3.Rdata')  # name of output file

windows()
par(mfrow=c(3,1),mar=c(2.5,4.5,1,2)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(tstep,Ft,type='l',lwd=2,col='darkgreen',xlab='time step',ylab='Planktivores')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
plot(tstep,Ht,type='l',lwd=2,col='blue',xlab='time step',ylab='Herbivores')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
plot(tstep,Pt,type='l',lwd=2,col='forestgreen',xlab='time step',ylab='Phytoplankton')
abline(v=dpoint,lty=2,lwd=2,col='darkred')

# select variate for forecast DLM analysis
Xsel = Pt  
#Fname = c('qE_eps1_Chl_EWS.Rdata')  # file name to save EWS for plotting
fcstvar = c('Chlorophyll')
print(c('epsilon ',epsilon,' day of crit trans ',round(dpoint,1),' variate= ',fcstvar),quote=F)

# run 1-dimensional DLM on selected variate
# DLM <- function(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat) 
Xcen = Xsel - mean(Xsel)  # center the data
npar = 1   # if matrix is centered we don't need the intercept
nX = length(Xcen)
unit = rep(1,nX)

#
# Quick parameter estimate using LS 
print('Multivariate least squares regression to check the data',quote=F)
Xmat1 = as.matrix(Xcen[2:nX])
Xmat0 = as.matrix(Xcen[1:(nX-1)])
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
delta = 0.95
n.gamma = 1
d.gamma = verr
dt.day = 1 # timestep is 1 for qE simulations
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
vstat = ifelse(abs(ar1^2) <= 1,modvar/(1 - ar1^2),0)  # stationary variance
plot(tstep[2:nX],log10(vstat),type='l',lwd=1,col='black',xlab='tstep',ylab='log10(stationary var.)',
     main='NA when ar1^2 > 1')
grid()

# reorganize key DLM results ===============================================================

windows(height=12,width=12)
par(mfrow=c(3,2),mar=c(4,4.5,1,2)+0.1,cex.axis=1.8,cex.lab=1.8)

plot(tstep,qEsim,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='qE Simulation')
abline(h=qE.crit,lty=3,lwd=3,col='orangered')
text(x=20,y=1.8,paste('epsilon = ',epsilon),cex=1.8)
grid()

plot(tstep,Xcen,type='l',lwd=4,col='skyblue',xlab='Day',
     ylab=fcstvar)
points(tstep[2:nX],yhat,type='l',lwd=1,col='black')
grid()
abline(v=dpoint,lty=2,lwd=2,col='red')
#text(x=100,y=0.65,paste('sigma = ',round(sigma,1)),cex=1.8)
legend('bottomleft',legend=c('simulation','DLM'),col=c('skyblue','black'),
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
tfore = tstep[(kh+1):(nX-2)]
Qmat = fcast[(kh+1):(nX-2),4]
# Plot forecast vs obs
#yrange = range(c(fest,fobs),na.rm=T)
#plot(fest,fobs,type='p',pch=20,cex=0.5,col='black',xlab='mean of 250-step forecast',
#     ylab='Observation')
#grid()
plot(tfore,Qmat,log='y',type='l',lwd=1,col='purple',xlab='Day',
     ylab='Forecast Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='red')
text(x=100,y=0.97*max(Qmat),paste('days ahead = ',kh.day),cex=1.8)
#text(x=100,y=0.85*max(Qmat),paste(c(kh/dt.day,' days')),cex=1.8)

# save the simulation ----------------------------------------------------------------------------------
save(epsilon,sigma,dpoint,tstep,qEsim,qE.crit,nX,At,Ft,Ht,Pt,Jt,Xcen,kh,kh.day,yhat,
     ar1,sd1,modvar,tfore,Qmat,file=Fname)
# ------------------------------------------------------------------------------------------------------
