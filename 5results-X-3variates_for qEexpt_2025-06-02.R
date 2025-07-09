# Build composite plots of 5 results X 3 speeds with layout() for 1 of the 3 experiments

rm(list = ls())
graphics.off()

# Descriptions of input files - examples from each run ALL RUNS EPSILON = 1 
# file name abbreviations:  Chl, Hvore, Plvore
# Fname = c('qE_eps1_Plvore_EWS.Rdata')  # file name example

# SAVE LINES
# save(epsilon,sigma,dpoint,tstep,qEsim,qE.crit,nX,At,Ft,Ht,Pt,Jt,Xcen,kh,kh.day,yhat,
#     ar1,sd1,modvar,tfore,Qmat,file=Fname)

# ------------------------------------------------------------------------------------------------------

# layout example
## divide device into two rows and two columns
## allocate figure 1 and figure 2 as above
## respect relations between widths and heights
#nf <- layout(matrix(c(1,1,0,2), 2, 2, byrow = TRUE), respect = T)
#layout.show(nf)

# mfrow is incompatible with layout!!

windows(width=10,height=15)
par(mar=c(4,4.5,2,1)+0.1,cex.axis=1.4,cex.lab=1.6,cex.main=1.6,font.axis=2,
    font.lab=2,font.main=2)
nf = layout(matrix(c(1:15),nr=5,nc=3, byrow = F),respect=T)
layout.show(nf)
#
load(file='qEsim_EWS_Chl_eps1.Rdata')
# 1, Chl
plot(tstep,At,type='l',lwd=2,col='forestgreen',xlab='Day',
     ylab='Chlorophyll',main='Phytoplankton')
#points(tstep[2:nX],yhat,type='l',lwd=2,col='black')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 2, AR(1)
plot(tstep[2:nX],ar1,type='l',lwd=2,col='darkblue',xlab='Day',
     ylab='AR(1) coef')
grid()
abline(h=1,lty=3,lwd=3,col='black')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 3, sd of AR1
plot(tstep[2:nX],sd1,log='y',type='l',lwd=2,col='blue',xlab='Day',
     ylab='s.e. of AR(1) coef)')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 4, model variance
plot(tstep[2:nX],modvar,log='y',type='l',lwd=2,col='magenta',xlab='Day',
     ylab='Model Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 5, prediction variance
plot(tfore,Qmat,log='y',type='l',lwd=2,col='purple',xlab='Day',
     ylab='Forecast Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
#text(x=100,y=0.99*max(Qmat),paste('days ahead = ',kh.day),cex=1.8)
#
# -------------------------------------------------------------------------------
load(file='qEsim_EWS_Zoop_eps1.Rdata')
# 2, zoopl
plot(tstep,Ht,type='l',lwd=2,col='forestgreen',xlab='Day',
     ylab='Biomass',main='Herbivore')
#points(tstep[2:nX],yhat,type='l',lwd=2,col='black')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 2, AR(1)
plot(tstep[2:nX],ar1,type='l',lwd=2,col='darkblue',xlab='Day',
     ylab='AR(1) coef')
grid()
abline(h=1,lty=3,lwd=3,col='black')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 3, sd of AR1
plot(tstep[2:nX],sd1,log='y',type='l',lwd=2,col='blue',xlab='Day',
     ylab='s.e. of AR(1) coef)')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 4, model variance
plot(tstep[2:nX],modvar,log='y',type='l',lwd=2,col='magenta',xlab='Day',
     ylab='Model Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 5, prediction variance
plot(tfore,Qmat,log='y',type='l',lwd=2,col='purple',xlab='Day',
     ylab='Forecast Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
#text(x=100,y=0.99*max(Qmat),paste('days ahead = ',kh.day),cex=1.8)
#
# -------------------------------------------------------------------------------
load(file='qEsim_EWS_Plvore_eps1.Rdata')
# 1, Chl
plot(tstep,Pt,type='l',lwd=2,col='forestgreen',xlab='Day',
     ylab='Density',main='Planktivore')
#points(tstep[2:nX],yhat,type='l',lwd=2,col='black')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 2, AR(1)
plot(tstep[2:nX],ar1,type='l',lwd=2,col='darkblue',xlab='Day',
     ylab='AR(1) coef')
grid()
abline(h=1,lty=3,lwd=3,col='black')
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 3, sd of AR1
plot(tstep[2:nX],sd1,log='y',type='l',lwd=2,col='blue',xlab='Day',
     ylab='s.e. of AR(1) coef)')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 4, model variance
plot(tstep[2:nX],modvar,log='y',type='l',lwd=2,col='magenta',xlab='Day',
     ylab='Model Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
# 5, prediction variance
plot(tfore,Qmat,log='y',type='l',lwd=2,col='purple',xlab='Day',
     ylab='Forecast Variance')
grid()
abline(v=dpoint,lty=2,lwd=2,col='darkred')
#text(x=100,y=0.99*max(Qmat),paste('days ahead = ',kh.day),cex=1.8)
