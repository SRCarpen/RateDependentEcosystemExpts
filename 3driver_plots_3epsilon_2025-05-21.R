# Build composite plots of drivers, 3 experiments & 3 speeds with layout()

rm(list = ls())
graphics.off()

# Descriptions of input files - examples from each experiment 
#  (3 files per expt, each for a different rate)

# save the simulation ----------------------------------------------------------------------------------
#save(epsilon,sigma,dpoint,Pload0,tstep,p0sim,p0.crit,nX,Achl,kh,kh.day,Acen,yhat,
#     ar1,sd1,modvar,tfore,Qmat,file='p0sim_EWS_eps3.Rdata')
# ------------------------------------------------------------------------------------------------------

# save the simulation ----------------------------------------------------------------------------------
#save(epsilon,sigma,dpoint,Pload0,tstep,kNCsim,kNC.crit,nX,Achl,kh,kh.day,Acen,yhat,
#     ar1,sd1,modvar,tfore,Qmat,file='kNCsim_EWS_eps0dot3.Rdata')
# ------------------------------------------------------------------------------------------------------

# save the simulation ----------------------------------------------------------------------------------
#save(epsilon,sigma,dpoint,tstep,qEsim,qE.crit,nX,At,Ft,Ht,Pt,Jt,Xcen,kh,kh.day,yhat,
#     ar1,sd1,modvar,tfore,Qmat,file='qEsim_EWS_eps0dot3.Rdata')
# ------------------------------------------------------------------------------------------------------

# layout example
## divide device into two rows and two columns
## allocate figure 1 and figure 2 as above
## respect relations between widths and heights
#nf <- layout(matrix(c(1,1,0,2), 2, 2, byrow = TRUE), respect = T)
#layout.show(nf)

# mfrow is incompatible with layout!!

windows(width=12,height=5)
par(mar=c(4,4.5,2,1)+0.1,cex.axis=1.6,cex.lab=1.6)
nf = layout(matrix(c(1,2,3),nr=1,nc=3, byrow = T),respect=T)
#
load(file='p0sim_EWS_eps0dot3.Rdata')
p0sim03 = p0sim
plot(tstep,p0sim03,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='p0 Simulation',
     main='epsilon = 0.3')
abline(h=p0.crit,lty=3,lwd=3,col='orangered')
grid()
#
load(file='p0sim_EWS_eps1.Rdata')
p0sim1 = p0sim
plot(tstep,p0sim1,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='p0 Simulation',
     main='epsilon = 1')
abline(h=p0.crit,lty=3,lwd=3,col='orangered')
grid()
#
load(file='p0sim_EWS_eps3.Rdata')
p0sim3 = p0sim
plot(tstep,p0sim3,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='p0 Simulation',
     main='epsilon = 3')
abline(h=p0.crit,lty=3,lwd=3,col='orangered')
grid()

# load kNC experiment 
windows(width=12,height=5)
par(mar=c(4,4.5,2,1)+0.1,cex.axis=1.6,cex.lab=1.6)
nf = layout(matrix(c(1,2,3),nr=1,nc=3, byrow = T),respect=T)
#
load(file='kNCsim_EWS_eps0dot3.Rdata')
kNCsim03 = kNCsim
plot(tstep,kNCsim03,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='kNC Simulation',
     main='epsilon = 0.3')
abline(h=p0.crit,lty=3,lwd=3,col='orangered')
grid()
#
load(file='kNCsim_EWS_eps1.Rdata')
kNCsim1 = kNCsim
plot(tstep,kNCsim1,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='kNC Simulation',
     main='epsilon = 1')
abline(h=p0.crit,lty=3,lwd=3,col='orangered')
grid()
#
load(file='kNCsim_EWS_eps3.Rdata')
kNCsim3 = kNCsim
plot(tstep,kNCsim3,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='kNC Simulation',
     main='epsilon = 3')
abline(h=p0.crit,lty=3,lwd=3,col='orangered')
grid()

# save dense tstep with a different name
tstepd = tstep

# load qE experiment
load(file='qEsim_EWS_eps0dot3.Rdata')
tstepq = tstep
qEsim03 = qEsim

windows(width=12,height=5)
par(mar=c(4,4.5,2,1)+0.1,cex.axis=1.6,cex.lab=1.6)
nf = layout(matrix(c(1,2,3),nr=1,nc=3, byrow = T),respect=T)
plot(tstepq,qEsim03,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='qE Simulation',
     main='epsilon = 0.3')
abline(h=qE.crit,lty=3,lwd=3,col='orangered')
grid()
#
load(file='qEsim_EWS_eps1.Rdata')
qEsim1 = qEsim
plot(tstepq,qEsim1,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='qE Simulation',
     main='epsilon = 1')
abline(h=qE.crit,lty=3,lwd=3,col='orangered')
grid()
#
load(file='qEsim_EWS_eps3.Rdata')
qEsim3 = qEsim
plot(tstepq,qEsim3,type='l',lwd=2,col='sienna',xlab = 'Day',ylab='qE Simulation',
     main='epsilon = 3')
abline(h=qE.crit,lty=3,lwd=3,col='orangered')
grid()

# Try all 3 experiments together
windows(width=12,height=5)
par(mar=c(4,4.5,2,1)+0.1,cex.axis=1.6,cex.lab=1.8,cex.main=1.8)
nf = layout(matrix(c(1,2,3),nr=1,nc=3, byrow = T),respect=T)
# p0 experiment
yrange = range(c(p0sim03,p0sim1,p0sim3))
plot(tstepd,p0sim03,ylim=yrange,type='l',lwd=2,col='blue',
     xlab = 'Day',ylab='Enrichment rate, mg P / (m^2 day)',
     main = 'Nutrient Enrichment Experiment')
points(tstepd,p0sim1,type='l',lwd=2,col='purple')
points(tstepd,p0sim3,type='l',lwd=2,col='red')
abline(h=p0.crit,lty=3,lwd=3,col='black')
grid()
legend('topleft',legend=c('epsilon = 0.3','epsilon = 1','epsilon=3'),
        lty=c(1,1,1),lwd=c(2,2,2),
       col=c('blue','purple','red'),cex=1.8,bty='n')
# kNC experiment
yrange = range(c(kNCsim03,kNCsim1,kNCsim3))
plot(tstepd,kNCsim03,ylim=yrange,type='l',lwd=2,col='blue',
     xlab = 'Day',ylab='kNC, 1/m',
     main = 'PAR Extinction Experiment')
points(tstepd,kNCsim1,type='l',lwd=2,col='purple')
points(tstepd,kNCsim3,type='l',lwd=2,col='red')
abline(h=kNC.crit,lty=3,lwd=3,col='black')
grid()
# qE Experiment
yrange = range(c(qEsim03,qEsim1,qEsim3))
plot(tstepq,qEsim03,ylim=yrange,type='l',lwd=2,col='blue',
     xlab = 'Day',ylab='qE, 1/d',
     main = 'Piscivore Harvest Experiment')
points(tstepq,qEsim1,type='l',lwd=2,col='purple')
points(tstepq,qEsim3,type='l',lwd=2,col='red')
abline(h=qE.crit,lty=3,lwd=3,col='black')
grid()
