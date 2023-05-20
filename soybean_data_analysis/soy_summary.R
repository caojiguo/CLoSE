library(latex2exp)
# generate Figure 3 of the manuscript: a sample of daily minimum and daily maximum temperature curves.
mtx <- readRDS('soy_x.RDS')
X2 = mtx[,30:335] # daily minimum temperature
X1 = mtx[,395:700] # daily maximum temperature
plot(c(30:335),X2[1,],type='l',xlab='Month',ylab='',cex.lab=2,cex.axis=2,cex.main=2,xaxt='n',lwd=2, col='grey10',ylim=c(min(X2),max(X1)))
for (k in 2:50){
  lines(c(30:335),X2[k,],col='grey10')
}
for (k in 1:50){
  lines(c(30:335),X1[k,],col='gray80')
}
lines(c(30:335),apply(X2,2,mean), lwd=3,col='red')
lines(c(30:335),apply(X1,2,mean), lwd=3,col='blue')
axis(1, at=c(30,61,91,122,152,183,213,243,273,304,335), lab=c("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),cex.axis=2)
legend(20,30,legend="mean of MIN",col='red', lwd=2 ,lty=1,cex=2,bty='n',x.intersp = 0.5, pt.cex=0.5)
legend(150,30,legend="mean of MAX",col='blue', lwd=2 ,lty=1,cex=2,bty='n',x.intersp = 0.5, pt.cex=0.5)


# generate Figure 4(a) and 4(b) of the manuscript
setwd('~/r_functions/soybean_data_analysis/0.25-quantile/')
# estimations obtained from bootstrap samples
bhat2s <- bhat1s <- matrix(NA,nrow=55,ncol=100)
for (k in 1:100){
  tmp <- readRDS(paste0(k,'bHat1.RDS'))
  bhat1s[,k]<- tmp
  tmp <- readRDS(paste0(k,'bHat2.RDS'))
  bhat2s[,k]<- tmp
}
thetahat <- matrix(0,nrow=2,ncol=100)
for (k in 1:100){
  tmp <- readRDS(paste0(k,'thetahat.RDS'))
  thetahat[,k] <- tmp
}

# calculate the estimated standard deviation based on the bootstrap samples for generating the SCB
beta1s <- beta.basismat%*%bhat1s
beta2s <- beta.basismat%*%bhat2s
sd1 <- apply(beta1s,1,sd)
sd2 <- apply(beta2s,1,sd)
beta.basismat <- readRDS('beta.basismat.RDS')
cvresult1 <- readRDS('cvresult1.RDS')
optid <- which(cvresult1[,7]==min(cvresult1[,7]))
tmpcv <- cvresult1[cvresult1[,7]==min(cvresult1[,7]),]
bhat1s <- readRDS('bhat1s.RDS')
bhat2s <- readRDS('bhat2s.RDS')
initial_est <- readRDS('initial_est.RDS')

# plot the estimated slope function for the daily maximum temperature obtained from CLoSE method
plot(c(30:335),beta.basismat%*%bhat1s[optid[1],],type='l',xlab='Month',ylab=TeX("$\\beta_{1,0.5}$"),cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=3,ylim=range(c(beta.basismat%*%bhat1s[optid[1],], beta.basismat%*%initial_est[1:ncol(bhat1s)+3])),col='red')
axis(1, at=c(30,61,91,122,152,183,213,243,273,304,335), lab=c("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),cex.axis=1.5)
# plot the estimated slope function associated the daily maximum temperature obtained from SQL method
lines(c(30:335),beta.basismat%*%initial_est[1:ncol(bhat1s)+3],col='blue',lty=2,lwd=3)

tmp <- bhat2s[optid[1],]
s <- sum(tmp!=0)
qtau <- (2*log(s))^0.5-(2*log(s))^(-0.5)*(log(-0.5*log(1-tau))+0.5*(log(log(s))+log(4*pi)))
# plot the estimated slope function associated with the daily minimum temperature obtained from CLoSE method
plot(c(30:335),beta.basismat%*%bhat2s[optid[1],],type='l',xlab='Month',ylab=TeX("$\\beta_{2,0.25}$"),cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=3,ylim=range(c(beta.basismat%*%bhat2s[optid[1],]+sd2*qtau, beta.basismat%*%bhat2s[optid[1],]-sd2*qtau, beta.basismat%*%initial_est[1:ncol(bhat1s)*2+3]))*1.25)
axis(1, at=c(30,61,91,122,152,183,213,243,273,304,335), lab=c("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),cex.axis=1.5)
index <- min(which(abs(beta.basismat%*%bhat2s[optid[1],])>=1e-1))
xx <- c(30:335)[c(index:306)]
yy0 <- (beta.basismat%*%bhat2s[optid[1],])[c(index:306)]
yy1 <- (beta.basismat%*%bhat2s[optid[1],]+sd2*qtau)[c(index:306)]
yy2 <- (beta.basismat%*%bhat2s[optid[1],]-sd2*qtau)[c(index:306)]
# generate SCB for the slope function associated with the daily minimum temperature
polygon(c(xx,rev(xx)),c(yy1,rev(yy0)),col='grey90',border='grey90')
polygon(c(xx,rev(xx)),c(yy2,rev(yy0)),col='grey90',border='grey90')
lines(c(30:335),beta.basismat%*%bhat2s[optid[1],],lwd=3,col='red')
# plot the estimated slope function associated the daily minimum temperature obtained from SQL method
lines(c(30:335),beta.basismat%*%initial_est[1:ncol(bhat1s)*2+3],col='blue',lty=2,lwd=3)

# Generate Figure 4(c) and 4(d) of the manuscript
setwd('~/r_functions/soybean_data_analysis/0.5-quantile/')
# estimations obtained from bootstrap samples fo
bhat2s <- bhat1s <- matrix(NA,nrow=55,ncol=100)
for (k in 1:100){
  tmp <- readRDS(paste0(k,'bHat1.RDS'))
  bhat1s[,k]<- tmp
  tmp <- readRDS(paste0(k,'bHat2.RDS'))
  bhat2s[,k]<- tmp
}
thetahat <- matrix(0,nrow=2,ncol=100)
for (k in 1:100){
  tmp <- readRDS(paste0(k,'thetahat.RDS'))
  thetahat[,k] <- tmp
}

# calculate the estimated standard deviation based on the bootstrap samples for generating the SCB
beta1s <- beta.basismat%*%bhat1s
beta2s <- beta.basismat%*%bhat2s
sd1 <- apply(beta1s,1,sd)
sd2 <- apply(beta2s,1,sd)
beta.basismat <- readRDS('beta.basismat.RDS')
cvresult1 <- readRDS('cvresult1.RDS')
optid <- which(cvresult1[,7]==min(cvresult1[,7]))
tmpcv <- cvresult1[cvresult1[,7]==min(cvresult1[,7]),]
bhat1s <- readRDS('bhat1s.RDS')
bhat2s <- readRDS('bhat2s.RDS')
initial_est <- readRDS('initial_est.RDS')

tmp <- bhat1s[optid[1],]
s <- sum(tmp!=0)
qtau <- (2*log(s))^0.5-(2*log(s))^(-0.5)*(log(-0.5*log(1-tau))+0.5*(log(log(s))+log(4*pi)))
# plot the estimated slope function for the daily maximum temperature obtained from CLoSE method
plot(c(30:335),beta.basismat%*%bhat1s[optid[1],],type='l',xlab='Month',ylab=TeX("$\\beta_{1,0.5}$"),cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=3,ylim=range(c(beta.basismat%*%bhat1s[optid[1],]+sd1*qtau, beta.basismat%*%bhat1s[optid[1],]-sd1*qtau, beta.basismat%*%initial_est[1:ncol(bhat1s)+3])),col='red')
axis(1, at=c(30,61,91,122,152,183,213,243,273,304,335), lab=c("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),cex.axis=1.5)
index <- min(which(abs(beta.basismat%*%bhat1s[optid[1],])>=1e-1))
xx <- c(30:335)[c(index:306)]
yy0 <- (beta.basismat%*%bhat1s[optid[1],])[c(index:306)]
yy1 <- (beta.basismat%*%bhat1s[optid[1],]+sd1*qtau)[c(index:306)]
yy2 <- (beta.basismat%*%bhat1s[optid[1],]-sd1*qtau)[c(index:306)]
# generate SCB for the slope function associated with the daily maximum temperature
polygon(c(xx,rev(xx)),c(yy1,rev(yy0)),col='grey90',border='grey90')
polygon(c(xx,rev(xx)),c(yy2,rev(yy0)),col='grey90',border='grey90')
lines(c(30:335),beta.basismat%*%bhat1s[optid[1],],lwd=3,col='red')
# plot the estimated slope function associated the daily maximum temperature obtained from SQL method
lines(c(30:335),beta.basismat%*%initial_est[1:ncol(bhat1s)+3],col='blue',lty=2,lwd=3)


# plot the estimated slope function associated with the daily minimum temperature obtained from CLoSE method
plot(c(30:335),beta.basismat%*%bhat2s[optid[1],],type='l',xlab='Month',ylab=TeX("$\\beta_{2,0.5}$"),cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=3,ylim=range(c(beta.basismat%*%bhat2s[optid[1],], beta.basismat%*%initial_est[1:ncol(bhat1s)*2+3]))*1.25)
axis(1, at=c(30,61,91,122,152,183,213,243,273,304,335), lab=c("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),cex.axis=1.5)

# plot the estimated slope function associated the daily minimum temperature obtained from SQL method
lines(c(30:335),beta.basismat%*%initial_est[1:ncol(bhat1s)*2+3],col='blue',lty=2,lwd=3)

# Generate Figure 4(e) and 4(f) of the manuscript
setwd('~/r_functions/soybean_data_analysis/0.75-quantile/')
# estimations obtained from bootstrap samples fo
bhat2s <- bhat1s <- matrix(NA,nrow=55,ncol=100)
for (k in 1:100){
  tmp <- readRDS(paste0(k,'bHat1.RDS'))
  bhat1s[,k]<- tmp
  tmp <- readRDS(paste0(k,'bHat2.RDS'))
  bhat2s[,k]<- tmp
}
thetahat <- matrix(0,nrow=2,ncol=100)
for (k in 1:100){
  tmp <- readRDS(paste0(k,'thetahat.RDS'))
  thetahat[,k] <- tmp
}

# calculate the estimated standard deviation based on the bootstrap samples for generating the SCB
beta1s <- beta.basismat%*%bhat1s
beta2s <- beta.basismat%*%bhat2s
sd1 <- apply(beta1s,1,sd)
sd2 <- apply(beta2s,1,sd)
beta.basismat <- readRDS('beta.basismat.RDS')
cvresult1 <- readRDS('cvresult1.RDS')
optid <- which(cvresult1[,7]==min(cvresult1[,7]))
tmpcv <- cvresult1[cvresult1[,7]==min(cvresult1[,7]),]
bhat1s <- readRDS('bhat1s.RDS')
bhat2s <- readRDS('bhat2s.RDS')
initial_est <- readRDS('initial_est.RDS')

tmp <- bhat1s[optid[1],]
s <- sum(tmp!=0)
qtau <- (2*log(s))^0.5-(2*log(s))^(-0.5)*(log(-0.5*log(1-tau))+0.5*(log(log(s))+log(4*pi)))
# plot the estimated slope function for the daily maximum temperature obtained from CLoSE method
plot(c(30:335),beta.basismat%*%bhat1s[optid[1],],type='l',xlab='Month',ylab=TeX("$\\beta_{1,0.75}$"),cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=3,ylim=range(c(beta.basismat%*%bhat1s[optid[1],]+sd1*qtau, beta.basismat%*%bhat1s[optid[1],]-sd1*qtau, beta.basismat%*%initial_est[1:ncol(bhat1s)+3])),col='red')
axis(1, at=c(30,61,91,122,152,183,213,243,273,304,335), lab=c("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),cex.axis=1.5)
index <- min(which(abs(beta.basismat%*%bhat1s[optid[1],])>=1e-1))
xx <- c(30:335)[c(index:306)]
yy0 <- (beta.basismat%*%bhat1s[optid[1],])[c(index:306)]
yy1 <- (beta.basismat%*%bhat1s[optid[1],]+sd1*qtau)[c(index:306)]
yy2 <- (beta.basismat%*%bhat1s[optid[1],]-sd1*qtau)[c(index:306)]
# generate SCB for the slope function associated with the daily maximum temperature
polygon(c(xx,rev(xx)),c(yy1,rev(yy0)),col='grey90',border='grey90')
polygon(c(xx,rev(xx)),c(yy2,rev(yy0)),col='grey90',border='grey90')
lines(c(30:335),beta.basismat%*%bhat1s[optid[1],],lwd=3,col='red')
# plot the estimated slope function associated the daily maximum temperature obtained from SQL method
lines(c(30:335),beta.basismat%*%initial_est[1:ncol(bhat1s)+3],col='blue',lty=2,lwd=3)


# plot the estimated slope function associated with the daily minimum temperature obtained from CLoSE method
plot(c(30:335),beta.basismat%*%bhat2s[optid[1],],type='l',xlab='Month',ylab=TeX("$\\beta_{2,0.75}$"),cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=3,ylim=range(c(beta.basismat%*%bhat2s[optid[1],], beta.basismat%*%initial_est[1:ncol(bhat1s)*2+3])),col='red')
axis(1, at=c(30,61,91,122,152,183,213,243,273,304,335), lab=c("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),cex.axis=1.5)
# plot the estimated slope function associated the daily minimum temperature obtained from SQL method
lines(c(30:335),beta.basismat%*%initial_est[1:ncol(bhat1s)*2+3],col='blue',lty=2,lwd=3)