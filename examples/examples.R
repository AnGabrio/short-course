#example corr
library(MASS)
library(R2jags)
n<-500
mu<-c(0,0)
r<-0.98
Sigma<-matrix(c(1,r,r,1),nrow = 2, ncol = 2)
outcome<-mvrnorm(n=n,mu=mu,Sigma = Sigma)
mean<-c(0,0)
precis<-diag(2)
omega<-diag(2)
datalist<-list("n","outcome","mean","precis","omega")
parm<-c("mu","mu.diff")

set.seed(1234)
model.corr<-jags(data = datalist,parameters.to.save = parm,n.chains = 2,n.iter = 20000,model.file = "model.corr.txt")
c(quantile(model.corr$BUGSoutput$sims.list$mu.diff,probs = 0.025),mean(model.corr$BUGSoutput$sims.list$mu.diff),quantile(model.corr$BUGSoutput$sims.list$mu.diff,probs = 0.975))

datalist<-list("n","outcome")
set.seed(1234)
model.ind<-jags(data = datalist,parameters.to.save = parm,n.chains = 2,n.iter = 20000,model.file = "model.ind.txt")
c(quantile(model.ind$BUGSoutput$sims.list$mu.diff,probs = 0.025),mean(model.ind$BUGSoutput$sims.list$mu.diff),quantile(model.ind$BUGSoutput$sims.list$mu.diff,probs = 0.975))

#compare estimates
#compute densities
dens_corr<-density(model.corr$BUGSoutput$sims.list$mu.diff, bw=0.02)
dens_ind<-density(model.ind$BUGSoutput$sims.list$mu.diff, bw=0.02)
plot(dens_ind,xlim=c(-0.4,0.4),ylim=c(0,20), xlab="treatment effect", ylab="Density", type="n",bty="n",axes=F, main = "")
axis(1)
axis(2)
abline(v=0,lty=2)
polygon(dens_corr,col="red",border = "black",density = 30,angle =135 )
polygon(dens_ind,col="blue",border = "black",density = 30,angle =45 )
legend("topright",legend = c("correlated","independent"),fill=c("red","blue"),border = c("black","black"),
       density = c(30,30),angle =c(135,45), bty = "n", cex = 1, x.intersp = 1, y.intersp = 1,seg.len = 1)
text(-0.3,20,"mean  (95% CI)")
text(-0.28,19,paste(round(mean(dens_corr$x),2)," (",round(quantile(dens_corr$x, probs = 0.025),2),";",round(quantile(dens_corr$x, probs = 0.975),2),")",sep = ""),col="red",cex = 1)
text(-0.28,18,paste(round(mean(dens_ind$x),2)," (",round(quantile(dens_ind$x, probs = 0.025),2),";",round(quantile(dens_ind$x, probs = 0.975),2),")",sep = ""),col="blue",cex = 1)

#saveRDS(dens_corr,"dens_corr.rds")
#saveRDS(dens_ind,"dens_ind.rds")


#example skewness
library(fitdistrplus)
n<-50
mu<-500
sigma<-500
shape<-mu^2/(sigma^2)
rate<-mu/(sigma^2)
set.seed(12345)
outcome<-rgamma(n=n,shape = shape, rate = rate)
fit1.gamma <- fitdist(outcome, distr = c("gamma"), method = "mme")
fit2.lnorm <- fitdist(outcome, distr = c("lnorm"), method = "mme")
fit3.norm <- fitdist(outcome, distr = c("norm"), method = "mme")
plot.legend <- c("Gamma", "LogNormal","Normal")
denscomp(list(fit1.gamma,fit2.lnorm,fit3.norm),legendtext = plot.legend,ylim = c(0,0.002),xlim = c(-500,2000),
         xlab = "Costs (£)", fitlty=c(1,1,1), fitcol = c("blue","darkgreen","red"), lwd=c(2,2,2), cex=1,
         main = "")


#example ones
library(fitdistrplus)
n<-100
mu<-0.8
sigma<-0.08
shape1<-mu*((mu*(1-mu)/sigma^2)-1)
shape2<-(1-mu)*((mu*(1-mu)/sigma^2)-1)
set.seed(12345)
outcome<-rbeta(n=n,shape1 = shape1,shape2 = shape2)
outcome<-c(outcome,rep(1,30))
fit1.beta <- fitdist(outcome, distr = c("beta"), method = "mme")
fit3.norm <- fitdist(outcome, distr = c("norm"), method = "mme")
plot.legend <- c("Beta","Normal")
denscomp(list(fit1.beta,fit3.norm),legendtext = plot.legend,ylim = c(0,6),xlim = c(0.5,1.1),
         xlab = "Utilities", fitlty=c(1,1), fitcol = c("blue","red"), lwd=c(2,2), cex=1,
         main = "")



#missing
library(MASS)
logistic <- function(x) exp(x)/(1 + exp(x))
x <- -0.85
y <- -0.85                           
R <- matrix(numeric(2*2), nrow=2)   
diag(R) <- 1                       
R[upper.tri(R)] <- c(x) 
R[lower.tri(R)] <- c(y) 
eigen(R)$values  
vars  <- c(0.1^2, 500^2)
Sigma <- diag(sqrt(vars)) %*% R %*% diag(sqrt(vars))
n<-1000
mu<-c(0.5,1500)
set.seed(80122)
outcome<-mvrnorm(n=n,mu=mu,Sigma = Sigma)
m.mcar<-rbinom(n=n,size = 1,prob = 0.5)
m.mar<-rbinom(n=n,size = 1,prob = logistic(-5+0.003*outcome[,2]))
m.mnar<-rbinom(n=n,size = 1,prob = logistic(-4+8*outcome[,1]))


#MCAR
par(mfrow=c(1,2))
plot(outcome,main = "",xlab = "utilities",ylab = "costs (£)",type = "n",ylim = c(100,3000),xlim = c(0,1))
points(outcome,pch=1,col="black")
points(outcome[m.mcar==0,1],outcome[m.mcar==0,2],pch=1,col="blue")
points(outcome[m.mcar==1,1],outcome[m.mcar==1,2],pch=1,col="red")
legend("topright",legend = c("observed","missing"),col=c("blue","red"),pch=c(1,1),bty="n",seg.len = 1,
       x.intersp = 0.5,y.intersp = 1, lwd=c(2,2), lty=c(0,0))

plot(density(outcome[,1],bw=0.1),main = "",xlab = "utilities",ylab = "density",type = "n",axes=F,ylim=c(0,4),xlim=c(0,1))
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1))
axis(2,at=c(0,1,2,3,4),labels = c(0,1,2,3,4))
lines(density(outcome[m.mcar==0,1],bw=0.1),col="blue",lty=1,lwd=2)
polygon(density(outcome[m.mcar==0,1],bw=0.1),col="blue",border = "black",density = 30,angle =135 )
lines(density(outcome[m.mcar==1,1],bw=0.1),col="red",lty=1,lwd=2)
polygon(density(outcome[m.mcar==1,1],bw=0.1),col="red",border = "black",density = 30,angle =45 )
lines(density(outcome[,1],bw=0.1),col="black",lty=2,lwd=2)
abline(v=mean(outcome[,1]),lty=2,lwd=2,col="black")
abline(v=mean(outcome[m.mcar==0,1]),lty=1,lwd=1,col="blue")
abline(v=mean(outcome[m.mcar==1,1]),lty=1,lwd=1,col="red")
text(0.82,4,"mean")
text(0.832,3.8,paste("all: ",round(mean(outcome[,1]),2),sep = " "),col="black",cex = 1)
text(0.805,3.6,paste("observed: ",round(mean(outcome[m.mcar==0,1]),2),sep = " "),col="blue",cex = 1)
text(0.81,3.4,paste("missing: ",round(mean(outcome[m.mcar==1,1]),2),sep = " "),col="red",cex = 1)


#MAR
plot(outcome,main = "",xlab = "utilities",ylab = "costs (£)",type = "n",ylim = c(100,3000),xlim = c(0,1))
points(outcome,pch=1,col="black")
points(outcome[m.mar==0,1],outcome[m.mar==0,2],pch=1,col="blue")
points(outcome[m.mar==1,1],outcome[m.mar==1,2],pch=1,col="red")
legend("topright",legend = c("observed","missing"),col=c("blue","red"),pch=c(1,1),bty="n",seg.len = 1,
       x.intersp = 0.5,y.intersp = 1, lwd=c(2,2), lty=c(0,0))

plot(density(outcome[,1],bw=0.1),main = "",xlab = "utilities",ylab = "density",type = "n",axes=F,ylim=c(0,4),xlim=c(0,1))
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1))
axis(2,at=c(0,1,2,3,4),labels = c(0,1,2,3,4))
lines(density(outcome[m.mar==0,1],bw=0.1),col="blue",lty=1,lwd=2)
polygon(density(outcome[m.mar==0,1],bw=0.1),col="blue",border = "black",density = 30,angle =135 )
lines(density(outcome[m.mar==1,1],bw=0.1),col="red",lty=1,lwd=2)
polygon(density(outcome[m.mar==1,1],bw=0.1),col="red",border = "black",density = 30,angle =45 )
lines(density(outcome[,1],bw=0.1),col="black",lty=2,lwd=2)
abline(v=mean(outcome[,1]),lty=2,lwd=2,col="black")
abline(v=mean(outcome[m.mar==0,1]),lty=1,lwd=1,col="blue")
abline(v=mean(outcome[m.mar==1,1]),lty=1,lwd=1,col="red")
text(0.82,4,"mean")
text(0.832,3.8,paste("all: ",round(mean(outcome[,1]),2),sep = " "),col="black",cex = 1)
text(0.805,3.6,paste("observed: ",round(mean(outcome[m.mar==0,1]),2),sep = " "),col="blue",cex = 1)
text(0.81,3.4,paste("missing: ",round(mean(outcome[m.mar==1,1]),2),sep = " "),col="red",cex = 1)


#MNAR
plot(outcome,main = "",xlab = "utilities",ylab = "costs (£)",type = "n",ylim = c(100,3000),xlim = c(0,1))
points(outcome,pch=1,col="black")
points(outcome[m.mnar==0,1],outcome[m.mnar==0,2],pch=1,col="blue")
points(outcome[m.mnar==1,1],outcome[m.mnar==1,2],pch=1,col="red")
legend("topright",legend = c("observed","missing"),col=c("blue","red"),pch=c(1,1),bty="n",seg.len = 1,
       x.intersp = 0.5,y.intersp = 1, lwd=c(2,2), lty=c(0,0))

plot(density(outcome[,1],bw=0.1),main = "",xlab = "utilities",ylab = "density",type = "n",axes=F,ylim=c(0,4),xlim=c(0,1))
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1))
axis(2,at=c(0,1,2,3,4),labels = c(0,1,2,3,4))
lines(density(outcome[m.mnar==0,1],bw=0.1),col="blue",lty=1,lwd=2)
polygon(density(outcome[m.mnar==0,1],bw=0.1),col="blue",border = "black",density = 30,angle =135 )
lines(density(outcome[m.mnar==1,1],bw=0.1),col="red",lty=1,lwd=2)
polygon(density(outcome[m.mnar==1,1],bw=0.1),col="red",border = "black",density = 30,angle =45 )
lines(density(outcome[,1],bw=0.1),col="black",lty=2,lwd=2)
abline(v=mean(outcome[,1]),lty=2,lwd=2,col="black")
abline(v=mean(outcome[m.mnar==0,1]),lty=1,lwd=1,col="blue")
abline(v=mean(outcome[m.mnar==1,1]),lty=1,lwd=1,col="red")
text(0.82,4,"mean")
text(0.832,3.8,paste("all: ",round(mean(outcome[,1]),2),sep = " "),col="black",cex = 1)
text(0.805,3.6,paste("observed: ",round(mean(outcome[m.mnar==0,1]),2),sep = " "),col="blue",cex = 1)
text(0.81,3.4,paste("missing: ",round(mean(outcome[m.mnar==1,1]),2),sep = " "),col="red",cex = 1)

















