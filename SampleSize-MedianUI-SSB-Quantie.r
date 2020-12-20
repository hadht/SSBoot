################################################################################
##   SSB Procedure to Compute Sample Size
##
## Author: Thiago Rezende
## Date: 12/20/2020
## All rights reserved
################################################################################
setwd("C://Users//Windows 8//Documents//Documentos//DocUFMG//ESAMP3//PosterEsamp-2017//Paper//AplicaçãoMedianaUI")

#functions:
#Compute BCA function:
bca <- function(theta, conf.level = .95){
  low <- (1 - conf.level)/2
  high <- 1 - low
  sims <- length(theta)
  z.inv <- length(theta[theta < mean(theta)])/sims
  z <- qnorm(z.inv)
  U <- (sims - 1) * (mean(theta, na.rm=TRUE) - theta)
  top <- sum(U^3)
  under <- 6 * (sum(U^2))^{3/2}
  a <- top / under
  lower.inv <-  pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
  lower <- quantile(theta, lower.inv, names=FALSE)
  upper.inv <-  pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper <- quantile(theta, upper.inv, names=FALSE)
  return(c(lower, upper))
}
# Standard Sample Size (S)
SSS<-function(m=NULL,d=NULL,alfa,q,Xp,N=NULL){
# m sample size of interest
# d error margin
# alfa siginificance level
# q order of quantile of interest to estimate
# Xp the pilot sample size
z=qnorm(1-alfa/2)
varr<-q*(1-q)
np<-length(Xp)
k<-round(sqrt(np))
#Compute the q-th sample quantile
sxnp<-sort(Xp)
r<-floor(np*q)+1
xnq<-sxnp[r]
#xr<-quantile(Xp,probs=c(q))
lsup<-sxnp[r+k]
linf<-sxnp[r-k]
if((r+k)>np){lsup<-sxnp[np]}
if((r-k)<1){linf<-sxnp[1]}
AA<-(lsup-linf)
smk<-(np*(AA))/(2*k)
if(is.null(m)){
outmd<- ceiling(((z^2)*varr*(smk^(2))/(d^2)))    #Erro aqui!!
}else{
if(is.null(N)){
outmd<-sqrt((z^2)*varr*(smk^(2))/m)
}else{
outmd<-sqrt((z^2)*((1-(m/N))*varr)*(smk^(2))/m)
}
}
return(outmd)
}


set.seed(100)
#set.seed(300)
#### Initial values
N=1000
sig=87
B=20000        # # Boot resamples
alpha=0.05
#Quantile order
q<-0.50
###############################################
nvar=c(20,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,150,300)
na=18
###############################################

### Declare Variables
mediaam=numeric(na)
damboot=numeric(na)
EAAS=numeric(na)

################################################################################
##
## UI DATA
##
################################################################################
set.seed(10)
n=65
x=numeric(n)
for(i in 1:n){
aux1=sample(1:4,replace=F)[1]
if(aux1==1){u1=runif(1,0,0.25);x[i]=(u1*164/0.25)}
if(aux1==2){u2=runif(1,0.25,0.5);x[i]=(u2*224/0.5)}
if(aux1==3){u3=runif(1,0.5,0.75);x[i]=(u3*286/0.75)}
if(aux1==4){u3=runif(1,0.75,1.0);x[i]=(u3*350/1)}
}
#xna=rnorm(nvar[1],mu,sig)
x
summary(x)
quantile(x,probs=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
hist(x)
y=x
################################################################################

### Main Program
for(i in 1:na){ # Loop amostrq
mediaam[i]= (median(y))
 xnamedboot=numeric(B)
for(j in 1:B){  # Loop bootstrap

ii=sample(1:length(y),nvar[i],replace=T)
auxbooty=y[ii]

xnamedboot[j]=(median(auxbooty))

}
## Compute do Error margin
results <- bca(xnamedboot, conf.level =1-alpha)
linfboot<- results[1]
lsupboot<- results[2]
Amplitude=lsupboot-linfboot
medboot=median(xnamedboot)
damboot[i]=(Amplitude/2)*sqrt(1-(nvar[i]/N))
z=qnorm(1-alpha/2)
s2r=var(y)*(pi/2)
EAAS[i]=SSS(nvar[i],d=NULL,alpha,q,y,N)

}

###############################################################################
## RESULTS
###############################################################################
nvar
damboot
EAAS
mean((damboot-EAAS)^2)
minn=min(c(min(damboot),min(EAAS)))
mann=max(c(max(damboot)+0.02,max(EAAS)+0.02))

#par(mfrow=c(1,1))
plot(nvar,damboot,ylim=c(minn,mann),xlab="n",type='p',ylab="d",pch=c(2))
lines(c(0,65),c(29.6,29.6),type='l',col="red",lty=c(2))
lines(c(65,65),c(0,29.6),type='l',col="red",lty=c(2))
par(new=T)
plot(nvar,EAAS,ylim=c(minn,mann),type='p',xlab="n",ylab="d",pch=c(1))
lines(c(0,57),c(29.6,29.6),type='l',col="red",lty=c(1))
lines(c(57,57),c(0,29.6),type='l',col="red",lty=c(1))
#legend(300,0.08,c("Boot","Usual Formula (Ap.)"),pch=c(2,1))
legend(200,40,c("SSB","S"),pch=c(2,1))
title("Median")








