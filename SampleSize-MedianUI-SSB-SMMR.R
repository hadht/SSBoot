################################################################################
##   SSB Procedure to Compute Sample Size
##
## Author: Thiago Rezende
## Date: 11/08/2022
## All rights reserved
################################################################################

#Libraries:
library("bootstrap")

#Function:
#Compute BCA bootstrap confidence interval:
#Begin Function
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
#End Function

################################################################################
#Function:
#Compute The Exact (Standard) Sample Size (S):
# m sample size of interest
# d error margin
# alfa siginificance level
# q order of quantile of interest to estimate
# Xp the pilot sample size
################################################################################
#Begin Function
SSS<-function(m=NULL,d=NULL,alfa,q,Xp,N=NULL){
 z=qnorm(1-alfa/2)
 varr<-q*(1-q)
 np<-length(Xp)
 k<-round(sqrt(np))
 #Compute the q-th sample quantile
 sxnp<-sort(Xp)
 r<-floor(np*q)+1
 xnq<-sxnp[r]
 lsup<-sxnp[r+k]
 linf<-sxnp[r-k]
 if((r+k)>np){lsup<-sxnp[np]}
 if((r-k)<1){linf<-sxnp[1]}
 AA<-(lsup-linf)
 smk<-(np*(AA))/(2*k)
 if(is.null(m)){
  outmd<- ceiling(((z^2)*varr*(smk^(2))/(d^2)))
 }else{
  if(is.null(N)){
 outmd<-sqrt((z^2)*varr*(smk^(2))/m)
 }else{
  outmd<-sqrt((z^2)*((1-(m/N))*varr)*(smk^(2))/m)
 }
 }
 return(outmd)
}
#End Function

################################################################################
#Function:
#The Sample Size using the procedure SSB:
################################################################################
#Begin Function
SSB<-function(y=y,m=nvar,errormarginset=NULL,B=20000,alpha=0.05,q=0.50,N=1000000){
 na<-length(m)
 ### Declare Variables
 mediaam=numeric(na)
 damboot=numeric(na)
 ### Loop sample size
 for(i in 1:na){ # Loop
   mediaam[i]= quantile(y,probs=c(0.5))
   xnamedboot=numeric(B)
   for(j in 1:B){  # Loop bootstrap
    ii=sample(1:length(y),nvar[i],replace=T)
    auxbooty=y[ii]
    xnamedboot[j]=quantile(auxbooty,probs=c(q))
   }
   ## Compute do Error margin
   #BCA bootstrap confidence interval
   results <- bca(xnamedboot, conf.level =1-alpha)
   linfboot<- results[1]
   lsupboot<- results[2]
   #Width
   Amplitude=lsupboot-linfboot
   medboot=median(xnamedboot)
   # Error margin
   damboot[i]=(Amplitude/2)*sqrt(1-(nvar[i]/N))
 }
 if(is.null(errormarginset)){
    out<-list(nvar,damboot)
    names(out)[1]="Sample size grid:"
    names(out)[2]="Error margin E grid:"
 }else{
    interplinear_n = approx(damboot, nvar, xout=errormarginset)
    interplinear_n$y=round(interplinear_n$y)
    out<-list(nvar,damboot,interplinear_n$y,interplinear_n$x)
    names(out)[1]="Sample size grid:"
    names(out)[2]="Error margin E grid:"
    names(out)[3]="Error margin desired E = "
    names(out)[4]="Sample size n = "
 }
 return(out)
}
#End Function
################################################################################
## Main Program: Results
################################################################################
#### Initial values
#set seed:
set.seed(100)
# Population Size:
N=1000
# Number of Bootstrap Samples:
B=5000
# Error Type I:
alpha=0.05
#Quantile to be estimated:
q<-0.50
#Load UI Data (pilot sample):
#n=65 obs.
y=scan("uidata.txt")

################################################################################
# Compute the standard or exact margin of error for several sample sizes
################################################################################
# Compute error margin with 18 sample sizes (m) - a grid:
# 18
nvar=c(20,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,150,300)
na=length(nvar)
EAAS=numeric(na)
for(i in 1:na){ # Loop amostra
# Call SSS function
EAAS[i]=SSS(nvar[i],d=NULL,alfa=alpha,q=q,y,N=N)
}
# sample sizes
nvar
# Error margin with the procedure S:
EAAS

################################################################################
#Compute The Sample Size using the procedure SSB:
################################################################################
# Call SSB function
# Error margin with the procedure SSB:
errormarginset=c(29)
damboot2=SSB(y=y,m=nvar,errormarginset,B=B,alpha=alpha,q=q,N=N)
damboot2
# sample sizes
nvar=damboot2[[1]]
# Error margin with the procedure SSB:
damboot=damboot2[[2]]

################################################################################
#Graphs:
################################################################################
#Median: #########
minn=min(c(min(damboot),min(EAAS)))
mann=max(c(max(damboot)+0.02,max(EAAS)+0.02))
plot(nvar,damboot,ylim=c(minn,mann),xlab="n",type='o',ylab="d",pch=c(2))
lines(c(0,43),c(29.6,29.6),type='l',col="red",lty=c(2))
lines(c(43,43),c(0,29.6),type='l',col="red",lty=c(2))
par(new=T)
plot(nvar,EAAS,ylim=c(minn,mann),type='o',xlab="n",ylab="d",pch=c(1))
lines(c(0,37),c(29.6,29.6),type='l',col="red",lty=c(1))
lines(c(37,37),c(0,29.6),type='l',col="red",lty=c(1))
legend(200,40,c("SSB","S"),pch=c(2,1))
title("Median")

#P_75: ###########
minn=min(c(min(damboot),min(EAAS)))
mann=max(c(max(damboot)+0.02,max(EAAS)+0.02))
plot(nvar,damboot,ylim=c(minn,mann),xlab="n",type='o',ylab="d",pch=c(2))
lines(c(0,106),c(29.6,29.6),type='l',col="red",lty=c(2))
lines(c(106,106),c(0,29.6),type='l',col="red",lty=c(2))
par(new=T)
plot(nvar,EAAS,ylim=c(minn,mann),type='o',xlab="n",ylab="d",pch=c(1))
lines(c(0,85),c(29.6,29.6),type='l',col="red",lty=c(1))
lines(c(85,85),c(0,29.6),type='l',col="red",lty=c(1))
legend(200,40,c("SSB","S"),pch=c(2,1))
title("75% Quantile")






