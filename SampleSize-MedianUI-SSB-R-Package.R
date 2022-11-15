################################################################################
##   SSB Procedure to Compute Sample Size
##
## Author: Thiago Rezende
## Date: 2022/11/08
## All rights reserved
################################################################################

#Libraries:
library("SSB")
library("bootstrap")

################################################################################
# Compute the standard or exact margin of error for several sample sizes
################################################################################
# Compute error margin with 18 sample sizes (m) - a grid:
# 18
################################################################################
## Main Program: Results
################################################################################
#### Initial values
#set seed:
set.seed(100)
# Population Size:
N=1000
# Number of Bootstrap Samples:
B=10000
# Error Type I:
alpha=0.05
#Quantile to be estimated:
q<-0.50
#Load UI Data (pilot sample):
#n=65 obs.
data(UI)
y=UI
nvar<-c(20,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,150,300)
na=length(nvar)

################################################################################
#Compute The Sample Size using the procedure SSB:
################################################################################
# Call SSB function
# Error margin with the procedure SSB:
errormarginset=c(29)
damboot2=SSB(y=y,m=nvar,errormarginset,B=B,alpha=alpha,q=q,N=N)
damboot2
# sample sizes
nvar1=damboot2[[1]]
# Error margin with the procedure SSB:
damboot1=damboot2[[2]]

################################################################################
#Compute The Sample Size using the procedure SSB:
################################################################################
# Call SSB function
# Error margin with the procedure SSB:
errormarginset=c(29)
q<-0.75
damboot2=SSB(y=y,m=nvar,errormarginset,B=B,alpha=alpha,q=q,N=N)
damboot2
# sample sizes
nvar2=damboot2[[1]]
# Error margin with the procedure SSB:
damboot2=damboot2[[2]]

################################################################################
#Graphs:
################################################################################
#Median: #########
minn=min(damboot)
mann=max(damboot)
plot(nvar1,damboot1,ylim=c(minn,mann),xlab="n",type='o',ylab="d",pch=c(2))
lines(c(0,43),c(29.6,29.6),type='l',col="red",lty=c(2))
lines(c(43,43),c(0,29.6),type='l',col="red",lty=c(2))
legend(200,40,c("SSB"),pch=c(2))
title("Median")

#P_75: ###########
minn=min(damboot2)
mann=max(damboot2)
plot(nvar2,damboot2,ylim=c(minn,mann),xlab="n",type='o',ylab="d",pch=c(2))
lines(c(0,106),c(29.6,29.6),type='l',col="red",lty=c(2))
lines(c(106,106),c(0,29.6),type='l',col="red",lty=c(2))
legend(200,40,c("SSB"),pch=c(2))
title("75 Quantile")







