# Calculation of L and ARLs for different values of rho

library(matrixStats)
library(matrixcalc)


rm(list=ls()) # clear workspace
cat("\014") # Clear console
alpha<-0.0027
m=100 # number of streams
n=1   # sample size
s<-1 #standard deviaion of each stream
al<-c(1, floor(0.25*m), floor(0.5*m), floor(0.75*m))
alpha2opt<-matrix(rep(0, 128),64,2) #optimal set of alpha2
j<-0
for (rho in seq(0, 0.75, by=0.25)){
for (r in al){ #number of shifted streams
  for (d in seq(0.5, 2, by=0.5)){ #amount of shift
j<-j+1
sa<-s*sqrt(rho/(1-rho)) #common standard deviation
st<-sqrt(s^2+sa^2) # total standard deviation
sigma <- (s^2)*diag(m)+(sa^2)*matrix(1,m,m) # covariance matrix
delta<-c(rep(d, r),rep(0,m-r)) #mu1-mu0
lambda<-n*t(delta)%*%solve(sigma) %*%(delta) # non- centrality parameter

# computing Power for the proposed method
PowP<-c(rep(0,1000)) #power for the proposed method
alpha2<-c(seq((10^(-3)), 1, (10^(-3))))
alpha1<-c(rep(0,1000))
a2<-0
for (i in 1:1000){
  
  f2<-function(x){exp(-x^2/2)*((pnorm(st*qnorm(1-alpha2[i]/2)/s-x*sa/s)-pnorm(-st*qnorm(1-alpha2[i]/2)/s-x*sa/s)))^m/((2*pi)^0.5)}
  I2<-integrate(f2, lower = -Inf, upper=Inf)
  a2<-I2[[1]]
  
  alpha1[i]<-alpha/(1-a2) # corresponding alpha1
  if (alpha1[i]<1){
    q=qchisq(1-alpha1[i], df=m, ncp = 0)
    PowS1<-(1-pchisq(q, df=m, ncp = lambda)) #power of first step
    f3<-function(x)((pnorm(st*qnorm(1-alpha2[i]/2)/s-x*sa/s-sqrt(n)*d/s)-pnorm(-st*qnorm(1-alpha2[i]/2)/s-x*sa/s-sqrt(n)*d/s)))^(r)
    f4<-function(x){exp(-x^2/2)*f3(x)*((pnorm(st*qnorm(1-alpha2[i]/2)/s-x*sa/s)-pnorm(-st*qnorm(1-alpha2[i]/2)/s-x*sa/s)))^(m-r)/((2*pi)^0.5)}
    I3<-integrate(f4, lower = -Inf, upper=Inf)
    PowS2<-(1-I3[[1]]) #power of second step
    PowP[i]<-(PowS1*PowS2)
  }
}
alpha2opt[j,1]<-alpha2[which.max(PowP)]
alpha2opt[j,2]<-rho
  }
}
}
lm1<-lm(formula = alpha2opt[, 1] ~ poly(alpha2opt[, 2], 3, raw = TRUE)) #polynomial regression fitting
plot(alpha2opt[, 2],alpha2opt[, 1], pch = 16, cex = .9,xlab="rho", ylab="alpha_2") #plot of data
f <- function(x) coef(lm1)[1] + coef(lm1)[2]*x + coef(lm1)[3]*x^2+coef(lm1)[4]*x^3
lines(alpha2,f(alpha2), col="red") 
print(c(coef(lm1)[1],coef(lm1)[2],coef(lm1)[3],coef(lm1)[4]))
print(summary(lm1)$ r.squared)
print(summary(lm1)$adj.r.squared)

