library(MASS)
library(matrixStats)
library(matrixcalc)
library(expm)
library(tictoc)
#library(xlsx)
library(psych)
library(ltm)
tic() # Start the clock!

# Pararmeters Setting -----------------------------------------------------

rm(list=ls()) # clear workspace
cat("\014") # Clear console
load("VDPdata.csv")
x0<-as.matrix(VDPdata)
m=NCOL(x0) # number of histoical (phase I) samples
p=NROW(x0)
alpha=0.005 #type I error probability
del<-alpha/2 #significance level in phase I
k=10; #initial subset for MDP (Phase I)
h=floor(m/2)+1; #subsets size for MDP (Phase I)

  # Phase I analysis (robust parameter estimation) ---------------------------------
  
  MinSet<-matrix(rep(0, k*(h+1)),k,h+1)
  d2raw<-c(rep(0,m))
  repp=0
  for (i in 1:k){
    y<-matrix(x0[,i:(i+h-1)],p,h) 
    Sety<-i:(i+h-1)
    T1<-rowMeans(y)
    S1<-matrix(rep(0, p*p),p,p)
    
    for (j in 1:h){
      S1<-S1+(y[,j]-T1)%*%t(y[,j]-T1)/h} # or cov()
    D1<-diag(diag(S1))
    d1<-matrix(rep(0, m*2),m,2)
    d2<-matrix(rep(0, m*2),m,2)
    sgn1=0;
    for (j in 1:m){
      d1[j,1]<-t(x0[,j]-T1)%*%solve(D1) %*%(x0[,j]-T1)
      d1[j,2]<-j
    }
    
    d2<-dfOrder(d1,1)
    a<-d2[1:h,2]        #a<-dfOrder(d2[1:h,],2)[,2]
    while (det(D1)!=0 & sgn1==0){
      repp<-repp+1
      if (identical(sort(Sety),sort(a))){ #if (isTRUE(all.equal(d2[1:h,2],a))){
        sgn1=1}
      else {for (j in 1:h){y[,j]<-x0[,d2[j,2]]
      Sety[j]<-d2[j,2]}
        T1<-rowMeans(y)
        S1<-matrix(rep(0, p*p),p,p)
        for (j in 1:h){
          S1<-S1+(y[,j]-T1)%*%t(y[,j]-T1)/h}
        D1<-diag(diag(S1)) 
        for (j in 1:m){
          d1[j,1]<-t(x0[,j]-T1)%*%solve(D1)%*%(x0[,j]-T1)
          d1[j,2]<-j
        }
        d2<-dfOrder(d1,1)
        a<-d2[1:h,2]
      } 
    }
    d<-d2
    
    MinSet[i,h+1]<-det(D1)
    MinSet[i,1:h]<-d[1:h,2]
  }
  MinSet<-dfOrder(MinSet,h+1)
  for (i in 1:h){
    y[,i]<-x0[,MinSet[1,i]]
  }
  muhat<-rowMeans(y) # mu vector of MDP
  Sigmaraw<-matrix(rep(0, p*p),p,p) 
  for (j in 1:h){
    Sigmaraw<-Sigmaraw+(y[,j]-muhat)%*%t(y[,j]-muhat)/h} #covariance matrix of MDP
  Rraw<-solve(sqrtm(diag(diag(Sigmaraw))))%*%Sigmaraw%*%solve(sqrtm(diag(diag(Sigmaraw))))
  trR2mdp<-matrix.trace(Rraw%^% 2)-(p^2)/h
  trR3mdp<-matrix.trace(Rraw%^%3)-3*p*trR2mdp/h-p^3/h^2
  for (i in 1:m){
    d2raw[i]<-t(x0[,i]-muhat)%*%solve(diag(diag(Sigmaraw)))%*%(x0[,i]-muhat)  
  }
  
  c<-median(d2raw)/p
  Dhatmdp<-c*(diag(diag(Sigmaraw)))
  chatpn<-1+2*(p/(m*trR2mdp^0.5))
  w<-c(rep(0,m)) #weight
  TW<-c(rep(0,p))
  SW<-matrix(rep(0, p*p),p,p)
  d2mdp<-c(rep(0,m))
  for (i in 1:m){
    d2mdp[i]<-t(x0[,i]-muhat)%*%solve(Dhatmdp)%*%(x0[,i]-muhat) 
    if ((d2mdp[i]-p)/((2*chatpn*trR2mdp)^0.5)-(8*trR3mdp)*(qnorm(1-del)^2-1)/(6*((2*trR2mdp)^1.5))<=qnorm(1-del)){w[i]=1}
    TW<-TW+(w[i]*x0[,i])
  }
  mutilde<-TW/sum(w)
  for (i in 1:m){
    SW<-SW+w[i]*(x0[,i]-mutilde)%*%t(x0[,i]-mutilde)/(sum(w)-1)}
  D0tilde<-diag(diag(SW))
  Rtilde<-solve(sqrtm(D0tilde))%*%SW%*%solve(sqrtm(D0tilde))
  trR2w<-matrix.trace(Rtilde%^% 2)-(p^2)/sum(w)
  trR3w<-matrix.trace(Rtilde%^%3)-3*p*trR2w/sum(w)-p^3/sum(w)^2
  d2muD0hat<-c(rep(0,m))
  for (i in 1:m){
    d2muD0hat[i]<-t(x0[,i]-mutilde)%*%solve(D0tilde)%*%(x0[,i]-mutilde)}
  d2ref<-c(rep(0,m)) #Refined distance
  dd2ref<-d2ref
  ctildepn<-1+2*(p/(sum(w)*trR2w^0.5))
  w<-c(rep(0,m)) 
  TW<-c(rep(0,p))
  SW<-matrix(rep(0, p*p),p,p)
  for (i in 1:m){
    d2ref[i]<-d2muD0hat[i]/(1+(dnorm(qnorm(1-del))*((2*trR2w)^0.5)/(p*(1-del))))
    dd2ref[i]<-(d2ref[i]-p)/((2*ctildepn*trR2w)^0.5)-(8*trR3w)*(qnorm(1-alpha)^2-1)/(6*((2*trR2w)^1.5))
     if ((d2ref[i]-p)/((2*ctildepn*trR2w)^0.5)-(8*trR3w)*(qnorm(1-alpha)^2-1)/(6*((2*trR2w)^1.5))<=qnorm(1-alpha)){w[i]=1}
    TW<-TW+(w[i]*x0[,i])
  }
  mu0hat<-TW/sum(w)
  for (i in 1:m){
    SW<-SW+w[i]*(x0[,i]-mu0hat)%*%t(x0[,i]-mu0hat)/(sum(w)-1)}
  Dhat<-diag(diag(SW))
  Rho<-solve(sqrtm(Dhat))%*%SW%*%solve(sqrtm(Dhat))
  trRho2<-matrix.trace(Rho%^% 2)-(p^2)/sum(w)
  ctilde2pn<-1+2*(p/(sum(w)*trRho2^0.5))
  UCL<-p+((2*ctildepn*trR2w)^0.5)*(qnorm(1-alpha)+(8*trR3w)*(qnorm(1-alpha)^2-1)/(6*((2*trR2w)^1.5)))
  #Dhat1<-diag(diag(sigma0))
  #Rho<-solve(sqrtm(Dhat1))%*%sigma0%*%solve(sqrtm(Dhat1))
  #trRho2<-matrix.trace(Rho%^% 2)
  trR3<-matrix.trace(Rho%^%3)-3*p*trRho2/sum(w)-p^3/sum(w)^2
  write.xlsx(dd2ref, "3.xlsx")
  getwd()