#without sorting data
# install packages "MASS", "matrixcalc", expm", "matrixStats", "psych" and"tictoc"
# Code for MDP estimator
# install packages "MASS" and "matrixStats" and "tictoc"
library(MASS)
library(matrixStats)
library(matrixcalc)
library(expm)
library(tictoc)
#library(xlsx)
library(ltm)
library(psych)

rm(list=ls()) # clear workspace
cat("\014") # Clear console
load("Example1.xlsx")
load("SelectedSemi.xlsx")
load("Example2.xlsx")
#Phase I and outlier detection
SV<-(as.matrix(SelectedSemi))
#SV<-seq(1,250,by=1)
x0<-Example2
x0<-t(x0)
rr<-floor(0.3*1567)
#

for (i in 1:400){
  f<-ecdf(x0[i, 1:rr])
  for (j in 1:1567){
  x0[i,j]<-qnorm(abs(f(x0[i,j])-10^(-5)))}}
p=250 # number of variables
m<-rr # number of histoical (phase I) samples
x1<-x0
x0<-x0[,1:rr]
x00<-matrix(rep(0, p*m),p,m)
for (i in 1:p){
  x00[i,]<-x0[SV[i],]
}
x0<-x00



#print(x0[,1])
#mydata <- read.table("Example.csv", header=TRUE, sep=",", row.names="id")

k=5; #initial subset
h=floor(m/2)+1; #subsets size
alpha=0.005
del<-alpha/2
MinSet<-matrix(rep(0, k*(h+1)),k,h+1)
d2raw<-c(rep(0,m))

for (i in 1:k){
  y<-matrix(x0[,i:(i+h-1)],p,h) 
  Sety<-i:(i+h-1)
  T1<-rowMeans(y)
  S1<-matrix(rep(0, p*p),p,p)
  repp=0
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
  iter<-0
  while (det(D1)!=0 & sgn1==0 & iter<5){
    iter<-iter+1
    print(iter)
    repp<-repp+1
    if (identical(sort(Sety),sort(a))){ #if (isTRUE(all.equal(d2[1:h,2],a))){
      sgn1=1}
    else {for (j in 1:h){y[,j]<-x0[,d2[j,2]]
    Sety[j]<-d2[j,2]}
      T1<-rowMeans(y)
      S1<-matrix(rep(0, p*p),p,p)
      for (j in 1:h){
        S1<-S1+(y[,j]-T1)%*%t(y[,j]-T1)/h}
      #for (j in 1:p){
      #if (S1[j,j]<0.00001){S1[j,j]<-0.001}}
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
#chatpn<-1+(matrix.trace(Rraw%^% 2)/(p^1.5))
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
#ctildepn<-1+(matrix.trace(Rtilde%^% 2)/(p^1.5))
ctildepn<-1+2*(p/(m*trR2w^0.5))
w<-c(rep(0,m)) 
TW<-c(rep(0,p))
SW<-matrix(rep(0, p*p),p,p)
for (i in 1:m){
  d2ref[i]<-d2muD0hat[i]/(1+(dnorm(qnorm(1-del))*((2*trR2w)^0.5)/(p*(1-del))))
  if ((d2ref[i]-p)/((2*ctildepn*trR2w)^0.5)-(8*trR3w)*(qnorm(1-alpha)^2-1)/(6*((2*trR2w)^1.5))<=qnorm(1-alpha)){w[i]=1}
  TW<-TW+(w[i]*x0[,i])
}
mu0hat<-TW/sum(w)
for (i in 1:m){
  SW<-SW+w[i]*(x0[,i]-mu0hat)%*%t(x0[,i]-mu0hat)/(sum(w)-1)}
Dhat<-diag(diag(SW))
Rho<-solve(sqrtm(Dhat))%*%SW%*%solve(sqrtm(Dhat))
trRho2<-matrix.trace(Rho%^% 2)-(p^2)/sum(w)
trR3<-matrix.trace(Rho%^%3)-3*p*trRho2/sum(w)-p^3/sum(w)^2
#ctilde2pn<-1+(matrix.trace(Rho%^% 2)/(p^1.5))
ctilde2pn<-1+2*(p/(m*trRho2^0.5))

qqnorm((d2ref[i]-p)/((2*ctildepn*trR2w)^0.5)-(8*trR3w)*(qnorm(1-alpha)^2-1)/(6*((2*trR2w)^1.5)))
qqline((d2ref[i]-p)/((2*ctildepn*trR2w)^0.5)-(8*trR3w)*(qnorm(1-alpha)^2-1)/(6*((2*trR2w)^1.5)))
#shapiro.test((d2ref[i]-p)/((2*ctildepn*trR2w)^0.5)-(8*trR3w)*(qnorm(1-alpha)^2-1)/(6*((2*trR2w)^1.5)))
ab<-0
#for (i in 1:p){ab[i]<-as.numeric(shapiro.test(x0[i,])[2])}
print(ab)
print(w)
aa<-(d2ref[1:m]-p)/((2*ctildepn*trR2w)^0.5)-(8*trR3w)*(qnorm(1-alpha)^2-1)/(6*((2*trR2w)^1.5))
write.xlsx(aa, "1.xlsx")

#Phase 2 monitoring
p=250 # number of variables
m=1567 # number of histoical (phase I) samples
alpha=0.005
ss<-sum(w[1:rr])
bb<-matrix(rep(0, p*ss),p,ss)
k=0
for (i in 1:rr){
  
  if (w[i]==1) {k<-k+1
  bb[,k]<-x0[,i]
    
  }
}
An1<-matrix(rep(0, p*p),p,p)
for (j in 1:k){
  An1<-An1+(bb[,j])%*%t(bb[,j])}
x2<-bb
An<-An1
#xbar2<-T1
m1<-ncol(x2)

x0<-x1
x00<-matrix(rep(0, p*m),p,m)
for (i in 1:p){
  x00[i,]<-x0[SV[i],]
}
x0<-x00
#S1<-cov(t(x1))
S1<-SW
#muhat<-rowMeans(x1)
muhat<-mu0hat
#D1<-diag(diag(S1))
D1<-diag(diag(SW))
#Rho<-solve(sqrtm(D1))%*%S1%*%solve(sqrtm(D1))
#trRho2<-matrix.trace(Rho%^% 2)-(p^2)/m
#trR3<-matrix.trace(Rho%^%3)-3*p*trRho2/m-p^3/m^2
w<-c(rep(0,m))
wzhang<-c(rep(0,m))
d<-p^2/trRho2
d2zhang<-c(rep(0,m))
rrr<-rr+1
for (i in rrr:1567){
  d2mdp[i]<-t(x0[,i]-muhat)%*%solve(D1)%*%(x0[,i]-muhat) 
  d2zhang[i]<-(1/(p))*d2mdp[i]
  if ((d2mdp[i]-p)/((2*trRho2)^0.5)-(8*trR3)*(qnorm(1-alpha)^2-1)/(6*((2*trRho2)^1.5))>qnorm(1-alpha)){w[i]=1}
  else {x2<-cbind(x2,x0[,i])
  An<-An+(x0[,i])%*%t(x0[,i])
  xbar2<-rowMeans(x2)
  m1<-ncol(x2)
  S2<-(An-xbar2%*%t(m1*xbar2)-(m1*xbar2)%*%t(xbar2)+(m1*xbar2%*%t(xbar2)))/(m1-1)
  D1<-diag(diag(S2))
  Rho<-solve(sqrtm(D1))%*%S2%*%solve(sqrtm(D1))
  trRho2<-matrix.trace(Rho%^% 2)-(p^2)/m1
  trR3<-matrix.trace(Rho%^%3)-3*p*trRho2/m1-p^3/m1^2
  muhat<-xbar2
  }
  
  if (d*d2zhang[i]>qchisq(1-alpha,d)){wzhang[i]=1}
  }
print(w[(rr+1):1567])
aa2<-(d2mdp[(rr+1):1567]-p)/((2*trRho2)^0.5)-(8*trR3)*(qnorm(1-alpha)^2-1)/(6*((2*trRho2)^1.5))
write.xlsx(aa2, "2.xlsx")
alpha=seq(1/1464,1,1/1464)
alpha1=seq(1/104,1,1/104)
plot(x=seq(-4,4,0.1),pnorm(seq(-4,4,0.1)),"l",lty=3,lwd=3,xlim=c(-5,6),ylim=c(0,1.1),xlab="X", ylab="Probability")
lines(ecdf(sort((d2mdp[1464:1567]-p)/((2*trRho2)^0.5),decreasing = TRUE)-(8*trR3)*(qnorm(1-alpha1)^2-1)/(6*((2*trRho2)^1.5))),col="blue")
lines(ecdf(sort((d2muD0hat-p)/((2*trRho2)^0.5),decreasing = TRUE)-(8*trR3)*(qnorm(1-alpha)^2-1)/(6*((2*trRho2)^1.5))),col="green")
ggcorrplot(Rho[1:50,1:50],sig.level=0.05 ,lab_size = 1.5, p.mat = NULL, insig = c("pch", "blank"), pch = 1, pch.col = "black", pch.cex =0.25,
           tl.cex = 9)
