
rm(list = ls())
qlist<-rep(c(.1,.1, .1, .1, .5, .5, .9, .9),2) 
qlistfalse<-c(qlist[1:8], .9, .9, .9, .5, .1, .9, .1, .5)

varlist=rep(c(3, .5, 1,1,1,1,1,1),2)
varlistfalse=rep(c(3, .5, 1,1,1,1,1,1),2)


neg_2_cor<-c()
neg_2_inc<-c()
pos_2_cor<-c()
pos_2_inc<-c()

#meanu=c(0)
#std=c(1)

treat_list<-c()
tau_est_list<-c()
n=25000
N=25000
mc=10000

for (i in 9:length(qlist)){
  b0=function(x1,x2,x3,x4,x5,x6,x7,x8){
    val = -1.75 + x5 + x1*sin(2*x6)
    return(val)
  }
  
  b1=function(x1,x2,x3,x4,x5,x6,x7,x8){
    val = 1.5 + b0(x1,x2,x3,x4,x5,x6,x7,x8)
    return(val)
  }
  
  
  g=function(x1,x2,x3,x4,x5,x6,x7,x8){
    val= 0.25 + 0.5*b0(x1,x2,x3,x4,x5,x6,x7,x8) + x2
    return(val)
  }
  set.seed(12345)
  a = 2
  x1 = runif(n,-a,a)
  x2 = runif(n,-a,a)
  x3 = runif(n,-a,a)
  x4 = runif(n,-a,a)
  x5 = runif(n,-a,a)
  x6 = runif(n,-a,a)
  x7 = runif(n,-a,a)
  x8 = runif(n,-a,a)
  

  
  ####given q and variance, find the sig we need####
  q=qlist[[i]]
  var=varlist[[i]]
  a=0.3633802
  b=0.7978846
  varexp=(1-q)*a*((1-q)/q)^2 + q*a + (b*(1 + (1-q)/q))^2*q*(1-q)
  
  sig<-sqrt(var/varexp)
  
  a=0.3633802
  b=0.7978846

  varexp=(1-q)*a*((1-q)/q)^2 + q*a + (b*(1 + (1-q)/q))^2*q*(1-q)
  ####given the false q and the false variance, find the false s parameter####
  qfalse=qlistfalse[[i]]
  varfalse=varlistfalse[[i]]
  a=0.3633802
  b=0.7978846
  varexpfalse=(1-qfalse)*a*((1-qfalse)/qfalse)^2 + qfalse*a + (b*(1 + (1-qfalse)/qfalse))^2*qfalse*(1-qfalse)
  
  sigfalse<-sqrt(varfalse/varexpfalse)
  
  
  
  ind = rbinom(n,1,q)
  u = rep(NA,n)
  
  
  
  library(truncnorm)
  u[ind==1] = rtruncnorm(sum(ind==1),a = -Inf, b = 0, sd = sig)
  u[ind==0] = rtruncnorm(sum(ind==0),a = 0, b = +Inf,sd = sig*(1-q)/q)  #prob 1-q over 0, 1/s and s cancel
  
  
  G = rbinom(n, 1,pnorm(g(x1,x2,x3,x4,x5,x6,x7,x8) + u))
  B = rep(NA,n)
  
  B[G == 1] = rbinom(sum(G==1),1, pnorm(b1(x1,x2,x3,x4,x5,x6,x7,x8)[G == 1] + u[G == 1]))
  B[G ==0] = rbinom(sum(G==0),1,pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8)[G == 0] + u[G == 0]))
  B1=pnorm(b1(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  B0=pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  B0prime=pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8) )
  #  G1=pnorm(g(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  print(table(B,G))

  
  fcorrect=function(u){
    val = rep(NA, length(u))
    val[u < 0] = 2*q*dnorm(u[u<0], sd = sig)
    val[u>0] = 2*(1-q)*dnorm(u[u>0],0,sd = sig*(1-q)/q)
    return(val)
  }
  f=function(u){
    val = rep(NA, length(u))
    val[u < 0] = 2*qfalse*dnorm(u[u<0], sd = sigfalse)
    val[u>0] = 2*(1-qfalse)*dnorm(u[u>0],0,sd = sigfalse*(1-qfalse)/qfalse)
    return(val)
  }
  print(q)
  ind = rbinom(mc,1,q)
  x = rep(NA,mc)
 
  library(truncnorm)

  x[ind==1] = rtruncnorm(sum(ind==1),a = -Inf, b = 0, sd = sig)
  x[ind==0] = rtruncnorm(sum(ind==0),a = 0, b = +Inf,sd = sig*(1-q)/q)  #prob 1-q over 0, 1/s and s cancel
  
  print(qfalse)
  indfalse = rbinom(mc,1,qfalse)
  xfalse = rep(NA,mc)
  
  library(truncnorm)
  xfalse[indfalse==1] = rtruncnorm(sum(indfalse==1),a = -Inf, b = 0, sd = sigfalse)
  xfalse[indfalse==0] = rtruncnorm(sum(indfalse==0),a = 0, b = +Inf,
                              sd = sigfalse*(1-qfalse)/qfalse)  #prob 1-q over 0, 1/s and s cancel

####double check
  print(paste0('true variance of f correct should be ', var, 'mc generated is ', var(x)))
  print(paste0('true variance of f incorrect should be ', var, 'mc generated is ', var(xfalse)))
  var(xfalse)
  ecdf_cor=ecdf(x)
  ecdf_incor=ecdf(xfalse)
  neg_2_cor[[i]]=ecdf_cor(-2)
  neg_2_inc[[i]]=ecdf_incor(-2)
  pos_2_cor[[i]]=ecdf_cor(2)
  pos_2_inc[[i]]=ecdf_incor(2)

  xgrid = seq(-4,4,length.out = 500)
  plot(xgrid,fcorrect(xgrid), type='l', col='red')
  lines(xgrid, f(xgrid), col='blue', legend=T)
  
  print(i)
}