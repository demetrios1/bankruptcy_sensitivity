##packages##
library(data.table)
library(readr)
library(dplyr)
library(tidyverse)


####init####

#make bigger later
N=25000
n = 25000
#### covariates####


####confounding on D and Y (private info analogue for audits)####
#can run this for various U's
 ufun=function(mean1, sd1){
   return(rnorm(n, mean=mean1, sd=sd1))
 }

ATE<-c()
means=c(-1,1)
sds=c(.5, .5,.5)
j=.9
treat_list<-c()
tau_est_list<-c()
naive_list<-c()
treat_list_bart<-c()
tau_est_list_bart<-c()
naive_list_bart<-c()
prob_treat<-c()
true_loc_ecdf<-c()
tau.true<-c()
treateffect<-c()
U<-c()
RR<-c()
RR_true<-c()
cor_list<-c()
RMSE_list<-c()
ate_2<-c()
#sanity check
#globally define the marginal distribution of U that we integrate over
meanu=c(rep(c(-3,-2,-1,0,1,2,3), 2))
std=c(1,1,1,0.1,1,1,1,rep(0.5, 7))

meanu=c(0,0,0,0, -1, 1, -2, 2)
std=c(1,1.5, 2, 2.5, 1,2, 2, 1)
for (i in 1:length(meanu)){
  n = 25000
  set.seed(12345)
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
  
  a = 1
  x1 = runif(n,-a,a)
  x2 = runif(n,-a,a)
  x3 = runif(n,-a,a)
  x4 = runif(n,-a,a)
  x5 = runif(n,-a,a)
  x6 = runif(n,-a,a)
  x7 = runif(n,-a,a)
  x8 = runif(n,-a,a)
  
  u = std[[i]]*rnorm(n, mean=meanu[[i]])
  
  G = rbinom(n, 1,pnorm(g(x1,x2,x3,x4,x5,x6,x7,x8) + u))
  B = rep(NA,n)
  
  B[G == 1] = rbinom(sum(G==1),1, pnorm(b1(x1,x2,x3,x4,x5,x6,x7,x8)[G == 1] + u[G == 1]))
  B[G ==0] = rbinom(sum(G==0),1,pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8)[G == 0] + u[G == 0]))
  # B1=pnorm(b1(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  #  B0=pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  #  G1=pnorm(g(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  print(table(B,G))
  f=function(u){
    
    dnorm(u, mean=meanu[[i]], sd=std[[i]])
  }
  
  U[[i]] = rnorm(n, mean=meanu[[i]], sd=std[[i]])
  u=U[[i]]
  G = rbinom(n, 1,pnorm(g(x1,x2,x3,x4,x5,x6,x7,x8) + u))
  B = rep(NA,n)
  
  B[G == 1] = rbinom(sum(G==1),1, pnorm(b1(x1,x2,x3,x4,x5,x6,x7,x8) + u))
  B[G ==0] = rbinom(sum(G==0),1,pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8) + u))
  covariates=data.frame(x1,x2,x3,x4,x5,x6,x7,x8,B, G)
  vars=c('x1','x2','x3','x4','x5','x6','x7','x8')
 

  

  tau.true = rep(NA,n)
  for (k in 1:n){
    #print(k)
    tau.true[k] = integrate(function(u) (pnorm(b1(x1[k], x2[k],x3[k],x4[k],x5[k],x6[k],x7[k],x8[k]) + u) - 
                                           pnorm(b0(x1[k], x2[k],x3[k],x4[k],x5[k],x6[k],x7[k],x8[k])+u))*f(u),lower = -Inf, upper = +Inf)$value
  }
  B1.true = rep(NA,n)
  for (k in 1:n){
    #print(k)
    B1.true[k] = integrate(function(u) (pnorm(b1(x1[k], x2[k],x3[k],x4[k],x5[k],x6[k],x7[k],x8[k]) + u) )*f(u),lower = -Inf, upper = +Inf)$value
  }
  
  B0.true = rep(NA,n)
  for (k in 1:n){
    #print(k)
    B0.true[k] = integrate(function(u) ( pnorm(b0(x1[k], x2[k],x3[k],x4[k],x5[k],x6[k],x7[k],x8[k])+u))*f(u),lower = -Inf, upper = +Inf)$value
  }
  treateffect[[i]]<-mean(tau.true)
  covariates=data.frame(x1,x2,x3,x4,x5,x6,x7,x8, B=B, G=G, U=unlist(U[[i]]))
  library(GJRM)
  
  ##no smoothing##
  out<-gjrm(list(G ~x1+x2+x3+x4+x5+x6+x7+x8, 
                 B ~ G+x1+x2+x3+x4+x5+x6+x7+x8), 
            data=covariates,
           margins =c("probit", "probit"),
            Model = "B"
  )
  
##with smoothing##
  #  out<-gjrm(list(G ~s(x1)+s(x2)+s(x3)+s(x4)+s(x5)+s(x6)+s(x7)+s(x8), 
  #               B ~ G+s(x1)+s(x2)+s(x3)+s(x4)+s(x5)+s(x6)+s(x7)+s(x8)), 
  #         data=covariates,
  #          margins =c("probit", "probit"),
  #          Model = "B"
# )
  
  ATE[[i]]<-AT(out, nm.end = "G", hd.plot = TRUE, n.sim = 500)  
  print(paste0('true treatment: ', round(treateffect[[i]],3), ' our estimate: ', round(ATE[[i]]$res[2],3 )))
 ####for the no smoothing case####
 ##replicates the docs of gjrm to see how they get their ITEs
  X.int <- as.matrix(out$X2)
  ind.int<-(1:out$X2.d2)+out$X1.d2
  coef.int <- as.numeric(out$coefficients[ind.int])
  d0 <- d1 <- X.int
  d0[, 'G'] <- 0
  d1[, 'G'] <- 1
  eti1 <- d1 %*% coef.int
  eti0 <- d0 %*% coef.int
   p.int1 <- probm(eti1, out$margins[2])$pr
  p.int0 <- probm(eti0, out$margins[2])$pr
  ####uncomment here if using the newer version
 # p.int1 <- probm(eti1, out$margins[2], min.dn = out$VC$min.dn, 
  #                min.pr = out$VC$min.pr, max.pr = out$VC$max.pr)$pr
  #p.int0 <- probm(eti0, out$margins[2],  min.dn = out$VC$min.dn, 
   #               min.pr = out$VC$min.pr, max.pr = out$VC$max.pr)$pr
  
  est.AT <- mean(p.int1, na.rm = TRUE) - mean(p.int0, na.rm = TRUE)
  n.sim=500
  bs <- rMVN(n.sim, mean = out$coefficients, sigma = out$Vb)
  eti1s <- d1 %*% t(bs[, ind.int])
  eti0s <- d0 %*% t(bs[, ind.int])
 

  ate_2[[i]]<-est.AT
  cor_list[[i]]<-cor(p.int1-p.int0, tau.true)
  RMSE <- function(m, o){
    sqrt(mean((m - o)^2))
  }
  RR_true[[i]]=B1.true/B0.true
  RR[[i]]=p.int1/p.int0
  RMSE_list[[i]]<-RMSE(p.int1-p.int0, tau.true)
  print(i)
}
RMSE <- function(m, o){
  sqrt(mean((m - o)^2))
}

RR_truth=sapply(1:8, function(i)mean(RR_true[[i]]))
RR_truth_log=sapply(1:8, function(i)mean(log(RR_true[[i]])))
RR_guess=sapply(1:8, function(i)mean(RR[[i]]))
RR_guess_log=sapply(1:8, function(i)mean(log(RR[[i]])))
RR_truth_log
RR_truth
RR_guess_log
RR_guess
RR_truth
RMSE_RR=sapply(1:8, function(i) RMSE(RR_true[[i]], RR[[i]]))
RMSE_RR_log=sapply(1:8, function(i) Metrics::rmsle(RR_true[[i]], RR[[i]]))
RR_guess
ATE_est<-sapply(1:8, function(i) ATE[[i]]$res[2])
ATE_est
View(data.frame(round(unlist(treateffect),3), round(ATE_est,5), round(unlist(cor_list), 3), 
                round(unlist(ate_2),5)))


sapply(1:8 ,function(i)mean(tau.true[[i]]))
trueate_list<-unlist(treateffect)


RMSE(ATE_est, trueate_list)

plot(ATE_est, trueate_list, main='non-lin DGP, BP reg')
abline(0,1, col='red')
unlist(RMSE_list)

