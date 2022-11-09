library(MASS)
library(tidyverse)
####now look using BART####
library(dplyr)




##Simulation for Audit code##

####helper functions####





mBG_out<-c()
#rho_list<-c(rep(0.25, 5), rep(0.4, 5), rep(0.6, 5), rep(0.8,5), 0.6)
#gamma_list<-c(rep(seq(from=1, to=5, length.out=5), 4), .7)
rho_list<-c(0.25, 0.25, 0.25, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8)
gamma_list<-c(1,1.75, 2.5, 1, 1.75, 2.5, 1, 1.75, 2.5, 1, 1.75, 2.5)
ATE<-c()
tau.true<-c()
rho_guess<-c()
cor_list<-c()
RR<-c()
gamma_guess<-c()
B1.true<-c()
B0.true<-c()
RR_true<-c()
cor_rr_list<-c()
RMSE_rr_list<-c()
cor_list<-c()
RMSE_list<-c()
for (m in 1:length(gamma_list)){
  
  N <- 100000# Number of random samples
  set.seed(0) #for variation in the x's
  # Target parameters for univariate normal distributions
  x1=runif(N, -1,1)
  x2=runif(N, -1,1)
  x3=runif(N,-1,1)
  x4=runif(N,-1,1)
  x5=runif(N, -1,1)

  beta1= -0.2
  alpha1= 0.7  
  beta0= -0
  alpha0= -0.5
  mu1 <- beta0+beta1*(x1+x2+x3+x4+x5)
  mu2 <- alpha0+alpha1*(x1+x2+x3+x5)
  
  mu<-matrix(c(mu1, mu2), nrow=N, ncol=2)
  rho=rho_list[[m]]
  gamma=gamma_list[[m]]
  sigma <- matrix(c(1, rho,rho,1),
                  2) # Covariance matrix
  
  sim_data=t(sapply(1:N, function(i)mvrnorm(1, mu = mu[i,], Sigma = sigma )))
  #sim_data=mvrnorm(N, c(-1.6,-2.6), sigma)
  G=sapply(1:N, function(i)ifelse(sim_data[i,1]>=0, 1,0))
  
  
  B=sapply(1:N, function(i)ifelse(sim_data[i,2]>=-1*gamma*G[i], 1,0))
  
  covariates=data.frame(x1,x2,x3,x4,x5,B, G)
  
  vars=c('x1','x2','x3','x4','x5')
  library(GJRM)
  out<-gjrm(list(G ~x1+x2+x3+x4+x5, 
                 B~ G+x1+x2+x3+x5), 
            data=covariates,
            margins =c("probit", "probit"),
            Model = "B"
  )
  
  AT(out, 'G')
  X.int <- as.matrix(out$X2)
  ind.int<-(1:out$X2.d2)+out$X1.d2
  coef.int <- as.numeric(out$coefficients[ind.int])
  d0 <- d1 <- X.int
  d0[, 'G'] <- 0
  d1[, 'G'] <- 1
  eti1 <- d1 %*% coef.int
  eti0 <- d0 %*% coef.int
  p.int1 <- probm(eti1, out$margins[2], min.dn = out$VC$min.dn, 
                  min.pr = out$VC$min.pr, max.pr = out$VC$max.pr)$pr
  p.int0 <- probm(eti0, out$margins[2],  min.dn = out$VC$min.dn, 
                  min.pr = out$VC$min.pr, max.pr = out$VC$max.pr)$pr

  est.AT <- mean(p.int1, na.rm = TRUE) - mean(p.int0, na.rm = TRUE)
  n.sim=500
  bs <- rMVN(n.sim, mean = out$coefficients, sigma = out$Vb)
  eti1s <- d1 %*% t(bs[, ind.int])
  eti0s <- d0 %*% t(bs[, ind.int])
  gamma=qnorm(p.int1)-qnorm(p.int0)
  
  rho_guess[[m]]=summary(out)[]$theta
  tau.true[[m]] = pnorm(gamma+mu2) - pnorm(mu2)
  B1.true=pnorm(gamma+mu2)
  B0.true=pnorm(mu2)
  gamma_guess[[m]]<-out$coefficients['G']
  cor_list[[m]]<-cor(p.int1-p.int0, tau.true[[m]])
  
  ATE[[m]]<-AT(out, nm.end = "G", hd.plot = TRUE, n.sim = 1000)
  
  #ate_2[[i]]<-est.AT
  cor_list[[m]]<-cor(p.int1-p.int0, tau.true[[m]])
  RMSE <- function(m, o){
    sqrt(mean((m - o)^2))
  }
  
  RR_true[[m]]=mean(B1.true/B0.true)
  RR[[m]]=mean(p.int1/p.int0)
  cor_rr_list[[m]]=cor(B1.true/B0.true, p.int1/p.int0)
  RMSE_rr_list[[m]]=RMSE(B1.true/B0.true, p.int1/p.int0)
  RMSE_list[[m]]<-RMSE(p.int1-p.int0, tau.true[[m]])
  
  print(m)
}

RMSE <- function(m, o){
  sqrt(mean((m - o)^2))
}
trueate_list<-sapply(1:length(rho_list), function(i) mean(tau.true[[i]]))
ATE_est<-sapply(1:length(rho_list), function(i) ATE[[i]]$res[2])

RR_cor=sapply(1:length(rho_list), function(i) cor(RR_true[[i]], RR[[i]]))
RR_rmse=sapply(1:length(rho_list), function(i) RMSE(RR_true[[i]], RR[[i]]))


View(data.frame(true_ate=round(trueate_list,3), 
                ate_est=round(ATE_est,3),
                ate_cor=round(cor_list,3),
                ate_rmse=round(RMSE_list, 3),
                RR_true=round(unlist(RR_true), 3), 
                RR_guess=round(RR, 3),
                RR_cor=round(cor_rr_list, 3), 
                RR_rmse=round(RMSE_rr_list ,3), 
                gamma_true=gamma_list, 
                gamme_guess=gamma_guess, 
                rho_true=rho_list, 
                rho_guess=round(unlist(rho_guess),3)))

library(xtable)
xtable(data.frame(true_ate=round(trueate_list,3), 
                  ate_est=round(ATE_est,3),
                  ate_cor=round(cor_list,3),
                  ate_rmse=round(RMSE_list, 3),
                  RR_true=round(unlist(RR_true), 3), 
                  RR_guess=round(RR, 3),
                  RR_cor=round(cor_rr_list, 3), 
                  RR_rmse=round(RMSE_rr_list ,3), 
                  gamma_true=gamma_list, 
                  gamme_guess=gamma_guess, 
                  rho_true=rho_list, 
                  rho_guess=round(unlist(rho_guess),3)),  
       include.rownames=FALSE)
round(unlist(RR_true), 3)
RMSE <- function(m, o){
  sqrt(mean((m - o)^2))
}


ATE_est

