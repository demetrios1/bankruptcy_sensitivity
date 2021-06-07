library(MASS)
library(tidyverse)
####now look using BART####
library(dplyr)
library(dbarts)
library(bcf)
library(fastbart)


monotone_bart = function(y, z, x, xpred, nskip=5000, ndpost=5000, m = 50, n=N) {
  
  sort_ix = order(z, y)
  x = x[sort_ix,]
  z = z[sort_ix]
  y = y[sort_ix]
  
  n0 = sum(z==0)
  n00 = sum(z==0 & y==0)
  yobs = y
  yobs0 = rbinom(n, 1, 0.5)
  yobs[1:n00] = rbinom(n00, 1, 0.5) # These are DA variables
  yobs0[1:n00] = rbinom(n00, 1, 0.5) # These are DA variables
  yobs0[1:n00][yobs[1:n00]==1] = 0 # To statisfy DA constraints
  yobs0[z==0 & y==1] = 1
  yobs0 = yobs0[1:n0] 
  
  offset =  0#qnorm(mean(y[z==1]))
  offset0 = 0# qnorm(mean(flu$grp==0 & flu$fluy2==1)/mean(flu$grp==1 & flu$fluy2==1)) #<- wtf
  
  zz = offset + 3*yobs - 3*(1-yobs)
  z0 = offset0 + 3*yobs0 - 3*(1-yobs0)
  
  ################################################################################
  # MCMC
  ################################################################################
  
  set.seed(1022)
  
  xi = lapply(1:ncol(x), function(i) bcf:::.cp_quantile(x[,i]))
  fit.mono = bartRcppMono(yobs, zz, t(as.matrix(x)), t(xpred),
                          yobs0, z0, t(as.matrix(x)),t(xpred),
                          n00,
                          xi,
                          nskip, ndpost, m, 3.0,
                          offset, offset0)
  
  #xpred.exp = rbind(data.frame(xpred, z=1), data.frame(xpred, z=0))
  #fit =  bart(cbind(x, z), y, xpred.exp,
  #               nskip=5000, ndpost=5000,
  #               ntree=50, usequants=T)
  
  
  # P(Y|Z=1, X)
  pr1 = pnorm(fit.mono$postpred)
  
  # P(Y|Z=0, X)
  pr0 = pr1*pnorm(fit.mono$postpred0)
  
  return(list(pr0 = pr0, pr1 = pr1))
  
}

BARTpred=function(df, treat=G, Outcome=B,vars, mono=T){
  covtreat=df%>%filter(G==1)
  covcontrol=df%>%filter(G==0)
  covtreat$G=as.factor(covtreat$G)
  #return(df[vars])
  #return(as.factor(covtreat$outcome))
  #case 1, the treatment
  if (mono==F){
   bart1=bart(covtreat[vars],as.factor(covtreat$B), df[vars],ndpost = 5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  
  #case 2 control
  bart2=bart(covcontrol[vars],as.factor(covcontrol$B),df[vars],ndpost =5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  
  #case 3 propensity
  bart3=bart(df[vars], as.factor(df$G),df[vars],ndpost = 5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  pred1 = colMeans(pnorm(bart1$yhat.test)) 
  pred2=colMeans(pnorm(bart2$yhat.test))
  pred3=colMeans(pnorm(bart3$yhat.test))
  }
  else if(mono==T){
  xtest=df[vars]
  xtraincontrol=covcontrol[vars]
  xtraintreat=covtreat[vars]
  ytrain0    = as.factor( covcontrol$B )
  ytrain1    = as.factor( covtreat$B )
  
  # mono fits
  
  bart_mono = monotone_bart(y = as.numeric(c(ytrain1, ytrain0)==1),
                            z = 1-c(rep(1, length(ytrain1)), rep(0, length(ytrain0))),
                           x = rbind(xtraintreat, xtraincontrol),
                            xpred = xtest, nskip = 1000, ndpost = 1000,m=100)
  
  
  
  bart3 = bart(df[vars], as.factor(df$G),df[vars], ndpost = 1000, nskip = 1000, ntree=100, usequants=T,  verbose=F)
  
  
  
  #use this method for prediction on binary
 
  pred1= 1-colMeans(bart_mono$pr0)
  pred2= 1-colMeans(bart_mono$pr1)
  pred3=colMeans(pnorm(bart3$yhat.test))
  }
else{
  print('mono argument is T or F Boolean')
}
  
  
  
  
  expoutcomesfun	   =cbind( df[,], data.frame(  treatpred = pred1, notreatpred=pred2), propensity=pred3 )
  
  expoutcomesfun2    = rbind( NULL, expoutcomesfun )
  outcomenew=expoutcomesfun2[,c('propensity','notreatpred','treatpred')]
  
  
  ####need the joints####
  
  eps     = 1e-6
  
  outcomenew[,c("prop", "Y_D0", "Y_D1")] = 0.5*eps + (1-eps)*outcomenew[,c('propensity','notreatpred','treatpred')]
  
  outcomenew[,"ProbY1D1"]=outcomenew[,"prop"]*outcomenew[,"treatpred"]
  outcomenew[,"ProbY1D0"] = (1-outcomenew[,"prop"])*outcomenew[,"notreatpred"]
  outcomenew[,"ProbY0D1"] = outcomenew[,"prop"]*(1-outcomenew[,"treatpred"])
  
  #for error analysis, not super necessary
  indexes=seq(from=1, to=length(outcomenew$treatpred), by=1)
  outcomenew=as.data.frame(cbind(indexes, outcomenew))
  outcomenew$indexes=as.numeric(outcomenew$indexes)
  row.names(outcomenew)<-NULL
  newoutcome4 =as.matrix(outcomenew)
  return(outcomenew)
}



##Simulation for Audit code##

library(foreach)
library(doParallel)
####helper functions####

np = detectCores()-1
cl = makeCluster(np)
n_cores <- detectCores() - 1

registerDoParallel(cores = 10)#n_cores)





mBG_out<-c()
rho_list<-c(rep(0.25, 5), rep(0.4, 5), rep(0.6, 5), rep(0.8,5), 0.6)
gamma_list<-c(rep(seq(from=1, to=5, length.out=5), 4), .7)
mu2_list<-c()
for (m in 1:length(gamma_list)){
  N <- 5000 # Number of random samples
  #set.seed(123) #for variation in the x's
  # Target parameters for univariate normal distributions
  x1=runif(N, -2,2)
  x2=runif(N, -2,2)
  x3=runif(N,-2,2)
  x4=runif(N,-2,2)
  x5=runif(N, -2,2)
  
  beta1= -0.2
  alpha1= 0.7
  beta0= -0
  alpha0= -0.5
  mu1 <- beta0+beta1*(x1+x2+x3+x4+x5)
  mu2 <- alpha0+alpha1*(x1+x2+x3+x4+x5)
  mu2_list[[m]]<-mu2
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
  intframe=BARTpred(covariates, treat=G, Outcome=B,vars, mono=T)
  

  #plot(pnorm((mu2+gamma)[sim_data[,1]>0]), intframe$Y_D1[sim_data[,1]>0])
  library(truncnorm)
  pB.G1 = sapply(1:N, function(j) mean(pnorm(mu2[j] + rho*rtruncnorm(10000,0,+Inf,mu1[j]) + gamma)))
  plot(pB.G1, intframe[, 'Y_D1'])
  plot(pnorm(mu2), intframe[,'Y_D0'])
  
  #cl = makeCluster(np)
  #registerDoParallel(cl)
  newoutcome=as.matrix(intframe)
  
  ptm = proc.time()
  
  mBG_out[[m]] = foreach( i = 1:dim(newoutcome)[1],  .combine=rbind )%dopar%{
    
    f=function(u){
      dnorm(u, mean=0, sd=sqrt(rho/(1-rho)))  #prob could just use var instead of the square root.
    }
    
    # Probability to fit.
    vProbBG     = newoutcome[i,c( "ProbY1D1", "ProbY1D0", "ProbY0D1" )] 
    
    
    # Starting value
    
    start0 = c( newoutcome[i,"ProbY1D1"], newoutcome[i,"ProbY1D0"],newoutcome[i,"ProbY0D1"] )
    
    names(start0) = c("y1","y0","d")
    
    
    
    #better global definition.
    optimdist=function(vals){
      y1 = vals[1]
      y0 = vals[2]
      d = vals[3]
      a=integrate( function(u) pnorm(y1+u)*pnorm(d+u)*f(u), lower = -Inf, upper = Inf )$value
      b=integrate(function(u) pnorm(y0+u)*(1-pnorm(d+u))*f(u), lower = -Inf, upper = Inf )$value
      c=integrate(  function(u) (1-pnorm(y1+u))*pnorm(d+u)*f(u), lower = -Inf, upper = Inf )$value 
      #return(c(a,b,c))
      return(sum((qnorm(c(a,b,c))-qnorm(vProbBG))^2))
    }
    
    #calculate the treatment
    treatment_Compute = function( vb){
      integrate( function(u) (pnorm(vb[1]+u) - pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value
      
    }
    
    counterfac_Compute = function(vb){
      y1 =(integrate( function(u) (pnorm(vb[1]+u))*f(u), lower = -Inf, upper = Inf )$value)
      y0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
      
      return(c(y1, y0))
      
    }
    vBG_out     = try(optim( qnorm(start0), optimdist,  control = list(maxit = 500) )	
                      , silent=TRUE)
    if ( class( vBG_out ) == "try-error" ){
      out        = as.numeric(newoutcome[i,c("indexes")])
      out        = c( out, rep(NA,5) ) 
    } else {
      out        = as.numeric(newoutcome[i,c("indexes")] )
      vBG        = vBG_out$par 
      tau      = treatment_Compute( vb = vBG[c(1,2)])
      
      counter=counterfac_Compute(vb = vBG[c(1,2)])
      
      out        = c( vBG, tau,vBG_out$convergence, vBG_out$value,counter) 
    }
    return( out)
  }
  
  print( proc.time() - ptm )

  colnames(mBG_out[[m]]) = c( "y1", "y0", "d", "tau", "convergence", "fnvalue" ,'B1','B0')
  print(m)
}



plot(mBG_out[[1]][, 'B1'], pnorm(gamma_list[[1]]+mu2_list[[1]]))
plot(mBG_out[[1]][, 'B0'], pnorm(mu2_list[[1]]))

p1 = mBG_out[[1]][, 'B1']
p0 = mBG_out[[1]][, 'B0']

tau.est = p1 - p0

tau.est_mean=sapply(1:length(gamma_list),function(i)mean(mBG_out[[i]][, 'B1']-mBG_out[[i]][, 'B0']))
tau.true = pnorm(gamma_list[1]+mu2_list[[1]]) - pnorm(mu2_list[[1]])
cor(tau.est,tau.true)


tau.true_mean=sapply(1:length(gamma_list),function(i)mean(pnorm(gamma_list[i]+mu2_list[[i]])-pnorm(mu2_list[[i]])))

plot(tau.true,tau.est,pch=20, col='dodgerblue4', main='Monotone BART Fit')
abline(0,1,col='red',lwd=3)
text(0.08, 0.45, paste0('correlation= ',round(cor(tau.est, tau.true),3)))
#write.csv(data.frame(tau.true, tau.est), 'justbart.csv')


#as.data.frame(mBG_out[[1]])%>%
#  ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),fill='dodgerblue4', color='white')+
#  geom_vline(aes(xintercept=mean(tau.true)))

tau_mean<-sapply(1:length(gamma_list), function(i)mean(mBG_out[[i]][,'tau']))
correlation<-sapply(1:4, 
        function(i)cor(pnorm(gamma_list[i]+mu2_list[[i]])-pnorm(mu2_list[[i]]),mBG_out[[i]][, 'B1']-mBG_out[[i]][, 'B0'] ))
ecdf_tau_true<-sapply(1:length(gamma_list),function(i)ecdf(mBG_out[[i]][, 'tau'])(tau_mean[[i]]))
data.frame(tau_est=tau.est_mean,tau_true=tau_mean, correlation,ecdf_tau, gamma=gamma_list,rho=rho_list)
library(xtable)
xtable(data.frame(tau_est=tau.est_mean,tau_true=tau_mean, correlation,ecdf_tau, gamma=gamma_list,rho=rho_list))
correlation
