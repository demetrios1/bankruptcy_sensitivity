
##packages##
library(data.table)
#install.packages('dbarts')
library(dbarts)
library(rmutil)
library(readr)
library(dplyr)
library(tidyverse)
library(foreach)
library(bcf)
library(fastbart)
#library(GJRM)
print('to avoid error message')
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

BARTpred=function(df, treat=G, Outcome=B,vars){
  covtreat=df%>%filter(G==1)
  covcontrol=df%>%filter(G==0)
  covtreat$treat=as.factor(covtreat$G)
  
  #return(df[vars])
  #return(as.factor(covtreat$outcome))
  #case 1, the treatment
  #bart1=bart(covtreat[vars],as.factor(covtreat$B), df[vars],ndpost = 5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  
  #case 2 control
  #bart2=bart(covcontrol[vars],as.factor(covcontrol$B),df[vars],ndpost =5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  
  #case 3 propensity
  bart3=bart(df[vars], as.factor(df$G),df[vars],ndpost = 2000, nskip = 2000,ntree=100,
             verbose=F,usequants = TRUE,numcut = 1000)
  
  
  xtest=df[vars]
  xtraincontrol=covcontrol[vars]
  xtraintreat=covtreat[vars]
  ytrain0    = as.factor( covcontrol$B )
  ytrain1    = as.factor( covtreat$B )
  
  # mono fits
  
  bart_mono = monotone_bart(y = as.numeric(c(ytrain1, ytrain0)==1),
                            z = 1-c(rep(1, length(ytrain1)), rep(0, length(ytrain0))),
                            x = rbind(xtraintreat, xtraincontrol),
                            xpred = xtest, nskip = 2000, ndpost = 2000,m=100)
  
  
  
  #bart3 = bart(df[vars], as.factor(df$G),df[vars], ndpost = 5000, nskip = 2000, ntree=100, usequants=T, numcuts=1000, verbose=F)
  
  
  
  #use this method for prediction on binary
  #pred1 = colMeans(pnorm(bart1$yhat.test)) 
  #pred2=colMeans(pnorm(bart2$yhat.test))
  pred1= 1-colMeans(bart_mono$pr0)
  pred2= 1-colMeans(bart_mono$pr1)
  pred3=colMeans(pnorm(bart3$yhat.test))
  
  # pred1 = colMeans(pnorm(bart1$yhat.test)) 
  #  pred2=colMeans(pnorm(bart2$yhat.test))
  
  
  
  
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





#sanity check
#globally define the marginal distribution of U that we integrate over
meanu=c(0,0,0,0,-1, 1,  -2, 2)
std=c(.1,1, 2, 2.5,  1,2,  2, 1)
misspecstd=c(1.2, 1.75, 2.5,2, 1.3, 2.4, 2.3, 1.3 )

qlist=c(0.75, 0.25)
siglist=c(1.25, 0.50)



treat_list<-c()
tau_est_list<-c()
n=25000
N=25000
  ##Simulation for Audit code##
  library(foreach)
  library(doParallel)
  
  #np = 7
  #cl = makeCluster(np)
  n_cores <- detectCores() - 1
  registerDoParallel(cores = n_cores)
for (i in 1:2){
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
  ##if sharkfin
  library(truncnorm)
  q=qlist[[i]]
  sig=siglist[[i]]
  ind = rbinom(n,1,q) 
  u[ind==1]=rtruncnorm(sum(ind==1), a=-Inf, b=0, sd=sig)
  u[ind==0]=rtruncnorm(sum(ind==0), a=0, b=+Inf, sd=sig*(1-q)/q) 
  G = rbinom(n, 1,pnorm(g(x1,x2,x3,x4,x5,x6,x7,x8) + u))
  B = rep(NA,n)
  
  B[G == 1] = rbinom(sum(G==1),1, pnorm(b1(x1,x2,x3,x4,x5,x6,x7,x8)[G == 1] + u[G == 1]))
  B[G ==0] = rbinom(sum(G==0),1,pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8)[G == 0] + u[G == 0]))
  B1=pnorm(b1(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  B0=pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  B0prime=pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8) )
  #  G1=pnorm(g(x1,x2,x3,x4,x5,x6,x7,x8) + u)
  print(table(B,G))
  
  #sharkfin=function( beta){
  #  ifelse(beta>0, (2*(1-q)*(1/s))*dnorm(beta/s),2*q*(dnorm(beta)))
  #}
  
  
  sharkfin = function(z){
    val = rep(NA, length(z))
    val[z < 0] = 2*q*dnorm(z[z<0], sd = sig)
    val[z>0] = 2*(1-q)*dnorm(z[z>0],0,sd = sig*(1-q)/q)
    return(val)
    
  }
  fcorrect=function(u){
    #ifelse(u>0, (2*(1-q)*(1/s))*dnorm(u/s),2*q*(dnorm(u)))
    val = rep(NA, length(u))
    val[u < 0] = 2*q*dnorm(u[u<0], sd = sig)
    val[u>0] = 2*(1-q)*dnorm(u[u>0],0,sd = sig*(1-q)/q)
    return(val)
  }
  f=function(u){
    #ifelse(u>0, (2*(1-qfalse)*(1/sfalse))*dnorm(u/sfalse),2*qfalse*(dnorm(u)))
    val = rep(NA, length(u))
    val[u < 0] = 2*q*dnorm(u[u<0], sd = sig)
    val[u>0] = 2*(1-q)*dnorm(u[u>0],0,sd = sig*(1-q)/q)
    return(val)
    # dnorm(u, 0,1)
  }
  
  
 # fcorrect=function(u){
#    dnorm(u, mean=meanu[[i]], sd=std[[i]])
#  }
#  f=function(u){
#    return(dnorm(u, mean=meanu[[i]], sd=std[[i]]))
   # dlaplace(u, meanu[[i]], misspecstd[[i]])
 #   dnorm(u, -1, 1)
#    dnorm(u, mean=meanu[[i]], sd=std[[i]])
#  }



  
  tau.true = rep(NA,n)
  for (k in 1:n){
    #print(k)
    tau.true[k] = integrate(function(u) (pnorm(b1(x1[k], x2[k],x3[k],x4[k],x5[k],x6[k],x7[k],x8[k]) + u) - 
                                           pnorm(b0(x1[k], x2[k],x3[k],x4[k],x5[k],x6[k],x7[k],x8[k])+u))*fcorrect(u),
                            lower = -Inf, upper = +Inf)$value
  }
  
  B1.true = rep(NA,n)
  for (k in 1:n){
    #print(k)
    B1.true[k] = integrate(function(u) (pnorm(b1(x1[k], x2[k],x3[k],x4[k],x5[k],x6[k],x7[k],x8[k]) + u) )*fcorrect(u),
                           lower = -Inf, upper = +Inf)$value
  }
  
  B0.true = rep(NA,n)
  for (k in 1:n){
    #print(k)
    B0.true[k] = integrate(function(u) ( pnorm(b0(x1[k], x2[k],x3[k],x4[k],x5[k],x6[k],x7[k],x8[k])+u))*fcorrect(u),
                           lower = -Inf, upper = +Inf)$value
  }
  #true_ate=mean(pnorm(b1(x1,x2,x3,x4,x5,x6,x7,x8) + U)-pnorm(b0(x1,x2,x3,x4,x5,x6,x7,x8) + U))
  
  covariates=data.frame(x1,x2,x3,x4,x5,x6,x7,x8,B, G)
  vars=c('x1','x2','x3','x4','x5','x6','x7','x8')
  
  intframe=BARTpred(covariates, treat=G, Outcome=B,vars)
  RR=intframe$Y_D1/intframe$Y_D0
  E_val=RR+sqrt(RR*(RR-1))
  
  #covariates=data.frame(x1,x2,x3,x4,x5,x6,x7,x8, epsy,outcome$Y, treat)
  #vars=c('x1','x2','x3','x4','x5','x6','x7','x8', 'epsy')
  #intframe=BARTpred(covariates, treat=treat, Outcome=outcome$Y,vars)
  
  
  
  #registerDoParallel(cl)
  # intframe$ProbY1D1=B1
  #  intframe$ProbY1D0=B0
  #  intframe$ProbY0D1=G1
  newoutcome=as.matrix(intframe)

  
  


  
  print( sd )
  ptm = proc.time()
  
  mBG_out = foreach( m = 1:dim(newoutcome)[1],  .combine=rbind )%dopar%{
    # Probability to fit.
    vProbBG     = newoutcome[m,c( "ProbY1D1", "ProbY1D0", "ProbY0D1" )] 
    
    
    # Starting value
    
    start0 = c( newoutcome[m,"ProbY1D1"], newoutcome[m,"ProbY1D0"],newoutcome[m,"ProbY0D1"] )
    
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
      #vProbBG is lhs, global variable
      return(sum((qnorm(c(a,b,c))-qnorm(vProbBG))^2))
    }
    optimdist_constraint=
      function(vals){
        y1 = vals[1]
        y0 = vals[2]
        d = vals[3]
        #replace y1 with y1+y0
        a=integrate( function(u) pnorm((exp(y1)+y0)+u)*pnorm(d+u)*f(u), lower = -Inf, upper = Inf )$value
        b=integrate(function(u) pnorm(y0+u)*(1-pnorm(d+u))*f(u), lower = -Inf, upper = Inf )$value
        c=integrate(  function(u) (1-pnorm((exp(y1)+y0)+u))*pnorm(d+u)*f(u), lower = -Inf, upper = Inf )$value 
        #return(c(a,b,c))
        #vProbBG is lhs, global variable
        return(sum((qnorm(c(a,b,c))-qnorm(vProbBG))^2))
      }
    #calculate the treatment
    treatment_Compute = function( vb){
      integrate( function(u) (pnorm(vb[1]+u) - pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value
      
    }
    treatment_Compute_constraint = function( vb){
      integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u) - pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value
      
    }
    counterfac_Compute = function(vb){
      y1 =(integrate( function(u) (pnorm(vb[1]+u))*f(u), lower = -Inf, upper = Inf )$value)
      y0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
      
      return(c(y1, y0))
      
    }
    
    counterfac_Compute_constraint = function(vb){
      y1 =(integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
      y0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
      
      return(c(y1, y0))
      
    }

    vBG_out     = try(optim( qnorm(start0), optimdist_constraint,  control = list(maxit = 500) 
    )
    , silent=TRUE)
    if ( class( vBG_out ) == "try-error" ){
      out        = as.numeric(newoutcome[i,c("indexes")])
      out        = c( out, rep(NA,5) ) 
    } else {
      out        = as.numeric(newoutcome[i,c("indexes")] )
      vBG        = vBG_out$par 
       #tau      = treatment_Compute( vb = vBG[c(1,2)])
      
        #counter=counterfac_Compute(vb = vBG[c(1,2)])
      
      tau      = treatment_Compute_constraint( vb = vBG[c(1,2)])
      
      counter=counterfac_Compute_constraint(vb = vBG[c(1,2)])
      
      out        = c( vBG, tau,vBG_out$convergence, vBG_out$value,counter) 
    }
    return( out)
  }
  
  print( proc.time() - ptm )
 
  
  
  colnames(mBG_out) = c( "y1", "y0", "d", "tau", "convergence", "fnvalue" ,'B1','B0')
  as.data.frame(mBG_out)%>%
    summarize(mean(tau,na.rm=T))
  tau_est=mean(mBG_out[, 'tau'], na.rm=T)
  as.data.frame(mBG_out)%>%
    ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),fill='dodgerblue4', color='white')
  print(paste0('true treatment: ', round(mean(tau.true),3), ' our estimate: ', round(tau_est,3 )))
  
  tau_est_list[[i]]<-tau_est
  treat_list[[i]]<- mean(tau.true)
  print(i)
  print(cor(mBG_out[, 'tau'], tau.true))
  #plot( tau.true, mBG_out[, 'tau'], main='shark fin q=0.25 true and same f(u)')
  #text(0.07,.5, paste0('cor=',round(cor(mBG_out[, 'tau'], tau.true),3)))
  #abline(0,1, col='red')
  saveframe=data.frame(mBG_out, B1.true, B0.true, tau.true, E_val)
  #setwd('/home/dpapakos/nonlin_sim_incorrect_holdcsv2/')
  setwd('/home/spange/Dropbox/audit_paper/audit_paper_aoas')
  write.csv(saveframe, paste0('sharkfin_q=_', q, 'sigma=_',sig, '.csv'))
  # write.csv(saveframe,paste0('correct_f_', i, '_sharkfin_',q,'incorrect_f_sharkfin_',qfalse,'.csv'))
}

mean(B1.true/B0.true)
mean(B1/B0)
RMSE <- function(m, o){
  sqrt(mean((m - o)^2))
}
getwd()
saveframe1=read.csv('/home/spange/Dropbox/audit_paper/audit_paper_aoas/norm0sd0.1eval_sim.csv')
saveframe2=read.csv('/home/spange/Dropbox/audit_paper/audit_paper_aoas/norm0sd1eval_sim.csv')
saveframe3=read.csv('/home/spange/Dropbox/audit_paper/audit_paper_aoas/sharkfin_q=_0.25sigma=_0.5.csv')
saveframe4=read.csv('/home/spange/Dropbox/audit_paper/audit_paper_aoas/sharkfin_q=_0.75sigma=_1.25.csv')
par(mfrow=c(2,2))
plot(saveframe4$E_val, saveframe4$B1/saveframe4$B0, xlab='E-value',ylab='Estimated ICRR',  main=expression(paste('Sharkfin q=0.75,', sigma, '=1.25')))
plot(saveframe3$E_val, saveframe3$B1/saveframe3$B0, xlab='E-value',ylab='Estimated ICRR',  main=expression(paste('Sharkfin q=0.25,', sigma, '=0.5')))
plot(saveframe2$E_val, saveframe2$B1/saveframe2$B0, xlab='E-value',ylab='Estimated ICRR',  main='f(u) ~ N(0, 1)')
plot(saveframe1$E_val, saveframe1$B1/saveframe1$B0, xlab='E-value',ylab='Estimated ICRR',  main='f(u) ~ N(0, 0.1)')

dev.off()
getwd()
setwd('/home/dpapakos/nonlin_sim_incorrect_codes/nonlin_sim_constrained')

file_list<-c()
list.files()
library(gtools)
new_files<-mixedsort(sort(list.files()))

for (i in 1:length(new_files)){
  file=new_files[[i]]
  file_list[[i]]<-data.table::fread(file)
}

new_files
#file_list<-rev(file_list)
cor(file_list[[7]]$tau, file_list[[7]]$tau.true, use="complete.obs")
sapply(1:8, function(i) mean(file_list[[i]]$tau.true))
sapply(1:8, function(i) mean(file_list[[i]]$tau))
sapply(1:8, function(i)cor( file_list[[i]]$tau, file_list[[i]]$tau.true, use="complete.obs"))
sapply(1:8, function(i)plot( file_list[[i]]$tau, file_list[[i]]$tau.true, use="complete.obs"))
####the ratios we wanna look at 
sapply(1:8, function(i) RMSE(file_list[[i]]$B1.true/file_list[[i]]$B0.true,file_list[[i]]$B1/file_list[[i]]$B0))

sapply(1:8, function(i)RMSE( file_list[[i]]$tau, file_list[[i]]$tau.true))
par(mfrow=c(1,1))
plot(file_list[[20]]$tau.true, file_list[[20]]$tau, main='q true, lower variance',
     xlab='true tau', ylab='our guess')
abline(0,1 , col='red')
text(0.07,.5, paste0('cor=',round(cor(file_list[[20]][, 'tau'], file_list[[20]]$tau.true),3)))

