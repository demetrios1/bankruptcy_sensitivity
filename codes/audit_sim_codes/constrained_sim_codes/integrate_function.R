#' A Cat Function
#'
#' This function does the integration from the paper
#' Here, we pass the mean of the posterior estimates of observed joint probabilities
#' rather than the full posterior estimates.  
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' integrate_function()
library(dplyr)
library(foreign)
library(tidyverse)
library(foreach)
library(bcf)
library(fastbart)
library(dbarts)
library(foreach)
library(doParallel)

integrate_function=function(intframe, constraint=T, f=f, n_cores=4){
  #n_cores <- detectCores() - 1
  n_cores=n_cores
  registerDoParallel(cores = n_cores)
  set.seed(12296)
  mBG_out<-c()
  #sdlist=c(0.1, 0.5, 1)
  q=0.25
  sig=.5
#for (m in sdlist){#:length(gamma_list)){
  #better global definition.
  optimdist=function(vals){
    b1 = vals[1]
    b0 = vals[2]
    g = vals[3]
    a=integrate( function(u) pnorm(b1+u)*pnorm(g+u)*f(u), lower = -Inf, upper = Inf )$value
    b=integrate(function(u) pnorm(b0+u)*(1-pnorm(g+u))*f(u), lower = -Inf, upper = Inf )$value
    c=integrate(  function(u) (1-pnorm(b1+u))*pnorm(g+u)*f(u), lower = -Inf, upper = Inf )$value 
    #return(c(a,b,c))
    #vProbBG is lhs, global variable
    return(sum((qnorm(c(a,b,c))-qnorm(vProbBG))^2))
  }
  optimdist_constraint=
    function(vals){
      b1 = vals[1]
      b0 = vals[2]
      g = vals[3]
      #replace b1 with b1+b0
      a=integrate( function(u) pnorm((exp(b1)+b0)+u)*pnorm(g+u)*f(u), lower = -Inf, upper = Inf )$value
      b=integrate(function(u) pnorm(b0+u)*(1-pnorm(g+u))*f(u), lower = -Inf, upper = Inf )$value
      c=integrate(  function(u) (1-pnorm((exp(b1)+b0)+u))*pnorm(g+u)*f(u), lower = -Inf, upper = Inf )$value 
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
    b1 =(integrate( function(u) (pnorm(vb[1]+u))*f(u), lower = -Inf, upper = Inf )$value)
    b0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
    
    return(c(b1, b0))
    
  }
  
  counterfac_Compute_constraint = function(vb){
    b1 =(integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
    b0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
    
    return(c(b1, b0))
    
  }
  newoutcome=as.matrix(intframe)
  
  ptm = proc.time()
  
  mBG_out = foreach( i = 1:dim(newoutcome)[1],  .combine=rbind )%dopar%{
   
 
    # Probability to fit.
    vProbBG     = newoutcome[i,c( "ProbB1G1", "ProbB1G0", "ProbB0G1" )] 
    
    
    # Starting value
    
    start0 = c( newoutcome[i,"ProbB1G1"], newoutcome[i,"ProbB1G0"],newoutcome[i,"ProbB0G1"] )
    
    names(start0) = c("b1","b0","g")
    
    
    
  
    if(constraint==T){
    vBG_out     = try(optim( qnorm(start0), optimdist_constraint,  control = list(maxit = 500) 
                      )
                      , silent=TRUE)
    if ( class( vBG_out ) == "try-error" ){
      out        = as.numeric(newoutcome[i,c("indexes")])
      out        = c( out, rep(NA,5) ) 
    } else {
      out        = as.numeric(newoutcome[i,c("indexes")] )
      vBG        = vBG_out$par 
     # tau      = treatment_Compute( vb = vBG[c(1,2)])
      
    #  counter=counterfac_Compute(vb = vBG[c(1,2)])
      
      tau      = treatment_Compute_constraint( vb = vBG[c(1,2)])
      
      counter=counterfac_Compute_constraint(vb = vBG[c(1,2)])
      
      out        = c( vBG, tau,vBG_out$convergence, vBG_out$value,counter) 
    }
    return( out)
}else{
vBG_out     = try(optim( qnorm(start0), optimdist,  control = list(maxit = 500) 
                      )
                      , silent=TRUE)
    if ( class( vBG_out ) == "try-error" ){
      out        = as.numeric(newoutcome[i,c("indexes")])
      out        = c( out, rep(NA,5) ) 
    } else {
      out        = as.numeric(newoutcome[i,c("indexes")] )
      vBG        = vBG_out$par 
      tau      = treatment_Compute( vb = vBG[c(1,2)])
      
     counter=counterfac_Compute(vb = vBG[c(1,2)])
      
     # tau      = treatment_Compute_constraint( vb = vBG[c(1,2)])
      
      #counter=counterfac_Compute_constraint(vb = vBG[c(1,2)])
      
      out        = c( vBG, tau,vBG_out$convergence, vBG_out$value,counter) 
    }
    return( out)
}
  }

  print( proc.time() - ptm )
  
  colnames(mBG_out) = c( "b1", "b0", "d", "tau", "convergence", "fnvalue" ,'B1','B0')
  setwd("/home/dpapakos/sensitivity_analysis/")


  return(mBG_out)
}




#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()





set.seed(12296)


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
var0<-vars
N=dim(data_sasub)[1]

BARTpred=function(df, treat=G, Outcome=B,vars, mono=T, nd_post=20, n_skip=20){
  df=data_sasub
  treat='going_concern'
  Outcome='bankrptobs'
  vars=var0
  
  covtreat=df%>%filter(going_concern==1)
  covcontrol=df%>%filter(going_concern==0)
  covtreat$going_concern=as.factor(covtreat$going_concern)
 covtreat=df[df[treat]==1,]
 covcontrol=df[df[treat]==0,]
 covtreat[treat]=as.factor(covtreat[treat])
  #return(df[vars])
  #return(as.factor(covtreat$outcome))
  #case 1, the treatment
  if (mono==F){
    bart1=bart(covtreat[vars],as.factor(unlist(covtreat[Outcome])), df[vars],ndpost = nd_post, nskip = n_skip,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
    
    #case 2 control
    bart2=bart(covcontrol[vars],as.factor(unlist(covcontrol[Outcome])),df[vars],ndpost =nd_post, nskip = n_skip,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
    
    #case 3 propensity
    bart3=bart(df[vars], as.factor(unlist(df[treat])),df[vars],ndpost = nd_post, nskip = n_skip,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
    pred1 = colMeans(pnorm(bart1$yhat.test)) 
    pred2=colMeans(pnorm(bart2$yhat.test))
    pred3=colMeans(pnorm(bart3$yhat.test))
  }
  else if(mono==T){
    xtest=df[vars]
    xtraincontrol=covcontrol[vars]
    xtraintreat=covtreat[vars]
    ytrain0    = as.factor( unlist(covcontrol[Outcome] ))
    ytrain1    =as.factor( unlist(covtreat[Outcome]))

    # mono fits

    bart_mono = monotone_bart(y = as.numeric(c(ytrain1, ytrain0)==1),
                              z = 1-c(rep(1, length(ytrain1)), rep(0, length(ytrain0))),
                              x = rbind(xtraintreat, xtraincontrol),
                              xpred = xtest, nskip = n_skip, ndpost = nd_post,m=100)
    

    bart3 = bart(df[vars], as.factor(unlist(df[treat])),df[vars], ndpost = nd_post, nskip = n_skip, ntree=100, usequants=T,  verbose=F)
   
  
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
  
  outcomenew[,c("prop", "B_G0", "B_G1")] = 0.5*eps + (1-eps)*outcomenew[,c('propensity','notreatpred','treatpred')]

  outcomenew[,"ProbB1G1"]=outcomenew[,"prop"]*outcomenew[,"treatpred"]
  outcomenew[,"ProbB1G0"] = (1-outcomenew[,"prop"])*outcomenew[,"notreatpred"]
  outcomenew[,"ProbB0G1"] = outcomenew[,"prop"]*(1-outcomenew[,"treatpred"])
  
  #for error analysis, not super necessary
  indexes=seq(from=1, to=length(outcomenew$treatpred), by=1)
  outcomenew=as.data.frame(cbind(indexes, outcomenew))
  outcomenew$indexes=as.numeric(outcomenew$indexes)
  row.names(outcomenew)<-NULL
  newoutcome4 =as.matrix(outcomenew)
  return(outcomenew)
}
