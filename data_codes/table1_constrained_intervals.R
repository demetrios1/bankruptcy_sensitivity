setwd("/home/dpapakos/moderating_variables/")

library(dplyr)
library(foreign)
library(tidyverse)
library(foreach)
library(bcf)
library(fastbart)
library(dbarts)
library(foreach)
library(doParallel)

rm(list = ls())
# load data

data      = read.dta( "AuditorFinal20160419.dta" )
#### updated data ####
library(haven)
data      = as.data.frame(read_dta( "auditor_with_stratio.dta"))
colnames(data)
set.seed(12296)

# define lgc

data[,"lgc"] = NA
data[!is.na( data[,"lgc_1"] ),"lgc"] = 1
data[!is.na( data[,"lgc_0"] ),"lgc"] = 0

# create auditor id

data[,"ey"] = data[,"auditor"] == "ey"
data[,"dt"] = data[,"auditor"] == "dt"
data[,"pwc"] = data[,"auditor"] == "pwc"
data[,"kpmg"] = data[,"auditor"] == "kpmg"
data[,"bdo"] = data[,"auditor"] == "bdo"
data[,"gt"] = data[,"auditor"] == "gt"

data      = data[!is.na( data[,"downgrade"] ), ]


#############################################################################

### select sample

sa        = "lgc_0"
data_sa = data[ !is.na( data[,sa] ), ]     # drop missing values of sample indicators
data_sa = data_sa[ data_sa[,sa] == 1 & data_sa$neg_ni==1, ]

new_vars = c('interest', 'stratio', 'sumlogret', 'retvol', 'year1', 
             'year2', 'year3', 'year4', 'year5', 'year6', 'year7', 
             'year8', 'year9','year10', 'year11', 
             'year12', 'year13', 'year14', 'year15')
var0      = c( "logassets_r", "lev_r", "investments_r", "cash_r", 
               "roa_r", "logprice_r", "Intangible_r", "RD_r", 
               "RDmissing", "fracnonfees_r", "feemissing", "NoRate",
               "RateC", "numyears", "downgrade", new_vars)#, "ey", "dt", "kpmg", "pwc", "gt", "bdo" )
vargc   = c( "going_concern", var0 )
vars      = c( "logassets_r", "lev_r", "investments_r", "cash_r", 
               "roa_r", "logprice_r", "Intangible_r", "RD_r", 
               "RDmissing", "fracnonfees_r", "feemissing", "NoRate", 
               "RateC", "numyears", "downgrade",new_vars)#, "ey", "dt", "kpmg", "pwc", "gt", "bdo" )
#subset of sample
data_sasub<-data_sa
data_sasub <- (data_sasub[complete.cases(data_sasub[, vars]), ])
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
BARTpred=function(df, treat=going_concern, Outcome=bankrptobs,vars){
  
  covtreat=df%>%filter(going_concern==1)
  covcontrol=df%>%filter(going_concern==0)
  covtreat$treat=as.factor(covtreat$going_concern)
  
  #return(df[vars])
  #return(as.factor(covtreat$outcome))
  #case 1, the treatment
  #bart1=bart(covtreat[vars],as.factor(covtreat$B), df[vars],ndpost = 5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  
  #case 2 control
  #bart2=bart(covcontrol[vars],as.factor(covcontrol$B),df[vars],ndpost =5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  
  #case 3 propensity
  bart3=bart(df[,vars], as.factor(df$going_concern),df[,vars],ndpost = 2000, nskip = 2000,ntree=100,
             verbose=F,usequants = TRUE,numcut = 1000)
  
  
  xtest=df[vars]
  xtraincontrol=covcontrol[vars]
  xtraintreat=covtreat[vars]
  ytrain0    = as.factor( covcontrol$bankrptobs )
  ytrain1    = as.factor( covtreat$bankrptobs )
  
  # mono fits
  
  bart_mono = monotone_bart(y = as.numeric(c(ytrain1, ytrain0)==1),
                            z = 1-c(rep(1, length(ytrain1)), rep(0, length(ytrain0))),
                            x = rbind(xtraintreat, xtraincontrol),
                            xpred = xtest, nskip = 2000, ndpost = 2000,m=100)
  
  
  
  #bart3 = bart(df[vars], as.factor(df$G),df[vars], ndpost = 5000, nskip = 2000, ntree=100, usequants=T, numcuts=1000, verbose=F)
  
  
  
  #use this method for prediction on binary
  #pred1 = colMeans(pnorm(bart1$yhat.test)) 
  #pred2=colMeans(pnorm(bart2$yhat.test))
  pred1= 1-bart_mono$pr0#1-colMeans(bart_mono$pr0)
  
  pred2= 1-bart_mono$pr1#1-colMeans(bart_mono$pr1)
  pred3=pnorm(bart3$yhat.test)#colMeans(pnorm(bart3$yhat.test))
  
  
  return(list(data.frame(pred1), data.frame(pred2), data.frame(pred3)))
}
var0<-vars
N=dim(data_sasub)[1]
BARTpred=function(df, treat=going_concern, Outcome=bankrptobs,vars){
  
  covtreat=df%>%filter(going_concern==1)
  covcontrol=df%>%filter(going_concern==0)
  covtreat$treat=as.factor(covtreat$going_concern)
  
  #return(df[vars])
  #return(as.factor(covtreat$outcome))
  #case 1, the treatment
  #bart1=bart(covtreat[vars],as.factor(covtreat$B), df[vars],ndpost = 5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  
  #case 2 control
  #bart2=bart(covcontrol[vars],as.factor(covcontrol$B),df[vars],ndpost =5000, nskip = 2000,ntree=100,verbose=F,usequants = TRUE,numcut = 1000)
  
  #case 3 propensity
  bart3=bart(df[,vars], as.factor(df$going_concern),df[,vars],ndpost = 2000, nskip = 2000,ntree=100,
             verbose=F,usequants = TRUE,numcut = 1000)
  
  
  xtest=df[vars]
  xtraincontrol=covcontrol[vars]
  xtraintreat=covtreat[vars]
  ytrain0    = as.factor( covcontrol$bankrptobs )
  ytrain1    = as.factor( covtreat$bankrptobs )
  
  # mono fits
  
  bart_mono = monotone_bart(y = as.numeric(c(ytrain1, ytrain0)==1),
                            z = 1-c(rep(1, length(ytrain1)), rep(0, length(ytrain0))),
                            x = rbind(xtraintreat, xtraincontrol),
                            xpred = xtest, nskip = 2000, ndpost = 2000,m=100)
  
  
  
  #bart3 = bart(df[vars], as.factor(df$G),df[vars], ndpost = 5000, nskip = 2000, ntree=100, usequants=T, numcuts=1000, verbose=F)
  
  
  
  #use this method for prediction on binary
  #pred1 = colMeans(pnorm(bart1$yhat.test)) 
  #pred2=colMeans(pnorm(bart2$yhat.test))
  pred1= 1-bart_mono$pr0#1-colMeans(bart_mono$pr0)
  
  pred2= 1-bart_mono$pr1#1-colMeans(bart_mono$pr1)
  pred3=pnorm(bart3$yhat.test)#colMeans(pnorm(bart3$yhat.test))
  
  
  return(list(data.frame(pred1), data.frame(pred2), data.frame(pred3)))
}
#intframe=BARTpred(data_sasub, treat=going_concern, Outcome=bankrptobs,vars=var0)
#intframe=BARTpred(data_sasub, treat=going_concern, Outcome=bankrptobs,vars=var0)
#write.csv(intframe, "/home/dpapakos/cred_intervals/mod_var_mbart_newvar.csv")
library(readr)
intframe_temp=data.table::fread("/home/dpapakos/cred_intervals/mod_var_mbart_newvar.csv")
intframe<-c()
divs = floor(dim(intframe_temp)[2]/3)
divs
intframe[[1]] <- intframe_temp[,2:(divs+1)]
intframe[[2]] <- intframe_temp[,(divs+2):(2*divs +1)]
intframe[[3]] <- intframe_temp[,(2*divs+2):(3*divs +1)]
rm(intframe_temp)



fill4=matrix(NA,nrow=2000, ncol=dim(intframe[[1]])[2])
library(foreach)
library(doParallel)
####helper functions####

np = detectCores()-1
#np=16
#cl = makeCluster(np)
#rm(cl)
n_cores <- detectCores() - 1
#n_cores=20
n_cores
registerDoParallel(cores = n_cores)
#registerDoMC(cores = n_cores)
parral_arrange <- split(1:2000, 1:n_cores)
parral_arrange[[1]]
rm(data)
rm(data_sa)
rm(data_sasub)
#mBG_out<-c()
n_cores <- detectCores() - 1
#n_cores=20
n_cores=16
#cl=makeCluster(n_cores)
shark_var_calc <- function(q, sig){
  a = 1-(((dnorm(0))/(1-pnorm(0)))^2)
  b = ((dnorm(0))/(1-pnorm(0)))
  
  sigminus=a*((sig*(1-q))/q)^2
  sigplus=a*sig^2
  mplus=b*((sig*(1-q))/q)
  mminus=-b*sig
  
  meq=q*mplus+(1-q)*mminus
  
  sd_est <- sqrt(q*(sigplus+mplus^2)+(1-q)*(sigminus+mminus^2)-meq^2)
  return(sd_est)
}
#for our mixture models, functiion to calculate variance
var_fun_mix <- function(w, means, sig){
  meq=sum(w*means)
  vareq=sum(w*(sig^2+means^2))-meq^2
  return(sqrt(vareq))
}

w1=c(0.01, 0.94, 0.05)
w2 = c(0.05, 0.90, 0.05)
mass = c(-2, 0, 2)
sigma_list =c(0.05, 0.05, 0.05)
assym_var = var_fun_mix(w1, mass, sigma_list)
sym_var = var_fun_mix(w2, mass, sigma_list)
shark_var = shark_var_calc(q,sigma_list)
functionlist <- c(function(u){dnorm(u, mean=0, sd=0.1)}, function(u){dnorm(u,mean=0,sd=0.5)},
                  function(u){dnorm(u,mean=0,sd=1)}, function(u){
                    .01*dnorm(u, mean=-2, sd=0.05)+0.94*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
                  },function(u){
                    .05*dnorm(u, mean=-2, sd=0.05)+0.90*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
                  }, function(z){
                    val = rep(NA, length(z))
                    val[z < 0] = 2*0.25*dnorm(z[z<0], sd = 0.5)
                    val[z>0] = 2*(1-0.25)*dnorm(z[z>0],0,sd = 0.5*(1-0.25)/0.25)
                    return(val)
                  },  function(z){
                    val = rep(NA, length(z))
                    val[z < 0] = 2*0.75*dnorm(z[z<0], sd = 1.25)
                    val[z>0] = 2*(1-0.75)*dnorm(z[z>0],0,sd = 1.25*(1-0.75)/0.75)
                    return(val)
                    
                  }
)
filenames_B0 <- c("mean0sd0.1_credint_B0.csv", "mean0sd0.5_credint_B0.csv",
                  "mean0sd1_credint_B0.csv","rightbumpsigma0.48_B0.csv", 
                  "symmetric90percmid_B0.csv", "sharkfinq0.25sig0.5_B0.csv",
                  "sharkfinq0.75sig1.25_B0.csv"
)
filenames_B1 <- c("mean0sd0.1_credint_B1.csv", "mean0sd0.5_credint_B1.csv",
                  "mean0sd1_credint_B1.csv","rightbumpsigma0.48_B1.csv", 
                  "symmetric90percmid_B1.csv",  "sharkfinq0.25sig0.5_B1.csv", 
                  "sharkfinq0.75sig1.25_B1.csv" 
)

registerDoParallel(cores = n_cores)
ptm = proc.time()
for (k in 5:length(functionlist)){
  set.seed(12296)
  post_num = 500
  post_list=sample.int(2000, post_num, replace=F)
  f=functionlist[k][[1]]
  comb_frame=foreach(i= 1:post_num, .combine=rbind) %dopar% {
    
    
    #merge_frame=data.frame(treatpred=t(intframe[[1]][parral_arrange[[i]],]), 
    #                       notreatpred=t(intframe[[2]][parral_arrange[[i]],]),
    #                          propensity=t(intframe[[3]][parral_arrange[[i]],]))
    merge_frame=data.frame(treatpred=t(intframe[[1]][post_list[[i]],]), 
                           notreatpred=t(intframe[[2]][post_list[[i]],]),
                           propensity=t(intframe[[3]][post_list[[i]],]))
    dim(merge_frame)
    #f=function(u){
    #  dnorm(u, mean=0, sd=0.1)
    #}
    
    
    
    # f = function(z){
    #   val = rep(NA, length(z))
    #    val[z < 0] = 2*q*dnorm(z[z<0], sd = sig)
    #    val[z>0] = 2*(1-q)*dnorm(z[z>0],0,sd = sig*(1-q)/q)
    #   return(val)
    # }
    #  }
    # f=function(u){
    #    .01*dnorm(u, mean=-2, sd=0.05)+0.94*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
    #  }
    #f=function(u){
    #  .01*dnorm(u, mean=-2, sd=0.05)+0.98*dnorm(u,mean=0,sd=0.05)+.01*dnorm(u,mean=2,sd=0.05)
    #  }
    colnames(merge_frame)=c('treatpred', 'notreatpred','propensity')
    # expoutcomesfun	   =cbind(data.frame(  treatpred = pred1, notreatpred=pred2), propensity=pred3 )
    #  expoutcomesfun
    # expoutcomesfun2    = rbind( NULL, expoutcomesfun )
    outcomenew=merge_frame[,c('propensity','notreatpred','treatpred')]
    
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
    #row.names(outcomenew)<-NULL
    newoutcome =as.matrix(outcomenew)
    
    dim(newoutcome)
    
    #ptm = proc.time()
    
    #newoutcome     = newoutcome[,c( "ProbY1D1", "ProbY1D0", "ProbY0D1" )] 
    print( sd )
    
    
    #mBG_out[parral_arrange[[i]]]
    #fill4[parral_arrange[[i]],]= 
    
    #fill4[i,]=foreach( m = 1:dim(newoutcome)[1],  .combine=rbind )%dopar%{
    
    mBG_out<-foreach( m = 1:dim(newoutcome)[1],  .combine=rbind ) %do% {
      
      # Probability to fit.
      vProbBG     = newoutcome[m,c( "ProbY1D1", "ProbY1D0", "ProbY0D1" )] 
      
      
      # Starting value
      
      start0 = c( newoutcome[m,"ProbY1D1"], newoutcome[m,"ProbY1D0"],newoutcome[m,"ProbY0D1"] )
      
      names(start0) = c("y1","y0","d")
      
      
      #better global definition.
      #this is the distance we wanna minimize
      #here, we note that b1>b0, so we set b1=exp(b1)+b0 to ensure this
      #convert back later
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
      treatment_Compute_constraint = function( vb){
        integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u) - pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value
        
      }
      
      #calculate the counter factuals
      counterfac_Compute_constraint = function(vb){
        y1 =(integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
        y0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
        
        return(c(y1, y0))
        
      }
      #?optim
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
        
        
        out        = counter#tau#c( vBG, tau,vBG_out$convergence, vBG_out$value,counter) 
      }
      return( out)
      
    }
    
    
    #stopImplicitCluster()
    #colnames(mBG_out[parral_arrange[[i]]]) = 'tau'#c( "y1", "y0", "d", "tau", "convergence", "fnvalue" ,'B1','B0')
    
    
    #fill4[parral_arrange[[i]], ]<-mBG_out[[parral_arrange[[i]]]][, 'tau']
    #saveframe=data.frame(mBG_out)
    #rm(saveframe)
    print(post_list[i])
    return(c(B1=unlist(mBG_out[, 1]),B0=mBG_out[,2]))
    
    # print(parral_arrange[[i]])
  }
  
  
  
  n_obs=dim(intframe[[1]])[2]
  
  B1=comb_frame[,1:n_obs]
  dim(comb_frame)
  
  B0=comb_frame[, (n_obs+1):(2*n_obs)]
  print( proc.time() - ptm )
  
  setwd('/home/dpapakos/cred_intervals/credinterval_constrained_integral_newvar/')
  #write.csv(fill, 'mean0sd0.5_largestATEsub.csv')
  write.csv(B1, filenames_B1[k],row.names = FALSE)
  write.csv(B0, filenames_B0[k],row.names = FALSE)
  print(k)
}
setwd('/home/dpapakos/cred_intervals/credinterval_constrained_integral_newvar/')

#dim(fill4)
fill4=data.table::fread('mean0sd1_credint_B1.csv')[-1]/data.table::fread('mean0sd1_credint_B0.csv')[-1]
fill3=data.table::fread('mean0sd0.5_credint_B1.csv')[-1]/data.table::fread('mean0sd0.5_credint_B0.csv')[-1]
fill2=data.table::fread('mean0sd0.1_credint_B1.csv')[-1]/data.table::fread('mean0sd0.1_credint_B0.csv')[-1]
fillshark1=data.table::fread('sharkfinq0.25sig0.5_B1.csv')[-1]/data.table::fread('sharkfinq0.25sig0.5_B0.csv')[-1]
fillshark2=data.table::fread('sharkfinq0.75sig1.25_B1.csv')[-1]/data.table::fread('sharkfinq0.75sig1.25_B0.csv')[-1]
fillrightbump=data.table::fread('rightbumpsigma0.48_B1.csv')[-1]/data.table::fread('rightbumpsigma0.48_B0.csv')[-1]
#fillsym98=data.table::fread('symmetric98percmid_B1.csv')[-1]-data.table::fread('symmetric98percmid_B0.csv')[-1]
fillsym90=data.table::fread('symmetric90percmid_B1.csv')[-1]/data.table::fread('symmetric90percmid_B0.csv')[-1]
gammasym90=rowMeans(fillsym90)
#gammasym98=rowMeans(fillsym98)
gammarightbump=rowMeans(fillrightbump)

#gamma=rowMeans(fill)
gammashark1=rowMeans(fillshark1)
gammashark2=rowMeans(fillshark2)

gamma2=rowMeans(fill2)
gamma3=rowMeans(fill3)
gamma4=rowMeans(fill4)


fill4_rd=data.table::fread('mean0sd1_credint_B1.csv')[-1]-data.table::fread('mean0sd1_credint_B0.csv')[-1]
fill3_rd=data.table::fread('mean0sd0.5_credint_B1.csv')[-1]-data.table::fread('mean0sd0.5_credint_B0.csv')[-1]
fill2_rd=data.table::fread('mean0sd0.1_credint_B1.csv')[-1]-data.table::fread('mean0sd0.1_credint_B0.csv')[-1]
fillshark1_rd=data.table::fread('sharkfinq0.25sig0.5_B1.csv')[-1]-data.table::fread('sharkfinq0.25sig0.5_B0.csv')[-1]
fillshark2_rd=data.table::fread('sharkfinq0.75sig1.25_B1.csv')[-1]-data.table::fread('sharkfinq0.75sig1.25_B0.csv')[-1]
fillrightbump_rd=data.table::fread('rightbumpsigma0.48_B1.csv')[-1]-data.table::fread('rightbumpsigma0.48_B0.csv')[-1]
#fillsym98=data.table::fread('symmetric98percmid_B1.csv')[-1]-data.table::fread('symmetric98percmid_B0.csv')[-1]
fillsym90_rd=data.table::fread('symmetric90percmid_B1.csv')[-1]-data.table::fread('symmetric90percmid_B0.csv')[-1]
gammasym90_rd=rowMeans(fillsym90_rd)
#gammasym98=rowMeans(fillsym98)
gammarightbump_rd=rowMeans(fillrightbump_rd)

#gamma=rowMeans(fill)
gammashark1_rd=rowMeans(fillshark1_rd)
gammashark2_rd=rowMeans(fillshark2_rd)

gamma2_rd=rowMeans(fill2_rd)
gamma3_rd=rowMeans(fill3_rd)
gamma4_rd=rowMeans(fill4_rd)

#mean(fill4[1, ], na.rm=T)
#sum(fill4[2000,])/ncol(fill4)
mean(gamma3, na.rm=T) #N(0, .5)
mean(gamma4, na.rm=T) #(N(0,1))
mean(gamma2, na.rm=T) #N(0, .1)

df2=data.frame(run=c('N(0,.1)','N(0, 0.5)', 'N(0,1)', 'shark q=0.25, s=0.5; v=1.11', 
                     'shark q=0.75, s=1.25; v=0.77', 
                     '90% peak', 'rightbump'), 
               mean=c(mean(gamma2),mean(gamma3), mean(gamma4), 
                      mean(gammashark1), mean(gammashark2),
                      mean(gammasym90),
                      mean(gammarightbump)), 
               quantiles_low=c(quantile(gamma2, c(0.025, .5, .975))[1],
                               quantile(gamma3, c(0.025, .5, .975))[1], 
                               quantile(gamma4, c(0.025, .5, .975))[1], 
                               quantile(gammashark1, c(0.025, .5, .975))[1], 
                               quantile(gammashark2, c(0.025, .5, .975))[1], 
                               quantile(gammarightbump, c(0.025, .5, .975))[1], 
                               # quantile(gammasym98, c(0.025, .5, .975))[1],
                               quantile(gammasym90, c(0.025, .5, .975))[1]), 
               quantile_mid=c(quantile(gamma2, c(0.025, .5, .975))[2],
                              quantile(gamma3, c(0.025, .5, .975))[2], 
                              quantile(gamma4, c(0.025, .5, .975))[2], 
                              quantile(gammashark1, c(0.025, .5, .975))[2], 
                              quantile(gammashark2, c(0.025, .5, .975))[2], 
                              quantile(gammasym90, c(0.025, .5, .975))[2], 
                              #quantile(gammasym98, c(0.025, .5, .975))[2],
                              quantile(gammarightbump, c(0.025, .5, .975))[2]), 
               quantile_high=c(quantile(gamma2, c(0.025, .5, .975))[3],
                               quantile(gamma3, c(0.025, .5, .975))[3], 
                               quantile(gamma4, c(0.025, .5, .975))[3], 
                               quantile(gammashark1, c(0.025, .5, .975))[3], 
                               quantile(gammashark2, c(0.025, .5, .975))[3], 
                               quantile(gammasym90, c(0.025, .5, .975))[3], 
                               #quantile(gammasym98, c(0.025, .5, .975))[3],
                               quantile(gammarightbump, c(0.025, .5, .975))[3]), 
               mean_rd=c(mean(gamma2_rd),mean(gamma3_rd), mean(gamma4_rd), 
                         mean(gammashark1_rd), mean(gammashark2_rd),
                         mean(gammasym90_rd),
                         mean(gammarightbump_rd)), 
               quantiles_low_rd=c(quantile(gamma2_rd, c(0.025, .5, .975))[1],
                                  quantile(gamma3_rd, c(0.025, .5, .975))[1], 
                                  quantile(gamma4_rd, c(0.025, .5, .975))[1], 
                                  quantile(gammashark1_rd, c(0.025, .5, .975))[1], 
                                  quantile(gammashark2_rd, c(0.025, .5, .975))[1], 
                                  quantile(gammarightbump_rd, c(0.025, .5, .975))[1], 
                                  # quantile(gammasym98, c(0.025, .5, .975))[1],
                                  quantile(gammasym90_rd, c(0.025, .5, .975))[1]), 
               quantile_mid_rd=c(quantile(gamma2_rd, c(0.025, .5, .975))[2],
                                 quantile(gamma3_rd, c(0.025, .5, .975))[2], 
                                 quantile(gamma4_rd, c(0.025, .5, .975))[2], 
                                 quantile(gammashark1_rd, c(0.025, .5, .975))[2], 
                                 quantile(gammashark2_rd, c(0.025, .5, .975))[2], 
                                 quantile(gammasym90_rd, c(0.025, .5, .975))[2], 
                                 #quantile(gammasym98, c(0.025, .5, .975))[2],
                                 quantile(gammarightbump_rd, c(0.025, .5, .975))[2]), 
               quantile_high_rd=c(quantile(gamma2_rd, c(0.025, .5, .975))[3],
                                  quantile(gamma3_rd, c(0.025, .5, .975))[3], 
                                  quantile(gamma4_rd, c(0.025, .5, .975))[3], 
                                  quantile(gammashark1_rd, c(0.025, .5, .975))[3], 
                                  quantile(gammashark2_rd, c(0.025, .5, .975))[3], 
                                  quantile(gammasym90_rd, c(0.025, .5, .975))[3], 
                                  #quantile(gammasym98, c(0.025, .5, .975))[3],
                                  quantile(gammarightbump_rd, c(0.025, .5, .975))[3]))
xtable::xtable(df2)
#View(df2)



#write.csv(df2, 'constrained_int_summary.csv')
#plot(density(gamma4), col='dodgerblue4', xlim=c(-.02, 0.02),main='', ylim=c(0, 180), lwd=2, xlab='Posterior ATE,f(u)~N(0,1)')
#lines(density(gamma2), col='firebrick4', lwd=2)
#lines(density(gamma3), col='darkorchid4', lwd=2)

#lines(density(gamma4), col='black', lwd=2)
#legend("topright", legend=c(expression(paste("f" [1],"(u)")), expression(paste("f" [2],"(u)")), 
#                            expression(paste("f" [3],"(u)")), expression(paste("f" [4],"(u)"))),
#       col=c("dodgerblue4", "firebrick4", 'darkorchid4', 'black'), lty=1, cex=1.1)

quantile(gamma4, c(0.025, .5, .975))[1]
#quantile(gammashark1, c(0.025, .5, .975))
#quantile(gammashark2, c(0.025, .5, .975))
#quantile(gamma3, c(0.025, .5, .975))
#quantile(gamma2, c(0.025, .5, .975))
#quantile(gammarightbump, c(0.025, .5, .975))
#quantile(gammasym98, c(0.025, .5, .975))
quantile(gammasym90, c(0.025, .5, .975))