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

var0      = c( "logassets_r", "lev_r", "investments_r", "cash_r", "roa_r", "logprice_r", "Intangible_r", "RD_r", "RDmissing", "fracnonfees_r", "feemissing", "NoRate", "RateC", "numyears", "downgrade")#, "ey", "dt", "kpmg", "pwc", "gt", "bdo" )
vargc   = c( "going_concern", var0 )
vars      = c( "logassets_r", "lev_r", "investments_r", "cash_r", "roa_r", "logprice_r", "Intangible_r", "RD_r", "RDmissing", "fracnonfees_r", "feemissing", "NoRate", "RateC", "numyears", "downgrade")#, "ey", "dt", "kpmg", "pwc", "gt", "bdo" )
#subset of sample
data_sasub<-data_sa

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
BARTpred=function(df, treat=G, Outcome=B,vars){
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
  
  set.seed(12296)
  post_list=sample.int(2000, 500, replace=F)
  
  #use this method for prediction on binary
  #pred1 = colMeans(pnorm(bart1$yhat.test)) 
  #pred2=colMeans(pnorm(bart2$yhat.test))
  pred1= 1-colMeans(bart_mono$pr0)
  pred2= 1-colMeans(bart_mono$pr1)
  #pred1= 1-colMeans(bart_mono$pr0[])#dont index at post_list
  #pred2= 1-colMeans(bart_mono$pr1[]) #dont index at post_list
  pred3=colMeans(pnorm(bart3$yhat.test[])) #dont index at post_list
  
  #pred1 = colMeans(pnorm(bart1$yhat.test)) 
  #pred2=colMeans(pnorm(bart2$yhat.test))
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
intframe=BARTpred(data_sasub, treat=going_concern, Outcome=bankrptobs,vars=var0)

#write.csv(intframe, "/home/dpapakos/moderating_variables/mod_var_mbart.csv")
library(readr)
#intframe=read_csv("/home/dpapakos/moderating_variables/mod_var_mbart.csv")

#indices=read_csv('mu0sd0.5largestsub.csv')[,2]



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

#mBG_out<-c()
n_cores <- detectCores() - 1
#n_cores=20
n_cores=10
#cl=makeCluster(n_cores)
registerDoParallel(cores = n_cores)
ptm = proc.time()
set.seed(12296)
mBG_out<-c()
sdlist=c(0.1, 0.5, 1)
q=0.25
sig=.5
for (m in sdlist){#:length(gamma_list)){


  #cl = makeCluster(np)
  #registerDoParallel(cl)
  newoutcome=as.matrix(intframe)
  
  ptm = proc.time()
  
  mBG_out = foreach( i = 1:dim(newoutcome)[1],  .combine=rbind )%dopar%{
    f=function(u){
      .01*dnorm(u, mean=-2, sd=0.05)+0.94*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
    }
    #f=function(u){
    #  .01*dnorm(u, mean=-2, sd=0.05)+0.98*dnorm(u,mean=0,sd=0.05)+.01*dnorm(u,mean=2,sd=0.05)
    #}
    #f=function(u){
  #    dnorm(u, mean=0, sd=m)
  #  }
   
  #  f = function(z){
  #    val = rep(NA, length(z))
  #    val[z < 0] = 2*q*dnorm(z[z<0], sd = sig)
  #    val[z>0] = 2*(1-q)*dnorm(z[z>0],0,sd = sig*(1-q)/q)
  #    return(val)
      
  #  }
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
      
      out        = c( vBG, tau,vBG_out$convergence, vBG_out$value,counter) 
    }
    return( out)
  }
  
  print( proc.time() - ptm )
  
  colnames(mBG_out) = c( "y1", "y0", "d", "tau", "convergence", "fnvalue" ,'B1','B0')
  setwd("/home/dpapakos/sensitivity_analysis/constrained_integration/")

  write.csv( mBG_out, file = paste("inducement_rightbump.csv",sep = ""), row.names = FALSE )
  print(m)

}

setwd("/home/dpapakos/sensitivity_analysis/constrained_integration/")
file_list<-c()
list.files()
library(gtools)
new_files<-mixedsort(sort(list.files()))

for (i in 1:length(new_files)){
  file=new_files[[i]]
  file_list[[i]]<-data.table::fread(file)
}
new_files
RMSE <- function(m, o){
  sqrt(mean((m - o)^2))
}





eps=1e-6
ratio_Compute_constraint = function( vb){
  integrate( function(u) ((pnorm(exp(vb[1])+vb[2]+u))/(pnorm(vb[2]+u)))*f(u), lower = -Inf, upper = Inf )$value
}
tau_ratio<-c()
sd_list=c(0.1, 0.5, 1)
for (i in 1:3){
vb=file_list[[i]][,1:2]
f=function(u){
  dnorm(u, mean=0, sd=sd_list[[i]])
}
tau_ratio[[i]]=foreach( i = 1:dim(file_list[[1]])[1],  .combine=rbind )%dopar%{

  tau_ratio=ratio_Compute_constraint( vb = unlist(vb[i,]))
}
}

mean(tau_ratio[[3]])
  f = function(z){
    val = rep(NA, length(z))
    val[z < 0] = 2*q*dnorm(z[z<0], sd = sig)
    val[z>0] = 2*(1-q)*dnorm(z[z>0],0,sd = sig*(1-q)/q)
    return(val)
  }
q=0.25
sig=0.5

vb=file_list[[7]][, 1:2]
tau_ratio=foreach( i = 1:dim(file_list[[1]])[1],  .combine=rbind )%dopar%{
  
  tau_ratio=ratio_Compute_constraint( vb = unlist(vb[i,]))
}

mean(tau_ratio)

q=0.75
sig=1.25

vb=file_list[[8]][, 1:2]
tau_ratio=foreach( i = 1:dim(file_list[[1]])[1],  .combine=rbind )%dopar%{
  
  tau_ratio=ratio_Compute_constraint( vb = unlist(vb[i,]))
}

mean(tau_ratio)


#rightbump
vb=file_list[[4]][, 1:2]
f=function(u){
  .01*dnorm(u, mean=-2, sd=0.05)+0.94*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
}
tau_ratio=foreach( i = 1:dim(file_list[[1]])[1],  .combine=rbind )%dopar%{
  
  tau_ratio=ratio_Compute_constraint( vb = unlist(vb[i,]))
}

mean(tau_ratio)
#90%
vb=file_list[[5]][, 1:2]
f=function(u){
  .05*dnorm(u, mean=-2, sd=0.05)+0.9*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
}
tau_ratio=foreach( i = 1:dim(file_list[[1]])[1],  .combine=rbind )%dopar%{
  
  tau_ratio=ratio_Compute_constraint( vb = unlist(vb[i,]))
}

mean(tau_ratio)

#98%
vb=file_list[[6]][, 1:2]
f=function(u){
  .01*dnorm(u, mean=-2, sd=0.05)+0.98*dnorm(u,mean=0,sd=0.05)+.01*dnorm(u,mean=2,sd=0.05)
}
tau_ratio=foreach( i = 1:dim(file_list[[1]])[1],  .combine=rbind )%dopar%{
  
  tau_ratio=ratio_Compute_constraint( vb = unlist(vb[i,]))
}

mean(tau_ratio)


new_files
sapply(1:length(file_list), function(i) mean(file_list[[i]]$tau))
sapply(1:length(file_list), function(i) mean(file_list[[i]]$B1))
#the ratio estimand
sapply(1:length(file_list), function(i) mean(file_list[[i]]$B1/file_list[[i]]$B0))
#sapply(1:length(file_list), function(i)cor( file_list[[i]]$tau, file_list[[i]]$tau.true, use="complete.obs"))
#file_list[[1]]$tau
#sapply(1:(length(file_list)), function(i)RMSE( file_list[[i]]$tau, file_list[[i]]$tau.true))
new_files
theme_set(theme_minimal(base_size = 16))
gammarightmore2<-file_list[[4]]%>%
  ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),color='white',fill='black')+
  ggtitle(expression(paste(mu,'=-2,0,2,',  sigma, '=0.05'))) +
  ylab('Density')+xlab('ITE')+
  theme(plot.title = element_text(hjust = 0.5))+xlim(-0.01, .14)#+annotate("text", x = .15, y = .089, label='mean inducement:',size=6)+annotate("text", x = .265, y = .089, label=100*round(rightmore2mean$meangam,4),size=6)+annotate("text", x = .16, y = .159, label='mean bankruptcy prob:',size=6)+annotate("text", x = .265, y = .159, label=100*round(rightmore2mean$meanB1,4),size=6)+xlim(-.06,.29)

gammarightmore2

setwd("/home/dpapakos/sensitivity_analysis/constraint_plots/")
ggsave(filename="gammarightmore2_constrained.pdf", height=7, width=7, plot=gammarightmore2)
multi1perc<-file_list[[6]]%>%
  ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),color='white',fill='black')+
  ggtitle('98 percent around main peak')+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab('Density')+xlab('ITE')+xlim(-0.02, .3)#+annotate("text", x = .15, y = .089, label='mean inducement:',size=6)+annotate("text", x = .285, y = .089, label=100*round(multi1percmean$meangam,4),size=6)+annotate("text", x = .16, y = .139, label='mean bankruptcy prob:',size=6)+annotate("text", x = .285, y = .139, label=100*round(multi1percmean$meanB1,4),size=6)+xlim(-.02,.4)
multi1perc
setwd("/home/dpapakos/sensitivity_analysis/constraint_plots/")
ggsave(filename="symmetric1percbump_constrained.pdf",height=7, width=7, plot=multi1perc)

new_files
multi5perc<-file_list[[8]]%>%
  ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),color='white',fill='black')+
  ggtitle('90 percent in main peak')+
  theme(plot.title = element_text(hjust = 0.5))+xlim(-.01, .12)+
  ylab('Density')+xlab('ITE')
multi5perc
#+annotate("text", x = .15, y = .089, label='mean inducement:',size=6)+annotate("text", x = .295, y = .089, label=100*round(multi5percmean$meangam,4),size=6)+annotate("text", x = .16, y = .139, label='mean bankruptcy prob:',size=6)+annotate("text", x = .295, y = .139, label=100*round(multi5percmean$meanB1,4),size=6)+xlim(-.1,.34)
setwd("/home/dpapakos/sensitivity_analysis/constraint_plots/")
ggsave(filename="symmetric5percbump_constrained.pdf",height=7, width=7,  plot=multi5perc)



print( proc.time() - ptm )
