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
#write.csv(intframe, "/home/dpapakos/moderating_variables/mod_var_mbart.csv")
#read_csv()
library(readr)
data_sasub
bartstuff=read.csv('RawScores4_monotone_bart_lgc_020160907.csv')
###all firms
#merging in the going concern=1, no going concern
data_all<-bartstuff %>%
  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
                RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade, bprob_gc0, bprob_gc1)%>%
  group_by(year,going_concern)%>%
  mutate(nogc=year<=year&going_concern==0,
         gc=year<=year&going_concern==1)

data_all_2<-as.data.frame(data_all[, c('bankrptobs',vargc)])
bankruptcy=as.data.frame(data_all[, c('bankrptobs')])
xcov=data_all[, c(vargc)]




dd_2<-subset(data_all_2)
library(rpart)
library(readr)

getwd()
outcome=read_csv("/home/dpapakos/sensitivity_analysis/constrained_integration/inducement_normal_mean0_sd0.5.csv")
outcome$RR=outcome$B1/outcome$B0
#outcome=read.csv('Inducement_lgc0_monotone_multimodes0.1.csv') 
dataset=data.frame(data_sasub, gamma=outcome$RR)
rrr<-rpart( gamma~logassets_r+ lev_r+ investments_r +cash_r+roa_r+logprice_r+Intangible_r+RD_r+
             RDmissing+fracnonfees_r+feemissing+NoRate+RateC+numyears+downgrade+
             ey+dt+kpmg+pwc+gt+bdo,data=dataset, minbucket=1000)
summary(rrr)
library(rattle)
library(partykit)
rrr$variable.importance
rrr$frame
rpart.plot::rpart.plot(rrr, extra=1)
party=as.party(rrr)
party
rrr$frame
?rpart.plot
bottom_row_index=unique(rrr$where)
bottom_row_index
####these are a bit off, the largest is actually second largest because
#it splits to left not right
largest_index=bottom_row_index[length(bottom_row_index)]
second_largest_index=bottom_row_index[length(bottom_row_index)-1]
smallest_index=bottom_row_index[1]
second_smallest_index=bottom_row_index[2]

##sanity check!##
length(which(rrr$where==smallest_index))
length(which(rrr$where==largest_index))


#write.csv(which(rrr$where==28), 'mu0sd0.5largestsub.csv')
#rrr
#the indices we need
#do the two largest, and two others

indices_large=which(rrr$where==largest_index)#read_csv('mu0sd0.5largestsub.csv')[,2]
indices_small=which(rrr$where==smallest_index)
##sdpt5
fill3_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd0.5_credint_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
fill3_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd0.5_credint_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

##sd1
sd1_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd1_credint_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
sd1_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd1_credint_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

#right bump
rightbump_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/rightbumpsigma0.05_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
rightbump_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/rightbumpsigma0.05_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

#90%
ninetyperc_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/symmetric90percmid_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
ninetyperc_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/symmetric90percmid_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

#the largest for the 0.5 case
sdpt5_large=t(fill3_B1)[indices_large,]/t(fill3_B0)[indices_large,]
dim(sdpt5_large)
#the smallest for the 0.5 case
sdpt5_small=t(fill3_B1)[indices_small,]/t(fill3_B0)[indices_small,]
#we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
#effect for each posterior trial across all firms in the subset
sdpt5_diff=colMeans(sdpt5_large)-colMeans(sdpt5_small)

#the largest for the 1 case
sd1_large=t(sd1_B1)[indices_large,]/t(sd1_B0)[indices_large,]

#the smallest for the 1 case
sd1_small=t(sd1_B1)[indices_small,]/t(sd1_B0)[indices_small,]
#we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
#effect for each posterior trial across all firms in the subset
sd1_diff=colMeans(sd1_large)-colMeans(sd1_small)

#the largest for the right bump case
rightbump_large=t(rightbump_B1)[indices_large,]/t(rightbump_B0)[indices_large,]

#the smallest for the right bump case
rightbump_small=t(rightbump_B1)[indices_small,]/t(rightbump_B0)[indices_small,]
#we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
#effect for each posterior trial across all firms in the subset
rightbump_diff=colMeans(rightbump_large)-colMeans(rightbump_small)

#the largest for the 90% case
ninetyperc_large=t(ninetyperc_B1)[indices_large,]/t(ninetyperc_B0)[indices_large,]

#the smallest for the 90% case
ninetyperc_small=t(ninetyperc_B1)[indices_small,]/t(ninetyperc_B0)[indices_small,]
#we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
#effect for each posterior trial across all firms in the subset
ninetyperc_diff=colMeans(ninetyperc_large)-colMeans(ninetyperc_small)



 mean(sdpt5_diff)
 mean(sd1_diff)
 mean(rightbump_diff)
 mean(ninetyperc_diff)
# 
 plot(density(sdpt5_diff), col='dodgerblue4',main='', ylim=c(0, .14),xlim=c(0, 120),
      lwd=2, xlab='Posterior Inducement')
 lines(density(sd1_diff), col='firebrick4', lwd=2)
 lines(density(rightbump_diff), col='darkorchid4', lwd=2)
 lines(density(ninetyperc_diff), col='black', lwd=2)
 legend("topright", legend=c(expression(paste("f" [1],"(u)")), expression(paste("f" [2],"(u)")), 
                        expression(paste("f" [3],"(u)")), expression(paste("f" [4],"(u)"))),
       col=c("dodgerblue4", "firebrick4", 'darkorchid4', 'black'), lty=1, cex=1.1)
# 
# fill4=matrix(NA,nrow=2000, ncol=length(indices))#dim(indices)[1])
# library(foreach)
# library(doParallel)
# ####helper functions####
# 
# np = detectCores()-1
# cl = makeCluster(np)
# n_cores <- detectCores() - 1
# registerDoParallel(cores = n_cores)
# for (i in 1:dim(intframe[[1]])[1]){
# 
# index_unl=unlist(indices)[]
# 
#   merge_frame=data.frame(treatpred=t(intframe[[1]][i,index_unl]), notreatpred=t(intframe[[2]][i,index_unl]),
#                          propensity=t(intframe[[3]][i,index_unl]))
# dim(merge_frame)
# #f=function(u){
# #  dnorm(u, mean=0, sd=1)
# #}
# #f=function(u){
# #  .01*dnorm(u, mean=-2, sd=0.05)+0.94*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
# #}
# f=function(u){
#   .05*dnorm(u, mean=-2, sd=0.05)+0.9*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
# }
# colnames(merge_frame)=c('treatpred', 'notreatpred','propensity')
#  # expoutcomesfun	   =cbind(data.frame(  treatpred = pred1, notreatpred=pred2), propensity=pred3 )
# #  expoutcomesfun
#  # expoutcomesfun2    = rbind( NULL, expoutcomesfun )
#   outcomenew=merge_frame[,c('propensity','notreatpred','treatpred')]
# 
#   ####need the joints####
#   
#   eps     = 1e-6
#   
#   outcomenew[,c("prop", "Y_D0", "Y_D1")] = 0.5*eps + (1-eps)*outcomenew[,c('propensity','notreatpred','treatpred')]
#   
#   outcomenew[,"ProbY1D1"]=outcomenew[,"prop"]*outcomenew[,"treatpred"]
#   outcomenew[,"ProbY1D0"] = (1-outcomenew[,"prop"])*outcomenew[,"notreatpred"]
#   outcomenew[,"ProbY0D1"] = outcomenew[,"prop"]*(1-outcomenew[,"treatpred"])
#  
#   #for error analysis, not super necessary
#   indexes=seq(from=1, to=length(outcomenew$treatpred), by=1)
#   outcomenew=as.data.frame(cbind(indexes, outcomenew))
#   outcomenew$indexes=as.numeric(outcomenew$indexes)
#   #row.names(outcomenew)<-NULL
#   newoutcome =as.matrix(outcomenew)
#  
#   dim(newoutcome)
#   ptm = proc.time()
# 
#   mBG_out = foreach( m = 1:dim(newoutcome)[1],  .combine=rbind )%dopar%{
#     # Probability to fit.
#     vProbBG     = newoutcome[m,c( "ProbY1D1", "ProbY1D0", "ProbY0D1" )] 
#     
#     
#     # Starting value
#     
#     start0 = c( newoutcome[m,"ProbY1D1"], newoutcome[m,"ProbY1D0"],newoutcome[m,"ProbY0D1"] )
#     
#     names(start0) = c("y1","y0","d")
#     
#       
#     #better global definition.
#     optimdist=function(vals){
#       y1 = vals[1]
#       y0 = vals[2]
#       d = vals[3]
#       a=integrate( function(u) pnorm(y1+u)*pnorm(d+u)*f(u), lower = -Inf, upper = Inf )$value
#       b=integrate(function(u) pnorm(y0+u)*(1-pnorm(d+u))*f(u), lower = -Inf, upper = Inf )$value
#       c=integrate(  function(u) (1-pnorm(y1+u))*pnorm(d+u)*f(u), lower = -Inf, upper = Inf )$value 
#       #return(c(a,b,c))
#       #vProbBG is lhs, global variable
#       return(sum((qnorm(c(a,b,c))-qnorm(vProbBG))^2))
#     }
#     optimdist_constraint=
#       function(vals){
#         y1 = vals[1]
#         y0 = vals[2]
#         d = vals[3]
#         lambda=0  #in case we want to regularize!
#         #replace y1 with y1+y0
#         a=integrate( function(u) pnorm((exp(y1)+y0)+u)*pnorm(d+u)*f(u), lower = -Inf, upper = Inf )$value
#         b=integrate(function(u) pnorm(y0+u)*(1-pnorm(d+u))*f(u), lower = -Inf, upper = Inf )$value
#         c=integrate(  function(u) (1-pnorm((exp(y1)+y0)+u))*pnorm(d+u)*f(u), lower = -Inf, upper = Inf )$value 
#         #return(c(a,b,c))
#         #vProbBG is lhs, global variable
#         return(sum((qnorm(c(a,b,c))-qnorm(vProbBG))^2)+lambda*((exp(y1)+y0) -y0)^2)
#       }
#     #calculate the treatment
#     treatment_Compute = function( vb){
#       integrate( function(u) (pnorm(vb[1]+u) - pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value
#       
#     }
#     treatment_Compute_constraint = function( vb){
#       integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u) - pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value
#       
#     }
#     counterfac_Compute = function(vb){
#       y1 =(integrate( function(u) (pnorm(vb[1]+u))*f(u), lower = -Inf, upper = Inf )$value)
#       y0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
#       
#       return(c(y1, y0))
#       
#     }
#     
#     counterfac_Compute_constraint = function(vb){
#       y1 =(integrate( function(u) (pnorm(exp(vb[1])+vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
#       y0 = (integrate( function(u) (pnorm(vb[2]+u))*f(u), lower = -Inf, upper = Inf )$value)
#       
#       return(c(y1, y0))
#       
#     }
#     
#     vBG_out     = try(optim( qnorm(start0), optimdist_constraint,  control = list(maxit = 500), method='Nelder-Mead' 
#     )
#     , silent=TRUE)
#     if ( class( vBG_out ) == "try-error" ){
#       out        = as.numeric(newoutcome[m,c("indexes")])
#       out        = c( out, rep(NA,5) ) 
#     } else {
#       out        = as.numeric(newoutcome[m,c("indexes")] )
#       vBG        = vBG_out$par 
#       tau      = treatment_Compute_constraint( vb = vBG[c(1,2)])
#       
#       counter=counterfac_Compute_constraint(vb = vBG[c(1,2)])
#       
#       out        = c( vBG, tau,vBG_out$convergence, vBG_out$value,counter) 
#     }
#     return( out)
#   }
# 
#   colnames(mBG_out) = c( "y1", "y0", "d", "tau", "convergence", "fnvalue" ,'B1','B0')
#   print( proc.time() - ptm )
# 
#  fill4[i, ]<-mBG_out[, 'B1']/mBG_out[, 'B0']
#   #saveframe=data.frame(mBG_out)
# 
#   
#   print(i)
# }
# setwd('/home/dpapakos/moderating_variables/RR_subs')
# #write.csv(fill, 'mean0sd0.5_largestATEsub.csv')
# write.csv(fill4, 'symm10percbump_largestCRRsub.csv')
# dim(fill)
# fill=read.csv('mean0sd0.5CRRsub.csv')
# fill2=read.csv('mean0sd1CRRsub.csv')
# fill3=read.csv('rightbump_largestCRRsub.csv')
# fill4=read.csv('symm10percbump_largestCRRsub.csv')
# 
# gamma=rowMeans(fill)
# dim(bartstuff)
# gamma2=rowMeans(fill2)
# gamma3=rowMeans(fill3)
# gamma4=rowMeans(fill4)
# mean(gamma)
# mean(gamma2)
# mean(gamma3)
# mean(gamma4)
# 
# plot(density(gamma), col='dodgerblue4',main='', ylim=c(0, .05),xlim=c(0, 200),
#      lwd=2, xlab='Posterior Inducement')
# lines(density(gamma2), col='firebrick4', lwd=2)
# lines(density(gamma3), col='darkorchid4', lwd=2)
# lines(density(gamma4), col='black', lwd=2)
# legend("topright", legend=c(expression(paste("f" [1],"(u)")), expression(paste("f" [2],"(u)")), 
#                        expression(paste("f" [3],"(u)")), expression(paste("f" [4],"(u)"))),
#        col=c("dodgerblue4", "firebrick4", 'darkorchid4', 'black'), lty=1, cex=1.1)

