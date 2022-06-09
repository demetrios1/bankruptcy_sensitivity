setwd("/home/dpapakos/moderating_variables/")

library(dplyr)
library(foreign)
library(tidyverse)
library(foreach)
library(bcf)
#library(fastbart)
library(dbarts)
library(foreach)
library(doParallel)

rm(list = ls())
# load data

data      = read.dta( "AuditorFinal20160419.dta" )
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
#intframe=BARTpred(data_sasub, treat=going_concern, Outcome=bankrptobs,vars=var0)
#write.csv(intframe, "/home/dpapakos/moderating_variables/mod_var_mbart.csv")
#read_csv()
library(readr)
data_sasub
#bartstuff=read.csv('RawScores4_monotone_bart_lgc_020160907.csv')
###all firms
#merging in the going concern=1, no going concern
#data_all<-bartstuff %>%
#  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
 #               RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade, bprob_gc0, bprob_gc1)%>%
#  group_by(year,going_concern)%>%
 # mutate(nogc=year<=year&going_concern==0,
#         gc=year<=year&going_concern==1)

#data_all_2<-as.data.frame(data_all[, c('bankrptobs',vargc)])
#bankruptcy=as.data.frame(data_all[, c('bankrptobs')])
#xcov=data_all[, c(vargc)]




#dd_2<-subset(data_all_2)
library(rpart)
library(readr)

getwd()
root_path <- "/home/dpapakos/sensitivity_analysis/constrained_integration_newvar/"
#outcome=read_csv("/home/dpapakos/sensitivity_analysis/constrained_integration/inducement_normal_mean0_sd0.5.csv")
outcome = read_csv(paste0(root_path,'inducement_normal_mean0_sd0.5.csv' ))
outcome$RR=outcome$B1/outcome$B0
#outcome=read.csv('Inducement_lgc0_monotone_multimodes0.1.csv') 
dataset=data.frame(data_sasub, gamma=outcome$RR, B1=outcome$B1, B0=outcome$B0, 
                   diff=outcome$B1-outcome$B0)
dataset=dataset[, c(vars, 'gamma', 'B1', 'B0', 'diff')]

colnames(dataset)=c('Logassets', 'Leverage', 'Investment', 
                    'Cash', 'ROA', 'Log Price', 'Intangible assets', 
                    'R&D', 'R&D missing',      'Non-audit fees', 'Non-audit fees missing', 
                    'No S&P Rating', 
                    'Rating below CCC+','Years client' ,'Rating downgrade', 
                    'Average short interest', 'Short interest ratio', 
                    'Sum of log returns',
                    'Return volatility', '2000', '2001', '2002', '2003', 
                    '2004', '2005', '2006', '2007',  '2008', '2009', 
                    '2010', '2011',  '2012', '2013', '2014',
                    'gamma', 
                    'B1', 'B0', 'diff')

rrr<-rpart( gamma~Logassets+Leverage+Investment+
            Cash+ROA+`Log Price`+`Intangible assets`+ 
            `R&D`+ `R&D missing`+ `No S&P Rating`+ 
            `Rating below CCC+`+ `Rating downgrade`+
            `Non-audit fees`+`Non-audit fees missing`+ 
            `Years client`+
              `Average short interest`+ `Short interest ratio`+
            `Sum of log returns`+
            `Return volatility`+ `2000`+ `2001`+ `2002`+ `2003`+
            `2004`+ `2005`+ `2006`+ `2007`+`2008`+ `2009`+ 
            `2010`+ `2011`+  `2012`+ `2013`+`2014`,data=dataset, minbucket=1000)
B1tree <- rpart( B1~Logassets+Leverage+Investment+
                   Cash+ROA+`Log Price`+`Intangible assets`+ 
                   `R&D`+ `R&D missing`+ `No S&P Rating`+ 
                   `Rating below CCC+`+ `Rating downgrade`+
                   `Non-audit fees`+`Non-audit fees missing`+ 
                   `Years client`+
                   `Average short interest`+ `Short interest ratio`+
                   `Sum of log returns`+
                   `Return volatility`+ `2000`+ `2001`+ `2002`+ `2003`+ 
                 `2004`+ `2005`+ `2006`+ `2007`+`2008`+ `2009`+ 
                   `2010`+ `2011`+  `2012`+ `2013`+`2014`,data=dataset, minbucket=1000)
B0tree <- rpart( B0~Logassets+Leverage+Investment+
                   Cash+ROA+`Log Price`+`Intangible assets`+ 
                   `R&D`+ `R&D missing`+ `No S&P Rating`+ 
                   `Rating below CCC+`+ `Rating downgrade`+
                   `Non-audit fees`+`Non-audit fees missing`+ 
                   `Years client`+
                   `Average short interest`+ `Short interest ratio`+
                   `Sum of log returns`+
                   `Return volatility`+ `2000`+ `2001`+ `2002`+ `2003`+
                 `2004`+ `2005`+ `2006`+ `2007`+`2008`+ `2009`+ 
                   `2010`+ `2011`+  `2012`+ `2013`+`2014`,data=dataset, minbucket=1000)
treattree <- rpart( diff~Logassets+Leverage+Investment+
                      Cash+ROA+`Log Price`+`Intangible assets`+ 
                      `R&D`+ `R&D missing`+ `No S&P Rating`+ 
                      `Rating below CCC+`+ `Rating downgrade`+
                      `Non-audit fees`+`Non-audit fees missing`+ 
                      `Years client`+
                      `Average short interest`+ `Short interest ratio`+
                      `Sum of log returns`+
                      `Return volatility`+ `2000`+ `2001`+ `2002`+ `2003`+
                    `2004`+ `2005`+ `2006`+ `2007`+`2008`+ `2009`+ 
                      `2010`+ `2011`+  `2012`+ `2013`+`2014`,data=dataset, minbucket=1000)


#+ey+dt+kpmg+pwc+gt+bdo
vars
summary(rrr)
#library(rattle)
library(partykit)
rrr$variable.importance
rrr$frame
#save as cart_tree_smallertree_newvar.pdf in the moderating variables/new_vars_trees folder
#save as 8 x 6
rpart.plot::rpart.plot(rrr, extra=1)
#save as cart_tree_B1_newvar.pdf in the moderating variables/new_vars_trees folder
rpart.plot::rpart.plot(B1tree, extra=1)
#save as cart_tree_B0_newvar.pdf in the moderating variables/new_vars_trees folder
rpart.plot::rpart.plot(B0tree, extra=1)
#save as cart_tree_treatment_newvar.pdf in the moderating variables/new_vars_trees folder
rpart.plot::rpart.plot(treattree, extra=1)

#### The RR case
bottom_nodes_RR=unique(rrr$where)
indices_of_interest_RR=sapply(1:length(bottom_nodes_RR), function(i)
  which(rrr$where==bottom_nodes_RR[i]))
index_RR_length <- sapply(1:length(indices_of_interest_RR), function(i) length(indices_of_interest_RR[[i]]))
index_RR_frame <- data.frame(n=rrr$frame$n,value=rrr$frame$yval)
min_RR <- index_RR_frame[which(index_RR_frame$value == min(index_RR_frame$value)), ]
max_RR <- index_RR_frame[which(index_RR_frame$value == max(index_RR_frame$value)), ]
smallest_index_RR<- indices_of_interest_RR[[which(index_RR_length == min_RR$n)]]
largest_index_RR <- indices_of_interest_RR[[which(index_RR_length == max_RR$n)]]
#### The treat
bottom_nodes_treat=unique(treattree$where)
indices_of_interest_treat=sapply(1:length(bottom_nodes_treat), function(i)
  which(treattree$where==bottom_nodes_treat[i]))
index_treat_length <- sapply(1:length(indices_of_interest_treat), function(i) length(indices_of_interest_treat[[i]]))
index_treat_frame <- data.frame(n=treattree$frame$n,value=treattree$frame$yval)
min_treat <- index_treat_frame[which(index_treat_frame$value == min(index_treat_frame$value)), ]
max_treat <- index_treat_frame[which(index_treat_frame$value == max(index_treat_frame$value)), ]
smallest_index_treat<- indices_of_interest_treat[[which(index_treat_length == min_treat$n)]]
largest_index_treat <- indices_of_interest_treat[[which(index_treat_length == max_treat$n)]]
####B1
bottom_nodes_B1=unique(B1tree$where)
indices_of_interest_B1=sapply(1:length(bottom_nodes_B1), function(i)
  which(B1tree$where==bottom_nodes_B1[i]))
index_B1_length <- sapply(1:length(indices_of_interest_B1), function(i) length(indices_of_interest_B1[[i]]))
index_B1_frame <- data.frame(n=B1tree$frame$n,value=B1tree$frame$yval)
min_B1 <- index_B1_frame[which(index_B1_frame$value == min(index_B1_frame$value)), ]
max_B1 <- index_B1_frame[which(index_B1_frame$value == max(index_B1_frame$value)), ]
smallest_index_B1<- indices_of_interest_B1[[which(index_B1_length == min_B1$n)]]
largest_index_B1 <- indices_of_interest_B1[[which(index_B1_length == max_B1$n)]]
####B0
bottom_nodes_B0=unique(B0tree$where)
indices_of_interest_B0=sapply(1:length(bottom_nodes_B0), function(i)
  which(B0tree$where==bottom_nodes_B0[i]))
index_B0_length <- sapply(1:length(indices_of_interest_B0), function(i) length(indices_of_interest_B0[[i]]))
index_B0_frame <- data.frame(n=B0tree$frame$n,value=B0tree$frame$yval)
min_B0 <- index_B0_frame[which(index_B0_frame$value == min(index_B0_frame$value)), ]
max_B0 <- index_B0_frame[which(index_B0_frame$value == max(index_B0_frame$value)), ]
smallest_index_B0<- indices_of_interest_B0[[which(index_B0_length == min_B0$n)]]
largest_index_B0 <- indices_of_interest_B0[[which(index_B0_length == max_B0$n)]]




root_path_post <- '/home/dpapakos/cred_intervals/credinterval_constrained_integral_newvar/'
##sdpt5
fill3_B1=data.table::fread(file = paste0(root_path_post, 'mean0sd0.5_credint_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill3_B0=data.table::fread(file = paste0(root_path_post, 'mean0sd0.5_credint_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]

##sd1
sd1_B1=data.table::fread(file = paste0(root_path_post, 'mean0sd1_credint_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
sd1_B0=data.table::fread(file = paste0(root_path_post, 'mean0sd1_credint_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]

#right bump
rightbump_B1=data.table::fread(file = paste0(root_path_post, 'rightbumpsigma0.48_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
rightbump_B0=data.table::fread(file = paste0(root_path_post, 'rightbumpsigma0.48_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]

#90%
ninetyperc_B1=data.table::fread(file = paste0(root_path_post, 'symmetric90percmid_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
ninetyperc_B0=data.table::fread(file = paste0(root_path_post, 'symmetric90percmid_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]

#the largest for the 0.5 case
sdpt5_large=t(fill3_B1)[largest_index_RR,]/t(fill3_B0)[largest_index_RR,]
dim(sdpt5_large)
#the smallest for the 0.5 case
sdpt5_small=t(fill3_B1)[smallest_index_RR,]/t(fill3_B0)[smallest_index_RR,]
#we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
#effect for each posterior trial across all firms in the subset
sdpt5_diff=colMeans(sdpt5_large)-colMeans(sdpt5_small)

#the largest for the 1 case
sd1_large=t(sd1_B1)[largest_index_RR,]/t(sd1_B0)[largest_index_RR,]

#the smallest for the 1 case
sd1_small=t(sd1_B1)[smallest_index_RR,]/t(sd1_B0)[smallest_index_RR,]
#we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
#effect for each posterior trial across all firms in the subset
sd1_diff=colMeans(sd1_large)-colMeans(sd1_small)

#the largest for the right bump case
rightbump_large=t(rightbump_B1)[largest_index_RR,]/t(rightbump_B0)[largest_index_RR,]

#the smallest for the right bump case
rightbump_small=t(rightbump_B1)[smallest_index_RR,]/t(rightbump_B0)[smallest_index_RR,]
#we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
#effect for each posterior trial across all firms in the subset
rightbump_diff=colMeans(rightbump_large)-colMeans(rightbump_small)

#the largest for the 90% case
ninetyperc_large=t(ninetyperc_B1)[largest_index_RR,]/t(ninetyperc_B0)[largest_index_RR,]

#the smallest for the 90% case
ninetyperc_small=t(ninetyperc_B1)[smallest_index_RR,]/t(ninetyperc_B0)[smallest_index_RR,]
#we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
#effect for each posterior trial across all firms in the subset
ninetyperc_diff=colMeans(ninetyperc_large)-colMeans(ninetyperc_small)



 mean(sdpt5_diff)
 mean(sd1_diff)
 mean(rightbump_diff)
 mean(ninetyperc_diff)

#  #name RR_post_diff_newvars.pdf in the new_vars_trees with 7x5 pdf
 plot(density(sdpt5_diff), col='dodgerblue4',main='', ylim=c(0, 0.1),xlim=c(0, 200),
      lwd=2,lty=1, xlab='Posterior Inducement')
 lines(density(sd1_diff),lty=2, col='firebrick4', lwd=2)
 lines(density(rightbump_diff),lty=3, col='darkorchid4', lwd=2)
 lines(density(ninetyperc_diff), lty=4,col='black', lwd=2)
 legend("topright", legend=c(expression(paste("f" [1],"(u)")), expression(paste("f" [2],"(u)")), 
                        expression(paste("f" [3],"(u)")), expression(paste("f" [4],"(u)"))),
       col=c("dodgerblue4", "firebrick4", 'darkorchid4', 'black'), lty=c(1,2,3,4), cex=1.1)



#### Repeat but for the treatment ####

 #the largest for the 0.5 case
 sdpt5_large=t(fill3_B1)[largest_index_treat,]-t(fill3_B0)[largest_index_treat,]
 dim(sdpt5_large)
 #the smallest for the 0.5 case
 sdpt5_small=t(fill3_B1)[smallest_index_treat,]-t(fill3_B0)[smallest_index_treat,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 sdpt5_diff=colMeans(sdpt5_large)-colMeans(sdpt5_small)
 
 #the largest for the 1 case
 sd1_large=t(sd1_B1)[largest_index_treat,]-t(sd1_B0)[largest_index_treat,]
 
 #the smallest for the 1 case
 sd1_small=t(sd1_B1)[smallest_index_treat,]-t(sd1_B0)[smallest_index_treat,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 sd1_diff=colMeans(sd1_large)-colMeans(sd1_small)
 
 #the largest for the right bump case
 rightbump_large=t(rightbump_B1)[largest_index_treat,]-t(rightbump_B0)[largest_index_treat,]
 
 #the smallest for the right bump case
 rightbump_small=t(rightbump_B1)[smallest_index_treat,]-t(rightbump_B0)[smallest_index_treat,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 rightbump_diff=colMeans(rightbump_large)-colMeans(rightbump_small)
 
 #the largest for the 90% case
 ninetyperc_large=t(ninetyperc_B1)[largest_index_treat,]-t(ninetyperc_B0)[largest_index_treat,]
 
 #the smallest for the 90% case
 ninetyperc_small=t(ninetyperc_B1)[smallest_index_treat,]-t(ninetyperc_B0)[smallest_index_treat,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 ninetyperc_diff=colMeans(ninetyperc_large)-colMeans(ninetyperc_small)
 
 
 
 
 mean(sdpt5_diff)
 mean(sd1_diff)
 mean(rightbump_diff)
 mean(ninetyperc_diff)
 #  #name treat_post_diff_newvars.pdf in the new_vars_trees with 7x5 pdf
 plot(density(sdpt5_diff), col='dodgerblue4',main='', ylim=c(0, 60),xlim=c(0, .2),
      lwd=2, lty=1,xlab='Posterior Treatment')
 lines(density(sd1_diff),lty=2, col='firebrick4', lwd=2)
 lines(density(rightbump_diff), lty=3,col='darkorchid4', lwd=2)
 lines(density(ninetyperc_diff), lty=4,col='black', lwd=2)
 legend("topright", legend=c(expression(paste("f" [1],"(u)")), expression(paste("f" [2],"(u)")), 
                             expression(paste("f" [3],"(u)")), expression(paste("f" [4],"(u)"))),
        col=c("dodgerblue4", "firebrick4",
              'darkorchid4', 'black'), lty=c(1,2,3,4), cex=1.1)
  
####B1 tree####
 #### Repeat but for the B1 tree ####

 #the largest for the 0.5 case
 sdpt5_large=t(fill3_B1)[largest_index_B1,]#/t(fill3_B0)[largest_index_B1,]
 dim(sdpt5_large)
 #the smallest for the 0.5 case
 sdpt5_small=t(fill3_B1)[smallest_index_B1,]#/t(fill3_B0)[smallest_index_B1,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 sdpt5_diff=colMeans(sdpt5_large)-colMeans(sdpt5_small)
 
 #the largest for the 1 case
 sd1_large=t(sd1_B1)[largest_index_B1,]#/t(sd1_B0)[largest_index_B1,]
 
 #the smallest for the 1 case
 sd1_small=t(sd1_B1)[smallest_index_B1,]#/t(sd1_B0)[smallest_index_B1,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 sd1_diff=colMeans(sd1_large)-colMeans(sd1_small)
 
 #the largest for the right bump case
 rightbump_large=t(rightbump_B1)[largest_index_B1,]#/t(rightbump_B0)[largest_index_B1,]
 
 #the smallest for the right bump case
 rightbump_small=t(rightbump_B1)[smallest_index_B1,]#/t(rightbump_B0)[smallest_index_B1,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 rightbump_diff=colMeans(rightbump_large)-colMeans(rightbump_small)
 
 #the largest for the 90% case
 ninetyperc_large=t(ninetyperc_B1)[largest_index_B1,]#/t(ninetyperc_B0)[largest_index_B1,]
 
 #the smallest for the 90% case
 ninetyperc_small=t(ninetyperc_B1)[smallest_index_B1,]#/t(ninetyperc_B0)[smallest_index_B1,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 ninetyperc_diff=colMeans(ninetyperc_large)-colMeans(ninetyperc_small)
 
 
 
 
 mean(sdpt5_diff)
 mean(sd1_diff)
 mean(rightbump_diff)
 mean(ninetyperc_diff)
 # 
 library(latex2exp)
 #name B1_post_diff_newvars.pdf in the new_vars_trees with 7x5 pdf
 plot(density(sdpt5_diff), col='dodgerblue4',main='', ylim=c(0, 60),xlim=c(0, .2),
      lwd=2, lty=1,xlab=TeX(sprintf("$Posterior\\; \\Pr(B=1 | x, do(G=1))$", alpha)))
 lines(density(sd1_diff), lty=2,col='firebrick4', lwd=2)
 lines(density(rightbump_diff), lty=3,col='darkorchid4', lwd=2)
 lines(density(ninetyperc_diff), lty=4,col='black', lwd=2)
 legend("topright", legend=c(expression(paste("f" [1],"(u)")), expression(paste("f" [2],"(u)")), 
                             expression(paste("f" [3],"(u)")), expression(paste("f" [4],"(u)"))),
        col=c("dodgerblue4", "firebrick4", 
              'darkorchid4', 'black'), lty=c(1,2,3,4), cex=1.1)
 
 
 ####B0 tree####
 #### Repeat but for the B0 tree ####
 #the largest for the 0.5 case
 sdpt5_large=t(fill3_B0)[largest_index_B0,]
 dim(sdpt5_large)
 #the smallest for the 0.5 case
 sdpt5_small=t(fill3_B0)[smallest_index_B0,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 sdpt5_diff=colMeans(sdpt5_large)-colMeans(sdpt5_small)

 #the largest for the 1 case
 sd1_large=t(sd1_B0)[largest_index_B0,]
 
 #the smallest for the 1 case
 sd1_small=t(sd1_B0)[smallest_index_B0,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 sd1_diff=colMeans(sd1_large)-colMeans(sd1_small)
 
 #the largest for the right bump case
 rightbump_large=t(rightbump_B0)[largest_index_B0,]
 
 #the smallest for the right bump case
 rightbump_small=t(rightbump_B0)[smallest_index_B0,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 rightbump_diff=colMeans(rightbump_large)-colMeans(rightbump_small)
 
 #the largest for the 90% case
 ninetyperc_large=t(ninetyperc_B0)[largest_index_B0,]
 
 #the smallest for the 90% case
 ninetyperc_small=t(ninetyperc_B0)[smallest_index_B0,]
 #we have n_obs x n_post, so we want the column means which indicate the posterior inducement 
 #effect for each posterior trial across all firms in the subset
 ninetyperc_diff=colMeans(ninetyperc_large)-colMeans(ninetyperc_small)
 
 
 
 mean(sdpt5_diff)
 mean(sd1_diff)
 mean(rightbump_diff)
 mean(ninetyperc_diff)
 # 
 library(latex2exp)
 #name B0_post_diff_newvars.pdf in the new_vars_trees with 7x5 pdf
 plot(density(sdpt5_diff), col='dodgerblue4',main='', ylim=c(0, 60),xlim=c(0, .250),
      lwd=2,lty=1,xlab=TeX(sprintf("$Posterior\\; \\Pr(B=1 | x, do(G=0))$", alpha)))
 lines(density(sd1_diff), lty=2,col='firebrick4', lwd=2)
 lines(density(rightbump_diff),lty=3, col='darkorchid4', lwd=2)
 lines(density(ninetyperc_diff), lty=4,col='black', lwd=2)
 legend("topright", legend=c(expression(paste("f" [1],"(u)")), expression(paste("f" [2],"(u)")), 
                             expression(paste("f" [3],"(u)")), expression(paste("f" [4],"(u)"))),
        col=c("dodgerblue4", "firebrick4", 
              'darkorchid4', 'black'), lty=c(1,2,3,4), cex=1.1)
 
 
 
 
 
 
 
 #### Recreate TABLE 2 FOR THE LARGEST RD SUBGROUP ####
 library(haven)
 data      = as.data.frame(read_dta( "/home/dpapakos/moderating_variables/auditor_with_stratio.dta"))
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
 data_sasub$bankrptobs
 #count number of going concerns for largest RD
 num_gc_treat <- sum(data_sasub[largest_index_treat, 'going_concern'])
 #count number of bankruptcies for largest RD
 num_banks_treat <- sum(data_sasub[largest_index_treat, 'bankrptobs'])
 #count number of going concerns for largest RR
 num_gc_RR <- sum(data_sasub[largest_index_RR, 'going_concern'])
 #count number of bankruptcies for largest RR
 num_banks_RR <- sum(data_sasub[largest_index_RR, 'bankrptobs'])
 
 #get the observed values of probabilities (fitted with bart at least)
 intframe_temp=data.table::fread(file = "/home/dpapakos/cred_intervals/mod_var_mbart_newvar.csv", na.strings = c("", "NA", "#N/A"))
 divs = floor(dim(intframe_temp)[2]/3)
 intframe <- c()
 intframe[[1]] <- intframe_temp[,2:(divs+1)]
 intframe[[2]] <- intframe_temp[,(divs+2):(2*divs +1)]
 intframe[[3]] <- intframe_temp[,(2*divs+2):(3*divs +1)]
 #First 25350 are B1G1
 intframe_B1G1=intframe[[1]]
 #second 25350 are B1G0
 intframe_B1G0=intframe[[2]]
 intframe_G=intframe[[3]]
 rm(intframe)
 prop_score_mean=colMeans(intframe_G)
 #risk ratios
 RR=intframe_B1G1/intframe_B1G0
 #read N(0,.5)
 root_path_post <- '/home/dpapakos/cred_intervals/credinterval_constrained_integral_newvar/'
 fill3_B1=data.table::fread(file = paste0(root_path_post, 'mean0sd0.5_credint_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
 fill3_B0=data.table::fread(file = paste0(root_path_post, 'mean0sd0.5_credint_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]
 B1_frame <- t(fill3_B1)
 B0_frame <- t(fill3_B0)
 #calculate the biggest risk difference subgroup values
 RD_B1 <- mean(rowMeans(B1_frame[largest_index_treat,]))
 RD_B0 <- mean(rowMeans(B0_frame[largest_index_treat,]))
 RD_tau <- mean(rowMeans(B1_frame[largest_index_treat,])/rowMeans(B0_frame[largest_index_treat,]))
 RD_quantile <- quantile(rowMeans(B1_frame[largest_index_treat,])/rowMeans(B0_frame[largest_index_treat,]),
                         c(0.025, .5, .975))

 RD_obs <- t(RR)[largest_index_treat,]
 
 RD_obs_mean <- mean(rowMeans(RD_obs))

 xtable::xtable(data.frame(RD_obs_mean, RD_B0, RD_B1,
                           RD_tau, RD_quantile[1], 
                           RD_quantile[3]))
 #calculate the biggest risk ratio subgroup values
 RR_B1 <- mean(rowMeans(B1_frame[largest_index_RR,]))
 RR_B0 <- mean(rowMeans(B0_frame[largest_index_RR,]))
 RR_tau <- mean(rowMeans(B1_frame[largest_index_RR,])/rowMeans(B0_frame[largest_index_RR,]))
 RR_quantile <- quantile(rowMeans(B1_frame[largest_index_RR,])/rowMeans(B0_frame[largest_index_RR,]),
                         c(0.025, .5, .975))
 
 RR_obs <- t(RR)[largest_index_RR,]
 
 RR_obs_mean <- mean(rowMeans(RR_obs))
 xtable::xtable(data.frame(RR_obs_mean, RR_B0, RR_B1,
                           RR_tau, RR_quantile[1], 
                           RR_quantile[3]))
 #### do the same for asymmetric mixture ####
 asym_B1=data.table::fread(file = paste0(root_path_post, 'rightbumpsigma0.48_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
 asym_B0=data.table::fread(file = paste0(root_path_post, 'rightbumpsigma0.48_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]
 B0_frame_asym=t(data.frame(asym_B0))
 B1_frame_asym=t(data.frame(asym_B1))

 
 RD_B1_asym <- mean(rowMeans(B1_frame_asym[largest_index_treat,]))
 RD_B0_asym <- mean(rowMeans(B0_frame_asym[largest_index_treat,]))
 RD_tau_asym <- mean(rowMeans(B1_frame_asym[largest_index_treat,])/rowMeans(B0_frame_asym[largest_index_treat,]))
 RD_quantile_asym  <- quantile(rowMeans(B1_frame_asym[largest_index_treat,])/rowMeans(B0_frame_asym[largest_index_treat,]),
         c(0.025, .5, .975))
 xtable::xtable(data.frame(RD_B0_asym, RD_B1_asym,
                           RD_tau_asym, RD_quantile_asym[1], 
                           RD_quantile_asym[3]))
 #same for the risk ratio largest
 RR_B1_asym <- mean(rowMeans(B1_frame_asym[largest_index_RR,]))
 RR_B0_asym <- mean(rowMeans(B0_frame_asym[largest_index_RR,]))
 RR_tau_asym <- mean(rowMeans(B1_frame_asym[largest_index_RR,])/rowMeans(B0_frame_asym[largest_index_RR,]))
 RR_quantile_asym  <- quantile(rowMeans(B1_frame_asym[largest_index_RR,])/rowMeans(B0_frame_asym[largest_index_RR,]),
                               c(0.025, .5, .975))
 
xtable::xtable(data.frame(RR_B0_asym, RR_B1_asym,
                          RR_tau_asym, RR_quantile_asym[1], 
                          RR_quantile_asym[3]))

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

 
 
 ####hacky but have to do, the ordering is weird!
 # #want it to line up with tree
 # subgroup_RR<-c()
 # subgroup_RR[[1]]=mean(dataset[indices_of_interest[[1]], 'gamma'])
 # subgroup_RR[[2]]=mean(dataset[indices_of_interest[[3]], 'gamma'])
 # subgroup_RR[[3]]=mean(dataset[indices_of_interest[[2]], 'gamma'])
 # subgroup_RR[[4]]=mean(dataset[indices_of_interest[[5]], 'gamma'])
 # subgroup_RR[[5]]=mean(dataset[indices_of_interest[[4]], 'gamma'])
 # 
 # subgroup_B1<-c()
 # subgroup_B1[[1]]=mean(dataset[indices_of_interest[[1]], 'B1'])
 # subgroup_B1[[2]]=mean(dataset[indices_of_interest[[3]], 'B1'])
 # subgroup_B1[[3]]=mean(dataset[indices_of_interest[[2]], 'B1'])
 # subgroup_B1[[4]]=mean(dataset[indices_of_interest[[5]], 'B1'])
 # subgroup_B1[[5]]=mean(dataset[indices_of_interest[[4]], 'B1'])
 # 
 # subgroup_B0<-c()
 # subgroup_B0[[1]]=mean(dataset[indices_of_interest[[1]], 'B0'])
 # subgroup_B0[[2]]=mean(dataset[indices_of_interest[[3]], 'B0'])
 # subgroup_B0[[3]]=mean(dataset[indices_of_interest[[2]], 'B0'])
 # subgroup_B0[[4]]=mean(dataset[indices_of_interest[[5]], 'B0'])
 # subgroup_B0[[5]]=mean(dataset[indices_of_interest[[4]], 'B0'])
 # round(subgroup_B1,3)
 # #subgroup_RR=sapply(1:length(indices_of_interest), function(i)mean(dataset[indices_of_interest[[i]], 'gamma']))
 # data.frame(subgroup_RR, subgroup_B1, subgroup_B0)
 # xtable::xtable(data.frame(subgroup_RR, B1=round(subgroup_B1,3), B0=round(subgroup_B0,3)))
 # #write.csv(data.frame(subgroup_RR, subgroup_B1, subgroup_B0), '/home/dpapakos/moderating_variables/treat_tree_RR_newvar.csv')
 # 
 # 
 # 
 # plot_treat_frame=data.frame(subgroup=c(rep('1',length(mean(dataset[indices_of_interest[[1]], 
 #                                                                 'gamma']))), 
 #   rep('2',length(mean(dataset[indices_of_interest[[3]], 'gamma']))), 
 #   rep('3',length(mean(dataset[indices_of_interest[[2]], 'gamma']))), 
 #   rep('4',length(mean(dataset[indices_of_interest[[5]], 'gamma']))),
 #   rep('5',length(mean(dataset[indices_of_interest[[4]], 'gamma'])))), 
 #   inducement=c(mean(dataset[indices_of_interest[[1]], 'gamma']), 
 #                mean(dataset[indices_of_interest[[3]], 'gamma']), 
 # mean(dataset[indices_of_interest[[2]], 'gamma']), 
 #  mean(dataset[indices_of_interest[[5]], 'gamma']),
 #  mean(dataset[indices_of_interest[[4]], 'gamma'])), 
 # treatment=c(mean(dataset[indices_of_interest[[1]], 'diff']), 
 #   mean(dataset[indices_of_interest[[3]], 'diff']), 
 #   mean(dataset[indices_of_interest[[2]], 'diff']), 
 #   mean(dataset[indices_of_interest[[5]], 'diff']),
 #   mean(dataset[indices_of_interest[[4]], 'diff'])))
 # plot_treat_frame
 # library(viridis)
 # plot_treat_frame%>%
 #   ggplot(aes(x=treatment, y=inducement, color=subgroup))+geom_point()+
 #   theme_minimal()+scale_color_viridis(discrete = TRUE, option = "D")+
 #   ggtitle('Inducement vs Treatment from treatment tree')+
 #   theme(plot.title = element_text(hjust = 0.5,size=14))
 # 
 # B1tree
 #' 
 #' ####now repeat for the B1 tree####
 #' bottom_nodes_B1=unique(B1tree$where)
 #' indices_B1 = sapply(1:length(bottom_nodes_B1), function(i)
 #'   which(B1tree$where==bottom_nodes_B1[i]))
 #' indices_B1
 #' sapply(1:length(indices_B1), function(i)
 #'   length(indices_B1[[i]]))
 #' ####the correct indexing####
 #' correct_B1 <- c(1, 5, 4, 2, 6, 3, 7)
 #' #just doing this manually...oh well!
 #' correct_B1 <- c(1,2,3,4)
 #' B1_max_index = which(correct_B1==max(correct_B1))
 #' ####now repeat for the B0 tree####
 #' bottom_nodes_B0=unique(B0tree$where)
 #' indices_B0 = sapply(1:length(bottom_nodes_B0), function(i)
 #'   which(B0tree$where==bottom_nodes_B0[i]))
 #' sapply(1:length(indices_B0), function(i)
 #'   length(indices_B0[[i]]))
 #' ####the correct indexing####
 #' correct_B0 <- c(1,3,2,5,4)
 #' corect_B0 <- c(1,2,3,4)
 #' B0_max_index = which(correct_B0==max(correct_B0))
 #' 
 #' 
 #' #' 
 #' #' ####do the same for the RR tree with the treats ####
 #' #' ##this is confusing, but we want the treatment for the risk ratio tree, should
 #' #' #have named it the other direction tbh
 #' #' bottom_nodes_RR=unique(rrr$where)
 #' indices_of_interest_treat=sapply(1:length(bottom_nodes_RR), function(i)
 #'    which(rrr$where==bottom_nodes_RR[i]))
 
 #' ##the ordering is wrong
 #' #1->1
 #' #2->3
 #' #3->4
 #' #4->2
 #' #5->5
 #' #6->6
 #' 
 #' 
 #' #### actually, now we have 4 ->1
 #' #2 -> 2
 #' # 3->1
 #' # 4->5
 #' #5 -> 3
 #' 
 #' subgroup_treat<-c()
 #' subgroup_treat[[1]]=mean(dataset[indices_of_interest_treat[[1]], 'diff'])
 #' subgroup_treat[[2]]=mean(dataset[indices_of_interest_treat[[3]], 'diff'])
 #' subgroup_treat[[3]]=mean(dataset[indices_of_interest_treat[[4]], 'diff'])
 #' subgroup_treat[[4]]=mean(dataset[indices_of_interest_treat[[2]], 'diff'])
 #' subgroup_treat[[5]]=mean(dataset[indices_of_interest_treat[[5]], 'diff'])
 #' subgroup_treat[[6]]=mean(dataset[indices_of_interest_treat[[6]], 'diff'])
 #' subgroup_B1_treat<-c()
 #' subgroup_B1_treat[[1]]=mean(dataset[indices_of_interest_treat[[1]], 'B1'])
 #' subgroup_B1_treat[[2]]=mean(dataset[indices_of_interest_treat[[3]], 'B1'])
 #' subgroup_B1_treat[[3]]=mean(dataset[indices_of_interest_treat[[4]], 'B1'])
 #' subgroup_B1_treat[[4]]=mean(dataset[indices_of_interest_treat[[2]], 'B1'])
 #' subgroup_B1_treat[[5]]=mean(dataset[indices_of_interest_treat[[5]], 'B1'])
 #' subgroup_B1_treat[[6]]=mean(dataset[indices_of_interest_treat[[6]], 'B1'])
 #' subgroup_B0_treat<-c()
 #' subgroup_B0_treat[[1]]=mean(dataset[indices_of_interest_treat[[1]], 'B0'])
 #' subgroup_B0_treat[[2]]=mean(dataset[indices_of_interest_treat[[3]], 'B0'])
 #' subgroup_B0_treat[[3]]=mean(dataset[indices_of_interest_treat[[4]], 'B0'])
 #' subgroup_B0_treat[[4]]=mean(dataset[indices_of_interest_treat[[2]], 'B0'])
 #' subgroup_B0_treat[[5]]=mean(dataset[indices_of_interest_treat[[5]], 'B0'])
 #' subgroup_B0_treat[[6]]=mean(dataset[indices_of_interest_treat[[6]], 'B0'])
 #' round(subgroup_B1_treat,3)
 #' #subgroup_RR=sapply(1:length(indices_of_interest), function(i)mean(dataset[indices_of_interest[[i]], 'gamma']))
 #' data.frame(subgroup_treat, subgroup_B1_treat, subgroup_B0_treat)
 #' xtable::xtable(data.frame(subgroup_treat, B1=round(subgroup_B1_treat,3), 
 #'                           B0=round(subgroup_B0_treat,3)))
 #' #write.csv(data.frame(subgroup_treat, subgroup_B1_treat, subgroup_B0_treat), 
 #'           #'/home/dpapakos/moderating_variables/RR_tree_treat.csv')
 #' plot_RR_frame=data.frame(subgroup=c(rep('1',1),#length(dataset[indices_of_interest_treat[[1]], 'gamma'])), 
 #'                                     rep('2',1),#length(dataset[indices_of_interest_treat[[3]], 'gamma'])), 
 #'                                     rep('3',1),#length(dataset[indices_of_interest_treat[[4]], 'gamma'])), 
 #'                                     rep('4',1),#length(dataset[indices_of_interest_treat[[2]], 'gamma'])),
 #'                                     rep('5',1),#length(dataset[indices_of_interest_treat[[5]], 'gamma'])), 
 #'                                     rep('6', 1)),#length(dataset[indices_of_interest_treat[[6]], 'gamma']))), 
 #'                          inducement=c(mean(dataset[indices_of_interest_treat[[1]], 'gamma']), 
 #'                                       mean(dataset[indices_of_interest_treat[[3]], 'gamma']), 
 #'                                       mean(dataset[indices_of_interest_treat[[4]], 'gamma']), 
 #'                                       mean(dataset[indices_of_interest_treat[[2]], 'gamma']),
 #'                                       mean(dataset[indices_of_interest_treat[[5]], 'gamma']), 
 #'                                       mean(dataset[indices_of_interest_treat[[6]], 'gamma'])), 
 #'                          treatment=c(mean(dataset[indices_of_interest_treat[[1]], 'diff']), 
 #'                                      mean(dataset[indices_of_interest_treat[[3]], 'diff']), 
 #'                                      mean(dataset[indices_of_interest_treat[[4]], 'diff']), 
 #'                                      mean(dataset[indices_of_interest_treat[[2]], 'diff']),
 #'                                      mean(dataset[indices_of_interest_treat[[5]], 'diff']),
 #'                                      mean(dataset[indices_of_interest_treat[[6]], 'diff'])
 #'                                      ))
 #' plot_RR_frame
 #' #library(viridis)
 #' plot_RR_frame%>%
 #'   ggplot(aes(x=treatment, y=inducement, color=subgroup))+geom_point()+
 #'   theme_minimal()+scale_color_viridis(discrete = TRUE, option = "D")+
 #'   ggtitle('Inducement vs Treatment from inducement tree')+
 #'   theme(plot.title = element_text(hjust = 0.5,size=14))
 #' 
 
 
 ##double check the size issue!
 #sapply(1:5, function(i)length(indices_of_interest[[i]]))
 # party=as.party(rrr)
 # party
 # 
 # rrr$frame
 # #?rpart.plot
 # bottom_row_index=unique(rrr$where)
 # bottom_row_index
 # sapply(1:length(bottom_nodes_RR), function(i)length(indices_of_interest_treat[[i]]))
 # correct_RR <- c(4,2,3,5, 3)
 # ####these are a bit off, the largest is actually second largest because
 # #it splits to left not right
 # largest_index=bottom_row_index[length(bottom_row_index)]
 # second_largest_index=bottom_row_index[length(bottom_row_index)-1]
 # smallest_index=bottom_row_index[1]
 # second_smallest_index=bottom_row_index[2]
 # 
 # ##sanity check!##
 # length(which(rrr$where==smallest_index))
 # length(which(rrr$where==largest_index))
 
 
 #write.csv(which(rrr$where==28), 'mu0sd0.5largestsub.csv')
 #rrr
 #the indices we need
 #do the two largest, and two others
 # 
 # indices_large=which(rrr$where==largest_index)#read_csv('mu0sd0.5largestsub.csv')[,2]
 # indices_small=which(rrr$where==smallest_index)
