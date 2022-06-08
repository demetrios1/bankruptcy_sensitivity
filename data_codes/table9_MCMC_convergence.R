
library(coda)

library(dplyr)
library(foreign)
library(tidyverse)
library(foreach)
library(bcf)
library(monbart)
library(dbarts)
library(foreach)
library(doParallel)

#rm(list = ls())




## simulate data (example from Friedman MARS paper)
## y = f(x) + epsilon , epsilon ~ N(0, sigma)
## x consists of 10 variables, only first 5 matter

f <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
    10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

## run BART
set.seed(99)
bartFit <- bart(x, y)
bartFit$n.chains
plot(bartFit)

library(here)
# load data
root_path <- here('Dropbox', 'audit_paper', 'data')
data      = read.dta( paste0(root_path, "AuditorFinal20160419.dta" ))
#get these data from the cik archive on the web
cik_ticker=read_delim(paste0(root_path, '/cik_ticker.csv'), 
                        "|", escape_double = FALSE, trim_ws = TRUE)
cik_ticker=cik_ticker%>%dplyr::select(Name, CIK)
data_CIK_fix=read_csv(paste0(root_path, '/data_CIK_fix.csv'))

#### updated data ####
library(haven)
data      = as.data.frame(read_dta(paste0(root_path,
                                          "auditor_with_stratio.dta")))
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

data_CIK_fix$gvkey <- data_sasub$gvkey
data_CIK_fix$datadate <- data_sasub$datadate
#dumb way to do this, but ensures matching of indices
data_sasub$CIK <- data_CIK_fix$CIK
data_sasub$cik <- data_CIK_fix$cik
data_sasub <- (data_sasub[complete.cases(data_sasub[, vars]), ])

#matching the size!

data_CIK_fix<-right_join(data_CIK_fix , data_sasub, by=c('gvkey', 'datadate'))%>%
  dplyr::select(gvkey, datadate, sig_date_of_op_s.x,CIK.x,cik.x,going_concern.x,bankrptobs.x)
colnames(data_CIK_fix) <- c('gvkey', 'datadate', 'sig_date_of_op_s','CIK','cik','going_concern','bankrptobs')
data_CIK <- right_join(cik_ticker,data_CIK_fix)
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
BARTpred=function(df, treat=going_concern, Outcome=bankrptobs,vars, model){
  covtreat=df%>%filter(going_concern==1)
  covcontrol=df%>%filter(going_concern==0)
  covtreat$treat=as.factor(covtreat$going_concern)
  
  #return(df[vars])
  #return(as.factor(covtreat$outcome))
  #case 1, the treatment
  if (model=='BART'){
  bart1=bart(covtreat[vars],as.factor(covtreat$bankrptobs), df[vars],ndpost = 2000, nskip = 2000,ntree=100,verbose=T,usequants = TRUE,numcut = 1000)
  
  
  #case 2 control
  bart2=bart(covcontrol[vars],as.factor(covcontrol$bankrptobs),df[vars],ndpost =2000, nskip = 2000,ntree=100,verbose=T,usequants = TRUE,numcut = 1000)
  pred1 = pnorm(bart1$yhat.test)#colMeans(pnorm(bart1$yhat.test))
  pred2 = pnorm(bart2$yhat.test)#colMeans(pnorm(bart1$yhat.test))
  #print('first time:', sigma_save_pr0)
  }else{
  #case 3 propensity
  #bart3=bart(df[,vars], as.factor(df$going_concern),df[,vars],ndpost = 2000, nskip = 2000,ntree=100,
  #           verbose=F,usequants = TRUE,numcut = 1000)
  
  
  xtest=df[vars]
  xtraincontrol=covcontrol[vars]
  xtraintreat=covtreat[vars]
  ytrain0    = as.factor( covcontrol$bankrptobs )
  ytrain1    = as.factor( covtreat$bankrptobs )
  
  # mono fits
  
  bart_mono = monbart::monotone_bart(y = as.numeric(c(ytrain1, ytrain0)==1),
                            z = 1-c(rep(1, length(ytrain1)), rep(0, length(ytrain0))),
                            x = rbind(xtraintreat, xtraincontrol),
                            xpred = xtest, nskip = 2000, ndpost = 2000,m=100)
  
  

  #bart3 = bart(df[vars], as.factor(df$G),df[vars], ndpost = 5000, nskip = 2000, ntree=100, usequants=T, numcuts=1000, verbose=F)
  
  set.seed(12296)
  post_list=sample.int(2000, 500, replace=F)
  
  #use this method for prediction on binary
  #pred1 = colMeans(pnorm(bart1$yhat.test)) 
  #pred2=colMeans(pnorm(bart2$yhat.test))
  pred1= 1-bart_mono$pr0#1-colMeans(bart_mono$pr0)
  pred2= 1-bart_mono$pr1#1-colMeans(bart_mono$pr1)
  }
  #pred1= 1-colMeans(bart_mono$pr0[])#dont index at post_list
  #pred2= 1-colMeans(bart_mono$pr1[]) #dont index at post_list
  #pred3=colMeans(pnorm(bart3$yhat.test[])) #dont index at post_list
  
  #pred1 = colMeans(pnorm(bart1$yhat.test)) 
  #pred2=colMeans(pnorm(bart2$yhat.test))
  #print('now we have' , sigma_save_pr0)
  return(list(data.frame(pred1), data.frame(pred2)))
  #expoutcomesfun	   = data.frame(  treatpred = pred1, notreatpred=pred2)#, propensity=pred3 )
  
  
  #return(expoutcomesfun)
}



BART_model = BARTpred(data_sasub, treat=going_concern, Outcome=bankrptobs,vars=var0, 
                      model="BART")
mono_model = BARTpred(data_sasub, treat=going_concern, Outcome=bankrptobs,vars=var0, 
                      model="monoBART")
#both of these are the observations across the columns and the rows are the posterior draws

dim(mono_model)


mcmc(mono_model[[2]][,1])

acf(mcmc(mono_model[[2]][,which(data_CIK$Name%in%'Apple Inc')[1]]))
plot(mcmc(mono_model[[2]][,which(data_CIK$Name%in%'Apple Inc')[1]]))
data.frame(Estimate=mono_model[[2]][,
              which(data_CIK$Name%in%'Apple Inc')[1]], 
           Iterations=seq(from=1, to=2000))%>%
  ggplot(aes(x = Iterations, 
          y = Estimate)) +
  geom_line() +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size=12))+
  ggtitle("Trace plot of Pr(B|G=0,x) for Apple-2001")

data.frame(Estimate=mono_model[[1]][,
                  which(data_CIK$Name%in%'Apple Inc')[1]], 
           Iterations=seq(from=1, to=2000))%>%
  ggplot(aes(x = Iterations, 
             y = Estimate)) +
  geom_line() +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size=12))+
  ggtitle("Trace plot of Pr(B|G=1,x) for Apple-2001")



acf(mono_model[[1]][,
                    which(data_CIK$Name%in%'Apple Inc')[1]],
    main='ACF Pr(B|G=1,x) Apple 2001')

acf(mono_model[[2]][,
                    which(data_CIK$Name%in%'Apple Inc')[1]],
    main='ACF Pr(B|G=0,x) Apple 2001')
effective_size_BARTpr0 <- c()
effective_size_BARTpr1 <- c()
effective_size_monopr0 <- c()
effective_size_monopr1 <- c()
diag_z_BART_pr1 <- c()
diag_z_BART_pr0 <- c()
diag_z_mono_pr1 <- c()
diag_z_mono_pr0 <- c()
index_list <- c()
set.seed(12296)
post_list=sample.int(dim(mono_model[[1]])[2], 1000, replace=F)
post_list
for (i in 1:dim(mono_model[[1]])[2]){
  #if (i%%50==0)  {
  if (i%in%post_list){
    index_list[[i]] <- i
    BART_mcmc_pr1 <- mcmc(BART_model[[1]][,i])
    BART_mcmc_pr0 <- mcmc(BART_model[[2]][,i])
    #### the monotone way
    mono_mcmc_pr1 <- mcmc(mono_model[[1]][,i])
    mono_mcmc_pr0 <- mcmc(mono_model[[2]][,i])
    
    
    #summary(BART_mcmc_pr0)
   
    effective_size_BARTpr0[[i]] <- effectiveSize(BART_mcmc_pr0)
    effective_size_BARTpr1[[i]] <- effectiveSize(BART_mcmc_pr1)
    effective_size_monopr0[[i]] <- effectiveSize(mono_mcmc_pr0)
    effective_size_monopr1[[i]] <- effectiveSize(mono_mcmc_pr1)
    #plot(mcmc(BART_model[[1]][,100]))
    BART_diag_pr1 <- geweke.diag(BART_mcmc_pr1)
    BART_diag_pr0 <- geweke.diag(BART_mcmc_pr0)
    diag_z_BART_pr1[[i]] <- pnorm(abs(BART_diag_pr1$z),lower.tail=F)*2
    diag_z_BART_pr0[[i]] <- pnorm(abs(BART_diag_pr0$z),lower.tail=F)*2
    mono_diag_pr1 <- geweke.diag(mono_mcmc_pr1)
    mono_diag_pr0 <- geweke.diag(mono_mcmc_pr0)
    diag_z_mono_pr1[[i]] <- pnorm(abs(mono_diag_pr1$z),lower.tail=F)*2
    diag_z_mono_pr0[[i]] <- pnorm(abs(mono_diag_pr0$z),lower.tail=F)*2

    #pnorm(abs(mono_diag$z),lower.tail=F)*2
    #geweke.plot(mcmc(mono_model[[1]][,3]))
  }
  print(i)
}
mono_mcmc_pr1
a0 <- data.frame(prob=mono_mcmc_pr0)%>%
  ggplot(aes(x=var1))+
  geom_histogram(aes(y=..count../sum(..count..)),
          color='white',fill='#1d2951', bins=25)+
  ggtitle('Apple 2001')+
  ylab('Density')+xlab('Pr(B|G=0, x)')+theme_minimal(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5,size=14))
a1 <-data.frame(prob=mono_mcmc_pr1)%>%
  ggplot(aes(x=var1))+
  geom_histogram(aes(y=..count../sum(..count..)),
                 color='white',fill='#1d2951', bins=25)+
  ggtitle('Apple 2001')+
  ylab('Density')+xlab('Pr(B|G=1, x)')+theme_minimal(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5,size=14))
a1

ggsave("C:/Users/demetri/Dropbox/audit_paper/audit_paper_aoas/extraneous_figures/ApplePrBG1hist.pdf", 
       a1, height=6, width=6)
ggsave("C:/Users/demetri/Dropbox/audit_paper/audit_paper_aoas/extraneous_figures/ApplePrBG0hist.pdf", 
       a0, height=6, width=6)
#unlist(index_list)
length(unlist(diag_z_BART_pr1))
diag_df <- data.frame(index_list = unlist(index_list), 
           diag_z_BART_pr1 = unlist(diag_z_BART_pr1),
           diag_z_BART_pr0 = unlist(diag_z_BART_pr0), 
           diag_z_mono_pr1 = unlist(diag_z_mono_pr1),
           diag_z_mono_pr0 = unlist(diag_z_mono_pr0), 
           effective_size_BARTpr0 = unlist(effective_size_BARTpr0), 
           effective_size_BARTpr1 = unlist(effective_size_BARTpr1), 
           effective_size_monopr0 = unlist(effective_size_monopr0), 
           effective_size_monopr1 = unlist(effective_size_monopr1)
           )
diag_df%>%
  summarize_all(mean)

