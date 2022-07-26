library(dplyr)
library(foreign)
library(tidyverse)
library(foreach)
library(bcf)
library(fastbart)
library(dbarts)
library(foreach)
library(doParallel)

#library(xgboost)

#bartstuff=read.csv('outcomemine.csv')
#bartstuff=read.csv('monobartstufflgc_020160907.csv')


##Now, we make our prediction matrix
  

set.seed(12296)
theme_set(theme_minimal(base_size = 20))

# load data
getwd()
root_path <- "C:/Users/demetri/Dropbox/audit_paper/data/"
data      = read.dta( 
  paste0(root_path,"AuditorFinal20160419.dta" ))
#### updated data ####
library(haven)
data      = as.data.frame(read_dta(paste0(
  root_path,"auditor_with_stratio.dta")))
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

mbart_fit <- data.table::fread(paste0(
  root_path, '/mbart_audit_newvars.csv'
))
bartstuff<-data.frame(data_sasub, 
                      bprob_gc0=mbart_fit$notreatpred,
                      bprob_gc1=mbart_fit$treatpred)
###all firms
#merging in the going concern=1, no going concern
data_all<-bartstuff %>%
  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
                RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade, bprob_gc0, bprob_gc1)%>%
  group_by(year,going_concern)%>%
  mutate(nogc=year<=year&going_concern==0,
         gc=year<=year&going_concern==1)

##Try to calculate how accurate auditor is, compare vs bart and random forest##
#merging in the going concern=1, no going concern
#data_EY<-data_sasub %>%
#  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
#             RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade,
#             ey)%>%
#    group_by(year,going_concern, bankrptobs)%>%
#  mutate(nogc=year<=year&going_concern==0,
#         gc=year<=year&going_concern==1)


##the various firms##
#Ernst&Yound
data_EY<-bartstuff %>%
  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
                RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade,
                ey, bprob_gc1, bprob_gc0)%>%
  group_by(year,going_concern, bankrptobs)%>%
  mutate(nogc=year<=year&going_concern==0,
         gc=year<=year&going_concern==1)%>%
  filter(ey==T)


##KPMG##
data_kpmg<-bartstuff %>%
  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
                RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade,
                kpmg,bprob_gc1, bprob_gc0)%>%
  group_by(year,going_concern, bankrptobs)%>%
  mutate(nogc=year<=year&going_concern==0,
         gc=year<=year&going_concern==1)%>%
  filter(kpmg==T)


##PWC##
data_pwc<-bartstuff %>%
  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
                RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade,
                pwc,bprob_gc1, bprob_gc0)%>%
  group_by(year,going_concern, bankrptobs)%>%
  mutate(nogc=year<=year&going_concern==0,
         gc=year<=year&going_concern==1)%>%
  filter(pwc==T)

##deloitte##
data_dt<-bartstuff %>%
  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
                RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade,
                dt,bprob_gc1, bprob_gc0)%>%
  group_by(year,going_concern, bankrptobs)%>%
  mutate(nogc=year<=year&going_concern==0,
         gc=year<=year&going_concern==1)%>%
  filter(dt==T)

####other####
data_other<-bartstuff %>%
  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
                RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade,
                ey,dt, pwc, kpmg, bprob_gc1, bprob_gc0)%>%
  group_by(year,going_concern, bankrptobs)%>%
  mutate(nogc=year<=year&going_concern==0,
         gc=year<=year&going_concern==1)%>%
  filter( ey==F & dt==F & pwc == F & kpmg==F)

dim(data_all)[1]-(dim(data_EY)[1]+dim(data_kpmg)[1]+dim(data_dt)[1]+dim(data_pwc)[1]+dim(data_other)[1])

alphabetafun<-function(data, gamma){
  B1G1input=data[data$going_concern==1,]
  
  B1G1output=B1G1input$bankrptobs
  
  
  B1G0input=data[data$going_concern==1,]
  B1G0output=B1G0input$bankrptobs
  
  # probB1G1class=randomForest(B1G1input[var0], as.factor(B1G1output))
  # probB1G1=predict(probB1G1class, data[var0], type='prob')
  # 
  # probB1G0class=randomForest(B1G0input[var0], as.factor(B1G0output))
  # probB1G0=predict(probB1G0class, data[var0], type='prob')
  # 
  # 
  # 
  # #column 1 is B=0, column 2 is B=1
  # empty=c()
  # for (i in 1:length(data$year)){
  #   empty[[i]]=matrix(c(probB1G1[,2][i], probB1G0[,2][i], probB1G1[,1][i], probB1G0[,1][i]), nrow=2,ncol=2,byrow=T)
  # }
  # lhsmatrix=do.call(rbind, empty)
  
  
  ####instead of ranfom forest, use the precomputed monotone barts from csv, because monotone bart is SLOW!
  
  probB1G0=matrix(nrow=nrow(data), ncol=2)
  probB1G1=matrix(nrow=nrow(data), ncol=2)
  
  probB1G1[,1]=1-data$bprob_gc1
  probB1G1[,2]=data$bprob_gc1
  probB1G0[,1]=1-data$bprob_gc0
  probB1G0[,2]=data$bprob_gc0
  
  
  
  #column 1 is B=0, column 2 is B=1
  empty=c()
  for (i in 1:length(data$year)){
    empty[[i]]=matrix(c(probB1G1[,2][i], probB1G0[,2][i], probB1G1[,1][i], probB1G0[,1][i]), nrow=2,ncol=2,byrow=T)
  }
  lhsmatrix=do.call(rbind, empty)
  #rhsmatrix=alpha*matrix(c(0,eta,1-eta,0), nrow=2*length(data_EY$year),ncol=2)+beta*matrix(c(gamma, 0,0,1-gamma),nrow=2*length(data_EY$year),ncol=2, byrow=T )+(1-alpha-beta)*matrix(rep(matrix(cbind(c(1,0),c(0,0)),nrow=2, ncol=2),length(data_EY$year)), nrow=2*length(data_EY$year),ncol=2, byrow=T)
  
  pi11=c()
  pi10=c()
  pi01=c()
  pi00=c()
  eta=c()
  #gamma=c()
  alpha=c()
  beta=c()
  #gamma=runif(1)
  #worry about speeding up later!
  
  for (i in 1:length(data$year)){
    pi11[[i]]=empty[[i]][1]
    pi10[[i]]=empty[[i]][3]
    pi01[[i]]=empty[[i]][2]
    pi00[[i]]=empty[[i]][4]
    #eta[[i]]=pi10[i]/(pi01[i]+pi10[i])
    #gamma[[i]]=1-(pi00[i]/(pi11[i]+(pi10[i]^2/(pi01[i]+pi10[i]))-1))
    
    #alpha
    alpha[[i]]=1-(pi11[[i]][1]+pi00[[i]][1])
    
    beta[[i]]=pi00[[i]][1]/(1-gamma)#(pi00[i]-pi11[i]+1-alpha[[i]])/(2*(1-gamma))
  }
  
  tau=1-unlist(alpha)-unlist(beta)
  
 
  return(tau)
}

#eta=pi10/(pi01+pi10)
#data.frame(1-alpha-beta)%>%
#  ggplot(aes(x=`X1...alpha...beta`))+geom_histogram(fill='dodgerblue4', color='black')

#eta=lm(pi10~alpha-1)

#eta=lm(head(alpha)-head(pi01)~head(alpha)-1)#$coefficients

#just generate a random gamma between 0 and 1
#but repeat this many times to get multiple betas, run histogram of treatments


#gamma=runif(1,0, .1)

#tau=alphabetafun(runif(1,0,.1))
tau<-c()
bad_gamma<-c()
gamma_list=seq(from=0.005, to=.35, length.out=100)
tau_mean<-c()
# for ( i in 1:length(gamma_list)){
#   
# 
#   #while(any(tau>0)){
#   
#   #gamma<-runif(1,0, .1)
#   tau[[i]]<-alphabetafun(data_EY,gamma_list[[i]])
#   tau_mean[[i]]<-mean(tau[[i]])
#   if (any (tau[[i]]>1)=='TRUE'){
#     bad_gamma[[i]]<-gamma_list[[i]]
#     print('gamma too big')
#   }
#   else{
#     bad_gamma[[i]]<-'NA'
#   }
#   
# }


gam_const<-data_sa%>%select(bankrptobs, going_concern, pwc, ey, kpmg, dt)%>%
  na.omit()%>%
  filter(pwc==T|ey==T|kpmg==T|dt==T)%>%
  mutate(xx=going_concern)%>%#+bankrptobs)%>%
  summarize(gamma_rat=sum(xx)/n())

gam_const=gam_const$gamma_rat

####ernst and young first ####
ey_mean<-data_sa%>%select(bankrptobs, going_concern, ey)%>%
  filter(ey==T)%>%
  mutate(xx=going_concern)%>%#+bankrptobs)%>%
  summarize(gamma_rat=sum(xx)/n())

ey_value<-alphabetafun(data_EY,ey_mean$gamma_rat)#gam_const)#
tau_ey<-mean(ey_value)
tau_ey
tau_ey_hist<-data.frame(ey_value)%>%
  ggplot(aes(x=ey_value))+geom_histogram(fill='dodgerblue4', color='black')+theme_minimal()+
  xlab('Inducement')+ ggtitle(expression(paste('Inducement Full-ID model: Ernst & Young'))) +xlim(-.2, .5)+
  theme(plot.title = element_text(hjust = 0.5))

####deloitte next ####
dt_mean<-data_sa%>%select(bankrptobs, going_concern, dt)%>%
  filter(dt==T)%>%
  mutate(xx=going_concern)%>%#+bankrptobs)%>%
  summarize(gamma_rat=sum(xx)/n())
dt_value<-alphabetafun(data_dt,dt_mean$gamma_rat)#gam_const)#
tau_dt<-mean(dt_value)

tau_dt_hist<-data.frame(dt_value)%>%
  ggplot(aes(x=dt_value))+geom_histogram(fill='dodgerblue4', color='black')+theme_minimal()+
  xlab('Inducement')+ ggtitle(expression(paste('Inducement Full-ID model: Deloitte'))) +xlim(-.2, .5)+
  theme(plot.title = element_text(hjust = 0.5))


####pwc next ####
pwc_mean<-data_sa%>%select(bankrptobs, going_concern, pwc)%>%
  filter(pwc==T)%>%
  mutate(xx=going_concern)%>%#+bankrptobs)%>%
  summarize(gamma_rat=sum(xx)/n())
pwc_value<-alphabetafun(data_pwc,pwc_mean$gamma_rat)#gam_const)#
tau_pwc<-mean(pwc_value)

tau_pwc_hist<-data.frame(pwc_value)%>%
  ggplot(aes(x=pwc_value))+geom_histogram(fill='dodgerblue4', color='black')+theme_minimal()+
  xlab('Inducement')+ ggtitle(expression(paste('Inducement Full-ID model: PWC'))) +xlim(-.2, .5)+
  theme(plot.title = element_text(hjust = 0.5))

####Finally KPMG####
kpmg_mean<-data_sa%>%select(bankrptobs, going_concern, kpmg)%>%
  filter(kpmg==T)%>%
  mutate(xx=going_concern)%>%#+bankrptobs)%>%
  summarize(gamma_rat=sum(xx)/n())

kpmg_value<-alphabetafun(data_kpmg,kpmg_mean$gamma_rat)#gam_const)#
tau_kpmg<-mean(kpmg_value)

tau_kpmg_hist<-data.frame(kpmg_value)%>%
  ggplot(aes(x=kpmg_value))+geom_histogram(fill='dodgerblue4', color='black')+theme_minimal()+
  xlab('Inducement')+ ggtitle(expression(paste('Inducement Full-ID model: KPMG'))) +xlim(-.2, .5)+
  theme(plot.title = element_text(hjust = 0.5))

####Finally the rest####
other_mean<-data_sa%>%select(bankrptobs, going_concern, kpmg, dt, pwc, ey)%>%
  filter(kpmg==F & dt==F & pwc == F & ey == F)%>%
  mutate(xx=going_concern)%>%#+bankrptobs)%>%
  summarize(gamma_rat=sum(xx)/n())

other_value<-alphabetafun(data_other,other_mean$gamma_rat)#gam_const)#
tau_other<-mean(other_value)

tau_other_hist<-data.frame(other_value)%>%
  ggplot(aes(x=other_value))+geom_histogram(fill='dodgerblue4', color='black')+theme_minimal()+
  xlab('Inducement')+ ggtitle(expression(paste('Inducement Full-ID model: Other Auditors'))) +xlim(-.2, .5)+
  theme(plot.title = element_text(hjust = 0.5))
tau_kpmg_hist



####all next ####
all_mean<-data_sa%>%select(bankrptobs, going_concern)%>%
  na.omit()%>%
  # filter(pwc==T)%>%
  mutate(xx=going_concern)%>%#+bankrptobs)%>%
  summarize(gamma_rat=sum(xx)/n())

#data_all
#sapply(data_all, function(x) sum(is.na(x)))
all_value<-alphabetafun(data_all,all_mean$gamma_rat)#gam_const)#
tau_all<-mean(all_value)

tau_all_hist<-data.frame(all_value)%>%
  ggplot(aes(x=all_value))+geom_histogram(fill='dodgerblue4', color='black')+theme_minimal()+
  xlab('Inducement')+ ggtitle(expression(paste('Inducement Full-ID model: All Auditors'))) +xlim(-.2, .5)+
  theme(plot.title = element_text(hjust = 0.5))


c(kpmg_mean$gamma_rat, ey_mean$gamma_rat, dt_mean$gamma_rat, pwc_mean$gamma_rat, all_mean$gamma_rat, other_mean$gamma_rat)


c(tau_all, tau_pwc, tau_ey, tau_kpmg, tau_dt, tau_other)


#}else{
#  tau<-1-alpha-beta
##  print('not nice')
#}

any(tau<0)
mean(tau)
data.frame(tau)%>%
  ggplot(aes(x=tau))+geom_histogram(fill='dodgerblue4', color='black')
#lhsmatrix


#iterate over different values of gamma, 
#which we are assuming constant.
#beta=(pi00+pi11+1-alpha)/(2*(1-gamma))

beta=unlist(pi00)/(1-gamma)

#beta=(pi11+alpha-1)/(gamma-1)
head(beta)
tau=1-unlist(alpha)-unlist(beta)

sum(tau[tau<0])





##sanity check##
#good
#head(matrix(rep(matrix(cbind(c(1,0),c(0,0)),nrow=2, ncol=2),length(data_EY$year)), nrow=2*length(data_EY$year),ncol=2, byrow=T)) 

#good
#head(2*matrix(c(2, 0,0,1-2),nrow=2*length(data_EY$year),ncol=2, byrow=T ))

#good
#head(2*matrix(c(0,3,1-3,0), nrow=2*length(data_EY$year),ncol=2, byrow=T))
#good
#head(matrix(c(probB1G1[,2][1], probB1G0[,2][1], probB1G1[,1][1], probB1G0[,1][1]), nrow=2*length(data_EY$year),ncol=2,byrow=T))





##Uber Ratings##
 

library(boot)
list2=c(4.56, 4.74, 4.76, 4.81, 4.93, 4.92, 4.75, 4.98,5)
set.seed(12296)
mean.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)
boot.out <- boot(list2$list2, mean.fun, R=10000, sim="ordinary")



bootframe=as.data.frame(boot.out$t)
ecdf(bootframe$V1)(5)
uberrate=bootframe%>%
  ggplot(aes(x=V1))+geom_histogram( aes(y=..count../sum(..count..)),fill='dodgerblue4', color='white')+geom_vline(aes(xintercept=4.56,linetype = "Matt"),lwd=1,show.legend = TRUE)+geom_vline(aes(xintercept=5.0,linetype = "Demetri"),lwd=1,show.legend = TRUE)+geom_vline(aes(xintercept=4.81,linetype = "Eddie"),lwd=1,show.legend = TRUE)+ggtitle('Swole Patrole and Liam Bootstrap')+xlab('Uber Rating')+ylab('Bootstrapped Frequency')+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5))
uberrate 
ggsave('uberrating.png',uberrate)

hist(boot.out$t)
ecdf(boot.out$t)(17230)


list2=data.frame(list2)
list2%>%
  ggplot(aes(x=list2))+geom_density(fill='dodgerblue4', color='white')+theme_minimal()
teststat=(4.56-mean(list2$list2))/(sd(list2$list2))
teststat
1-pt(-teststat, length(list2$list2)-1)



