setwd("C:/Users/demetri/Documents/ASU/RTG/monotone_bart/auditor_stuff")

library(dplyr)
library(foreign)
library(EValue)
library(rattle)
library(partykit)
setwd("C:/Users/demetri/Documents/ASU/RTG/monotone_bart/auditor_stuff")
rm(list = ls())
# load data

data      = read.dta( "AuditorFinal20160419.dta" )
bartstuff=read.csv('RawScores4_monotone_bart_lgc_020160907.csv')
#bartstuff=read.csv('RawScores4_bart_lgc_020160907.csv')
setwd("C:/Users/demetri/Documents/ASU/RTG/")

#bartstuff=read.csv('outcomemine.csv')
#bartstuff=read.csv('monobartstufflgc_020160907.csv')



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

#subset of sample
data_sasub<-data_sa


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
outcome=read_csv("Inducement_lgc0_monotone_mean0sd0.5.csv")
#outcome=read.csv('Inducement_lgc0_monotone_multimodes0.1.csv') 
dataset=data.frame(data_sasub, gamma=outcome$gamma)
rrr<-rpart(gamma ~logassets_r+ lev_r+ investments_r +cash_r+roa_r+logprice_r+Intangible_r+RD_r+
             RDmissing+fracnonfees_r+feemissing+NoRate+RateC+numyears+downgrade+
             ey+dt+kpmg+pwc+gt+bdo,data=dataset)
summary(rrr)
#rrr



rrr$variable.importance
rrr$frame
rpart.plot::rpart.plot(rrr)
party=as.party(rrr)
write.csv(which(rrr$where==18), 'mu0sd0.5largestsub.csv')

library(tidyverse)
unique(rrr$where)
hist(dataset[which(rrr$where==18),'gamma'])
data.frame(gamma=dataset[which(rrr$where==18),'gamma'])%>%
  ggplot(aes(x=gamma))+geom_histogram(aes(y=..count../sum(..count..)),color='white',fill='dodgerblue4')+
  ggtitle(paste0('f(u)~ \u03BC=0, \u03C3=0.5'))+theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+xlab('ATE post., lev_r>0.77, log_assets>-0.5, downgrade=T')+ylab('density')
quantile(dataset[which(rrr$where==18),'gamma'], c(0.05, .5, .95))
rrr$frame
rrr$where==19
library(fmsb)
res <- riskratio(5, 10, 90, 90)
str(res)
print(res)
riskratio(12, 5, 18, 17)
riskratio(12, 5, 18, 17, p.calc.by.independence=FALSE)
sub_group_18=bartstuff$gprob[which(rrr$where==18)]
sub_group_18_gc=bartstuff$going_concern[which(rrr$where==18)]
sub_group_18_bank=bartstuff$bankrptobs[which(rrr$where==18)]


sub_group_19=bartstuff$gprob[which(rrr$where==19)]
sub_group_19_gc=bartstuff$going_concern[which(rrr$where==19)]
sub_group_19_bank=bartstuff$bankrptobs[which(rrr$where==19)]
#3 is the farthest down left
sub_group_3=bartstuff$gprob[which(rrr$where==3)]
sub_group_3_gc=bartstuff$going_concern[which(rrr$where==3)]
sub_group_3_bank=bartstuff$bankrptobs[which(rrr$where==3)]
#second farthest left
sub_group_6=bartstuff$gprob[which(rrr$where==6)]
sub_group_6_gc=bartstuff$going_concern[which(rrr$where==6)]
sub_group_6_bank=bartstuff$bankrptobs[which(rrr$where==6)]


cuts=table(cut(sub_group_18, 5, include.lowest=TRUE))
cuts
bank_size<-c()
Eval<-c()
Eval_low<-c()
Eval_high<-c()
bins=cut(sub_group_18, 5, include.lowest=TRUE, labels=c("1", "2", "3", "4", "5"))

bins

for ( i in 1:5){

  group=sub_group_18[which(bins==i)]
  conf_table=table(prop_score=sub_group_18_gc[which(bins==i)], #group>splits[i], 
                   bankrupcty=factor(sub_group_18_bank[which(bins==i)], levels=c('1', '0')))
  conf_table
  bank_size[[i]]<-conf_table[1]+conf_table[2]
  no_gc=conf_table[1]/(conf_table[1]+conf_table[3])
  yes_gc=conf_table[2]/(conf_table[2]+conf_table[4])
  RR=yes_gc/no_gc
  RR2=riskratio(conf_table[2], conf_table[1], conf_table[2]+conf_table[4], conf_table[1]+conf_table[3])
  RR=RR2$estimate
  RRlow=RR2$conf.int[1]
  RRhigh=RR2$conf.int[2]
  Eval_low[[i]]=ifelse(RRlow>1, RRlow+sqrt(RRlow*(RRlow-1)), 1)
  Eval_high[[i]]=RRhigh+sqrt(RRhigh*(RRhigh-1))
  Eval[[i]]=RR+sqrt(RR*(RR-1))
  Eval
}
RR
Eval
Eval_low
Eval_high
conf_table


bank_size<-c()
Eval<-c()
Eval_low<-c()
Eval_high<-c()
bins=cut(sub_group_19, 5, include.lowest=TRUE, labels=c("1", "2", "3", "4", "5"))
bins

cuts=table(cut(sub_group_19, 5, include.lowest=TRUE))
cuts
rrr$frame

for ( i in 1:5){
  
  group=sub_group_19[which(bins==i)]
  conf_table=table(prop_score=sub_group_19_gc[which(bins==i)], #group>splits[i], 
                   bankrupcty=factor(sub_group_19_bank[which(bins==i)], levels=c('1', '0')))
  conf_table
  bank_size[[i]]<-conf_table[1]+conf_table[2]
  no_gc=conf_table[1]/(conf_table[1]+conf_table[3])
  yes_gc=conf_table[2]/(conf_table[2]+conf_table[4])
  RR=yes_gc/no_gc
  RR2=riskratio(conf_table[2], conf_table[1], conf_table[2]+conf_table[4], conf_table[1]+conf_table[3])
  RR=RR2$estimate
  RRlow=RR2$conf.int[1]
  RRhigh=RR2$conf.int[2]
  Eval_low[[i]]=ifelse(RRlow>1, RRlow+sqrt(RRlow*(RRlow-1)), 1)
  Eval_high[[i]]=RRhigh+sqrt(RRhigh*(RRhigh-1))
  Eval[[i]]=RR+sqrt(RR*(RR-1))
  Eval
}
RR
Eval
Eval_low
Eval_high
conf_table
bank_size<-c()
Eval<-c()
Eval_low<-c()
Eval_high<-c()
bins=cut(sub_group_6, 5, include.lowest=TRUE, labels=c("1", "2", "3", "4", "5"))
bins
cuts=table(cut(sub_group_6, 5, include.lowest=TRUE))
cuts
for ( i in 1:5){
  
  group=sub_group_6[which(bins==i)]
  conf_table=table(prop_score=sub_group_6_gc[which(bins==i)], #group>splits[i], 
                   bankrupcty=factor(sub_group_6_bank[which(bins==i)], levels=c('1', '0')))
  conf_table
  bank_size[[i]]<-conf_table[1]+conf_table[2]
  no_gc=conf_table[1]/(conf_table[1]+conf_table[3])
  yes_gc=conf_table[2]/(conf_table[2]+conf_table[4])
  RR=yes_gc/no_gc
  RR2=riskratio(conf_table[2], conf_table[1], conf_table[2]+conf_table[4], conf_table[1]+conf_table[3])
  RR=RR2$estimate
  RRlow=RR2$conf.int[1]
  RRhigh=RR2$conf.int[2]
  Eval_low[[i]]=ifelse(RRlow>1, RRlow+sqrt(RRlow*(RRlow-1)), 1)
  Eval_high[[i]]=RRhigh+sqrt(RRhigh*(RRhigh-1))
  Eval[[i]]=RR+sqrt(RR*(RR-1))
  Eval
}
RR
Eval
Eval_low
Eval_high
conf_table
bank_size
Eval<-c()
bank_size<-c()
Eval_low<-c()
Eval_high<-c()
bins=cut(sub_group_3, 5, include.lowest=TRUE, labels=c("1", "2", "3", "4", "5"))

cuts
cuts=table(cut(sub_group_3, 5, include.lowest=TRUE))
cuts
for ( i in 1:5){
  
  group=sub_group_3[which(bins==i)]
  conf_table=table(prop_score=sub_group_3_gc[which(bins==i)], #group>splits[i], 
                   bankrupcty=factor(sub_group_3_bank[which(bins==i)], levels=c('1', '0')))
  bank_size[[i]]<-conf_table[1]+conf_table[2]
  no_gc=conf_table[1]/(conf_table[1]+conf_table[3])
  yes_gc=conf_table[2]/(conf_table[2]+conf_table[4])
  RR=yes_gc/no_gc
  RR2=riskratio(conf_table[2], conf_table[1], conf_table[2]+conf_table[4], conf_table[1]+conf_table[3])
  RR=RR2$estimate
  RRlow=RR2$conf.int[1]
  RRhigh=RR2$conf.int[2]
  Eval_low[[i]]=ifelse(RRlow>1, RRlow+sqrt(RRlow*(RRlow-1)), 1)
  Eval_high[[i]]=RRhigh+sqrt(RRhigh*(RRhigh-1))
  Eval[[i]]=RR+sqrt(RR*(RR-1))
  Eval
}
RR
Eval
Eval_low
Eval_high
conf_table




tot_eval=(Eval[[1]]*length(bin1))/dim(data_all)[1]+
  (Eval[[2]]*length(bin2))/dim(data_all)[1]+
  (Eval[[3]]*length(bin3))/dim(data_all)[1]+
  (Eval[[4]]*length(bin4))/dim(data_all)[1]+
  (Eval[[5]]*length(bin5))/dim(data_all)[1]
tot_eval



# length(dataset[which(rrr$where==6), 'gamma'])
# dataset[which(rrr$where==6), 'gamma']
# rrr$frame
# unique(rrr$where)


