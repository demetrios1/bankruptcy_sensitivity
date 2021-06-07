
library(dplyr)
library(foreign)
library(EValue)
setwd("C:/Users/demetri/Documents/ASU/RTG/monotone_bart/auditor_stuff")
rm(list = ls())
# load data

data      = read.dta( "AuditorFinal20160419.dta" )
View(head(data))
length(unique(data$sic_code_descrip))
bartstuff=read.csv('RawScores4_monotone_bart_lgc_020160907.csv')
#bartstuff=read.csv('RawScores4_bart_lgc_020160907.csv')
setwd("C:/Users/demetri/Documents/ASU/RTG/")

#bartstuff=read.csv('outcomemine.csv')
#bartstuff=read.csv('monobartstufflgc_020160907.csv')


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
# var0
# dd_2$bankrptobs=as.factor(dd_2$bankrptobs)



treat.eq <- import ~ coop + cost
out.eq <- export ~ cost + target
f.list <- list(treat.eq, out.eq)
mr <- c("probit", "probit")

library(GJRM)


## CLASSIC RECURSIVE BIVARIATE PROBIT

out <- gjrm(list(as.factor(going_concern) ~ s(logassets_r)+s(lev_r)+s(investments_r)+s(cash_r)+s(roa_r)+
                   s(logprice_r) +s(Intangible_r)+s(RD_r)+RDmissing+fracnonfees_r+feemissing+
                   NoRate+RateC+numyears+downgrade, 
                 as.factor(bankrptobs) ~ going_concern +s(logassets_r)+s(lev_r)+s(investments_r)+s(cash_r)+s(roa_r)+
                   s(logprice_r) +s(Intangible_r)+s(RD_r)+RDmissing+fracnonfees_r+feemissing+
                   NoRate+RateC+numyears+downgrade ), 
            data = dd_2,
            margins = c("probit", "probit"),
            Model = "B")
conv.check(out)                        
summary(out)
AIC(out); BIC(out)


X.int <- as.matrix(out$X2)

ind.int<-(1:out$X2.d2)+out$X1.d2

coef.int <- as.numeric(out$coefficients[ind.int])

d0 <- d1 <- X.int
d0[, 'going_concern'] <- 0
d1[, 'going_concern'] <- 1
eti1 <- d1 %*% coef.int
eti0 <- d0 %*% coef.int
  p.int1 <- probm(eti1, out$margins[2], min.dn = out$VC$min.dn, 
                  min.pr = out$VC$min.pr, max.pr = out$VC$max.pr)$pr
  p.int0 <- probm(eti0, out$margins[2],  min.dn = out$VC$min.dn, 
                  min.pr = out$VC$min.pr, max.pr = out$VC$max.pr)$pr

gamma=qnorm(p.int1)-qnorm(p.int0)
gamma=mean(gamma)
est.AT <- mean(p.int1, na.rm = TRUE) - mean(p.int0, na.rm = TRUE)

n.sim=500
bs <- rMVN(n.sim, mean = out$coefficients, sigma = out$Vb)
eti1s <- d1 %*% t(bs[, ind.int])
eti0s <- d0 %*% t(bs[, ind.int])

peti1s <- probm(eti1s, out$margins[2])$pr
peti0s <- probm(eti0s, out$margins[2])$pr
est.ATb <- colMeans(peti1s, na.rm = TRUE) - colMeans(peti0s, 
                                                     na.rm = TRUE)
mult <- 100
hist(est.ATb * mult, freq = FALSE, main = 'nice', xlab = 'treatment', 
     ylim = c(0, max(density(est.ATb * mult)$y, hist(est.ATb * 
                                                       mult, plot = FALSE)$density)))
lines(density(est.ATb * mult))

ITE=p.int1 - p.int0
rho=summary(out)$theta
rho
##double check using the source code
ate_2<-est.AT

summary(out)[]$theta

beta0=out$coefficients[1]
beta1=out$coefficients[2:80]

beta1=c(beta1[1:7], sum(beta1[8:16]),
sum(beta1[17:25]),sum(beta1[26:34]),sum(beta1[35:43]),sum(beta1[44:52]), 
sum(beta1[53:61]), sum(beta1[62:70]), sum(beta1[71:79]))
alpha0=out$coefficients[81]
alpha1=out$coefficients[83:161]

#out$coefficients[82] is gamma
alpha1=c(alpha1[1:7], sum(alpha1[8:16]),
        sum(alpha1[17:25]),sum(alpha1[26:34]),sum(alpha1[35:43]),sum(alpha1[44:52]), 
        sum(alpha1[53:61]), sum(alpha1[62:70]), sum(alpha1[71:79]))
XX=data_all[, var0]
mean(beta0)
beta0
mean(alpha1)
mean1=beta0+as.matrix(XX)%*%beta1
mean2=alpha0+as.matrix(XX)%*%alpha1

#mean 1 is for going concern
mean(mean1)
#mean2 is for bankruptcy
mean(mean2)

head(dd_2)
#no smoothing
out_2<-gjrm(list(going_concern ~logassets_r+lev_r+investments_r+cash_r+roa_r+
                     logprice_r+Intangible_r+RD_r+as.factor(RDmissing)+fracnonfees_r+as.factor(feemissing)+
                   as.factor(NoRate)+as.factor(RateC)+numyears+as.factor(downgrade), 
           bankrptobs~ going_concern+logassets_r+lev_r+investments_r+cash_r+roa_r+
             logprice_r+Intangible_r+RD_r +as.factor(RDmissing)+fracnonfees_r+as.factor(feemissing)+
             as.factor(NoRate)+as.factor(RateC)+numyears+as.factor(downgrade)), 
           data=dd_2,
           margins = c("probit", "probit"),
           Model = "B"
           )
#logassets_r+lev_r+investments_r+cash_r+roa_r+
#  logprice_r+Intangible_r+RD_r+RDmissing+fracnonfees_r+feemissing+
#  NoRate+RateC+numyears+downgrade
library(pROC)
auc( dd_2$bankrptobs,out$fit$p11)
#

#
## treatment effect, risk ratio and odds ratio with CIs

mb(dd_2$going_concern,dd_2$bankrptobs, Model = "B")
ATE=AT(out, nm.end = "going_concern", hd.plot = TRUE, n.sim = 100)
ATE
sd(ATE$sim.AT)/10
summary(out)
GJRM::RR(out, nm.end = "going_concern") 
GJRM::OR(out, nm.end = "going_concern") 
AT(out_2, nm.end = "going_concern", type = "joint", n.sim=100) 
re.imp <- imputeCounter(out, m = 10, "y1")
re.imp$AT

####calculate risk ratio ####
summary.gjrm(out)
aa<-GJRM::OR(out, 'going_concern')
RR<-GJRM::RR(out, 'going_concern')
library(EValue)
evalues.RR(RR$res[2],RR$res[1], RR$res[3], rare=T)
bias_plot(RR$res[2], xmax=15)
