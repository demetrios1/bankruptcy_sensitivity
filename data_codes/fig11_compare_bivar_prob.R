library(readr)
library(data.table)
library(haven)
getwd()
root_path<- "/home/dpapakos/moderating_variables/"
data      = as.data.frame(read_dta( paste0(root_path, "auditor_with_stratio_biprob_s25678.dta")))

rho = 0.420 
sqrt(rho*(1-rho))
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
subset_bivar <- data_sasub[, c('gvkey', 'datadate', 'prob_b1g1', 'prob_b1g0')]
mean(subset_bivar$prob_b1g1, na.rm=T)/mean(subset_bivar$prob_b1g0, na.rm=T)
bivar_RR <- subset_bivar$prob_b1g1/subset_bivar$prob_b1g0
bivar_RD <- subset_bivar$prob_b1g1-subset_bivar$prob_b1g0
set.seed(12296)

root_path<- '/home/dpapakos/cred_intervals/credinterval_constrained_integral_newvar/'
fill3_B1=data.table::fread(file = paste0(root_path, 'sharkfinq0.75sig1.25_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill3_B0=data.table::fread(file = paste0(root_path, 'sharkfinq0.75sig1.25_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]
#B1 =E(G=1) and B0=E(G=0).  So we are taking the average of E(G=1)/E(G=0), so this is still okay 
ITE_sharkq75s125=colMeans(fill3_B1[, ]-fill3_B0[, ])
ratio_sharkq75s125=colMeans(fill3_B1[, ]-fill3_B0[, ])


#fill_sharkfinq0.25sig0.5_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.25sig0.5_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
#fill_sharkfinq0.25sig0.5_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.25sig0.5_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

fill_sharkfinq0.25sig0.5_B1=data.table::fread(file = paste0(root_path, 'sharkfinq0.25sig0.5_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill_sharkfinq0.25sig0.5_B0=data.table::fread(file = paste0(root_path, 'sharkfinq0.25sig0.5_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]


ITE_sharkq25s5=colMeans(fill_sharkfinq0.25sig0.5_B1[, ]-fill_sharkfinq0.25sig0.5_B0[, ])
ratio_sharkq25s5=colMeans(fill_sharkfinq0.25sig0.5_B1[, ]/fill_sharkfinq0.25sig0.5_B0[, ])



#fill_normsd1_B1<-data.table::fread(file='/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd1_credint_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
#fill_normsd1_B0<-data.table::fread(file='/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd1_credint_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

fill_normsd1_B1<-data.table::fread(file=paste0(root_path, 'mean0sd1_credint_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill_normsd1_B0<-data.table::fread(file=paste0(root_path,'mean0sd1_credint_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]


ITE_normsd1<-colMeans(fill_normsd1_B1[, 1:divs]-fill_normsd1_B0[, 1:divs])
ratio_normsd1=colMeans(fill_normsd1_B1[, 1:divs]/fill_normsd1_B0[, 1:divs])

#fill_normpt1_B1<-data.table::fread(file='/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd0.1_credint_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
#fill_normpt1_B0<-data.table::fread(file='/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd0.1_credint_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]
fill_normpt1_B1<-data.table::fread(file=paste0(root_path,'mean0sd0.1_credint_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill_normpt1_B0<-data.table::fread(file=paste0(root_path,'mean0sd0.1_credint_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]
ITE_normpt1<-colMeans(fill_normpt1_B1[, ]-fill_normpt1_B0[, ])
ratio_normpt1<-colMeans(fill_normpt1_B1[, ]/fill_normpt1_B0[, ])
par(mfrow=c(2,2))
plot(bivar_RD, ITE_sharkq75s125,xlab='Bivariate Probit Risk Difference',
     ylab='Mean Posterior Ratio',
     main=expression(paste('Sharkfin q=0.75,s=1.25 ', sigma, '=0.88'), cex=0.75))


plot(bivar_RD, ITE_sharkq25s5,xlab='Bivariate Probit Risk Difference',
     ylab='Mean Posterior Ratio',
     main=expression(paste('Sharkfin q=0.25,s=0.5 ', sigma, '=1.05')), cex=0.75)

plot(bivar_RD, ITE_normpt1,xlab='Bivariate Probit Risk Difference',
     ylab='Mean Posterior Ratio',
     main='f(u) ~ N(0, 0.1)', cex=0.75)

plot(bivar_RD,ITE_normsd1,xlab='Bivariate Probit Risk Difference',
     ylab='Mean Posterior Risk Differences',
     main='f(u) ~ N(0, 1)', cex=0.75)




fill_normpt5_B1=data.table::fread(file = paste0(root_path, 'mean0sd0.5_credint_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill_normpt5_B0=data.table::fread(file = paste0(root_path, 'mean0sd0.5_credint_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]


ITE_normpt5=colMeans(fill_normpt5_B1[, ]-fill_normpt5_B0[, ])
ratio_normpt5=colMeans(fill_normpt5_B1[, ]/fill_normpt5_B0[, ])


fill_rightbumpsigma0.48_B1=data.table::fread(file = paste0(root_path, 'rightbumpsigma0.48_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill_rightbumpsigma0.48_B0=data.table::fread(file = paste0(root_path, 'rightbumpsigma0.48_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]


ITE_rightbumpsigma0.48=colMeans(fill_rightbumpsigma0.48_B1[, ]-fill_rightbumpsigma0.48_B0[, ])
ratio_rightbumpsigma0.48=colMeans(fill_rightbumpsigma0.48_B1[, ]/fill_rightbumpsigma0.48_B0[, ])




par(mfrow=c(1,2))
plot(bivar_RD, ITE_normpt5,xlab='Bivariate Probit Risk Difference',
     ylab='Mean Posterior Risk Difference',
     main='f(u) ~ N(0, 0.5)', cex=0.75)
#abline(0,1,col='red',lwd=3)

plot(bivar_RD, ITE_rightbumpsigma0.48,xlab='Bivariate Probit Risk Difference',
     ylab='Mean Posterior Risk Difference',
     main=expression(paste('Asymmetric Mixture (', sigma, '=0.49)')), cex=0.75)
#abline(0,1,col='red',lwd=3)

plot(bivar_RR, ratio_normpt5,xlab='Bivariate Probit Risk Ratio',
     ylab='Mean Posterior Ratio',
     main='f(u) ~ N(0, 0.5)', cex=0.75)

plot(bivar_RR,ratio_rightbumpsigma0.48,xlab='Bivariate Probit Risk Ratio',
     ylab='Asymmetric Mixture',
     main=expression(paste('Asymmetric Mixture (', sigma, '=0.49)')), cex=0.75)