library(readr)
library(data.table)
#intframe=read_csv("/home/dpapakos/moderating_variables/mod_var_mbart.csv")
#intframe=data.table::fread(file = "/home/dpapakos/moderating_variables/mod_var_mbart.csv", na.strings = c("", "NA", "#N/A"))
####the new variables
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
RR=intframe_B1G1/intframe_B1G0
head(RR[,1:10])
dim(RR)
set.seed(12296)
post_list=sample.int(2000, 500, replace=F)
#RR=RR[post_list, ]

#RR_lower=log(RR)+qnorm(0.025)*
EVal=RR+sqrt(RR*(RR-1))
quantile(Eval_mean, c(0.025, .5, .975))
mean()
Eval_mean=colMeans(EVal[post_list, ])
mean(Eval_mean)

post_list=sample.int(2000, 500, replace=F)
Eval_mean[post_list]
#fill3_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.75sig1.25_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
#fill3_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.75sig1.25_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]
root_path<- '/home/dpapakos/cred_intervals/credinterval_constrained_integral_newvar/'
fill3_B1=data.table::fread(file = paste0(root_path, 'sharkfinq0.75sig1.25_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill3_B0=data.table::fread(file = paste0(root_path, 'sharkfinq0.75sig1.25_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]

ITE_sharkq75s125=colMeans(fill3_B1[, ]-fill3_B0[, ])
ratio_sharkq75s125=colMeans(fill3_B1[, ]/fill3_B0[, ])



I2=ITE_sharkq75s125[which(prop_score_mean<0.05&prop_score_mean>0.005)]
E2=Eval_mean[which(prop_score_mean<0.05&prop_score_mean>0.005)]
prop_score_mean2=prop_score_mean[which(prop_score_mean<0.05&prop_score_mean>0.005)]
#plot(log(E2), I2/prop_score_mean2)



#fill_sharkfinq0.25sig0.5_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.25sig0.5_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
#fill_sharkfinq0.25sig0.5_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.25sig0.5_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

fill_sharkfinq0.25sig0.5_B1=data.table::fread(file = paste0(root_path, 'sharkfinq0.25sig0.5_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill_sharkfinq0.25sig0.5_B0=data.table::fread(file = paste0(root_path, 'sharkfinq0.25sig0.5_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]


ITE_sharkq25s5=colMeans(fill_sharkfinq0.25sig0.5_B1[, ]-fill_sharkfinq0.25sig0.5_B0[, ])
ratio_sharkq25s5=colMeans(fill_sharkfinq0.25sig0.5_B1[, ]/fill_sharkfinq0.25sig0.5_B0[, ])


prop_score_mean2=prop_score_mean[which(prop_score_mean<0.05&prop_score_mean>0.005)]
#plot(log(E2), I2/prop_score_mean2)




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
plot(log(Eval_mean), ITE_sharkq75s125,xlab='Mean E-value (log scale)',
     ylab='Mean posterior ITE',
     main=expression(paste('Sharkfin q=0.75, ', sigma,'=1.25')), cex=.75,col=ifelse(prop_score_mean>0.1, "firebrick4",
                                                  ifelse(prop_score_mean<0.005, 'dodgerblue4', "black")))


plot(log(Eval_mean), ITE_sharkq25s5,xlab='Mean E-value (log scale)',
     ylab='Mean posterior ITE',
     main=expression(paste('Sharkfin q=0.25, ', sigma,'=0.5')), cex=.75,col=ifelse(prop_score_mean>0.1, "firebrick4",
                                                 ifelse(prop_score_mean<0.005, 'dodgerblue4', "black")))


plot(log(Eval_mean), ITE_normsd1, xlab='Mean E-value (log scale)',
     ylab='Mean posterior ITE',
     main='f(u) ~ N(0, 1)', cex=.75, col=ifelse(prop_score_mean>0.1, "firebrick4",
                                         ifelse(prop_score_mean<0.005, 'dodgerblue4', "black")))
plot(log(Eval_mean), ITE_normpt1, cex=.75, xlab='Mean E-value (log scale)',
     ylab='Mean posterior ITE', main='f(u) ~ N(0, 0.1)', col=ifelse(prop_score_mean>0.1, "firebrick4",
                                        ifelse(prop_score_mean<0.005, 'dodgerblue4', "black")))


dev.off()


par(mfrow=c(2,2))

plot(Eval_mean, ratio_sharkq75s125,xlab='Mean E-value',
     ylab='Mean posterior Ratio',
     main=expression(paste('Sharkfin q=0.75,s=1.25 ', sigma, '=0.88'), cex=0.75))#,col=ifelse(prop_score_mean>0.1, "firebrick4",
                                                           #ifelse(prop_score_mean<0.005, 'dodgerblue4', "black")))


plot(Eval_mean, ratio_sharkq25s5,xlab='Mean E-value',
     ylab='Mean posterior Ratio',
     main=expression(paste('Sharkfin q=0.25,s=0.5 ', sigma, '=1.05')), cex=0.75)#,col=ifelse(prop_score_mean>0.1, "firebrick4",
                                                          #ifelse(prop_score_mean<0.005, 'dodgerblue4', "black")))


plot(Eval_mean, ratio_normsd1, xlab='Mean E-value',
     ylab='Mean posterior Ratio',
     main='f(u) ~ N(0, 1)', cex=0.75)#, col=ifelse(prop_score_mean>0.1, "firebrick4",
                                              #   ifelse(prop_score_mean<0.005, 'dodgerblue4', "black")))
plot(Eval_mean, ratio_normpt1, cex=0.75, xlab='Mean E-value',
     ylab='Mean posterior Ratio', main='f(u) ~ N(0, 0.1)')#, col=ifelse(prop_score_mean>0.1, "firebrick4",
                                                                            # ifelse(prop_score_mean<0.005, 'dodgerblue4', "black")))

#### save in cred_intervals as Eval_vs_U_ratio_constrained_newvar_png 600 x 600




data.frame(run=c('sharkq.25sig.5', 'sharkq.75sig1.25', 'N(0, 1)', 'N(0, .1'),
           cors=c(cor(Eval_mean,ratio_sharkq25s5), 
           cor(Eval_mean,ratio_sharkq75s125), 
           cor(Eval_mean, ratio_normsd1), 
           cor(Eval_mean, ratio_normpt1)))
dev.off()
plot(log(Eval_mean), log(ratio_normpt1),  cex=0.25)
abline(0,1, col='red')

fill3_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.75sig1.25_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
fill3_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.75sig1.25_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

shark25_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.25sig0.5_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
shark25_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/sharkfinq0.25sig0.5_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

sdpt1_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd0.1_credint_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
sdpt1_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd0.1_credint_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

sdpt5_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd0.5_credint_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
sdpt5_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd0.5_credint_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

sd1_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd1_credint_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
sd1_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/mean0sd1_credint_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

sym98_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/symmetric98percmid_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
sym98_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/symmetric98percmid_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

sym90_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/symmetric90percmid_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
sym90_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/symmetric90percmid_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]

rightbump_B1=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/rightbumpsigma0.05_B1.csv', na.strings = c("", "NA", "#N/A"))[-1]
rightbump_B0=data.table::fread(file = '/home/dpapakos/cred_intervals/credinterval_constrained_integral/rightbumpsigma0.05_B0.csv', na.strings = c("", "NA", "#N/A"))[-1]




quantile(colMeans(fill3_B1/fill3_B0), c(0.025, .5, .975))
mean(colMeans(fill3_B1/fill3_B0))
quantile(colMeans(shark25_B1/shark25_B0), c(0.025, .5, .975))
mean(colMeans(shark25_B1/shark25_B0))
quantile(colMeans(sdpt1_B1/sdpt1_B0), c(0.025, .5, .975))
mean(colMeans(sdpt1_B1/sdpt1_B0))
quantile(colMeans(sdpt5_B1/sdpt5_B0), c(0.025, .5, .975))
mean(colMeans(sdpt5_B1/sdpt5_B0))
quantile(colMeans(sd1_B1/sd1_B0), c(0.025, .5, .975))
mean(colMeans(sd1_B1/sd1_B0))
quantile(colMeans(sym98_B1/sym98_B0), c(0.025, .5, .975))
mean(colMeans(sym98_B1/sym98_B0))
quantile(colMeans(sym90_B1/sym90_B0), c(0.025, .5, .975))
mean(colMeans(sym90_B1/sym90_B0))
quantile(colMeans(rightbump_B1/rightbump_B0), c(0.025, .5, .975))
mean(colMeans(rightbump_B1/rightbump_B0))
quantile(colMeans(fill3_B1/fill3_B0), c(0.025, .5, .975))
tot_inducement=c(
                 colMeans(sdpt1_B1/sdpt1_B0), 
                 colMeans(sdpt5_B1/sdpt5_B0), 
                 colMeans(sd1_B1/sd1_B0), 
                 colMeans(shark25_B1/shark25_B0),
                 colMeans(fill3_B1/fill3_B0) 
)
 bump_inducement=c(                colMeans(sym98_B1/sym98_B0), 
                 colMeans(sym90_B1/sym90_B0),
                 colMeans(rightbump_B1/rightbump_B0))
tot_inducement                

NN=length(colMeans(sym90_B0))
library(viridis)   
library(ggsci)
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp1<-c('black',
        'firebrick4',
        'darkgreen',
        'darkorchid',
        'dodgerblue4')
cbp2<-c('black', 'dodgerblue4', 'firebrick4' )
data.frame(tot_inducement, distribution=c(rep('Normal 0.1', NN),rep('Normal 0.5',NN),
                                   rep('Normal 1',NN), rep('Shark q=0.25',NN), 
                                   rep('Shark q=0.75',NN)))%>%
ggplot(aes(tot_inducement, colour=distribution)) +xlim(.7,10)+
        #geom_histogram(aes(y=..count../sum(..count..), stat='identity'))+
        geom_density(aes(y=..count../sum(..count..)),alpha=0., lwd=.75) +
        ylab('Density')+
        xlab('Posterior Inducement')+
        #scale_color_viridis(discrete = TRUE, option = "D")+
        #scale_fill_viridis(discrete = TRUE)+
        #scale_color_jco()+
        #scale_fill_jco()+
        scale_color_manual(values = cbp1)+
        theme_minimal(base_size = 12)
data.frame(bump_inducement, distribution=c(rep('98% Peak',NN), 
                              rep('90% Peak',NN), rep('Right Bump',NN)))%>%
        ggplot(aes(bump_inducement, colour=distribution)) +xlim(.7,10)+
        #geom_histogram(aes(y=..count../sum(..count..), stat='identity'))+
        geom_density(aes(y=..count../sum(..count..)),alpha=0., lwd=.75) +ylab('Density')+
        xlab('Posterior Inducement')+
        #scale_color_jco()+
        #scale_fill_jco()+
        scale_color_manual(values = cbp2)+
        theme_minimal(base_size = 12)
plot(density(colMeans(shark25_B1/shark25_B0)), col='dodgerblue4',main='',
     xlim=c(0, 10),
     lwd=2, xlab='Posterior Inducement')
lines(density(colMeans(fill3_B1/fill3_B0)), col='firebrick4', lwd=2)
lines(density(colMeans(sdpt1_B1/sdpt1_B0)), col='darkorchid4', lwd=2)
lines(density(colMeans(sdpt5_B1/sdpt5_B0)), col='firebrick4', lwd=2)
lines(density(colMeans(sd1_B1/sd1_B0)), col='darkorchid4', lwd=2)
lines(density(colMeans(sym98_B1/sym98_B0)), col='firebrick4', lwd=2)
lines(density(colMeans(sym90_B1/sym90_B0)), col='darkorchid4', lwd=2)
lines(density(colMeans(rightbump_B1/rightbump_B0)), col='black', lwd=2)
legend("topright", legend=c(expression(paste("f" [1],"(u)")), expression(paste("f" [2],"(u)")), 
                            expression(paste("f" [3],"(u)")), expression(paste("f" [4],"(u)"))),
       col=c("dodgerblue4", "firebrick4", 'darkorchid4', 'black'), lty=1, cex=1.1)

