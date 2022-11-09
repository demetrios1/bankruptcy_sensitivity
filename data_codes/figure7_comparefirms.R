library(readr)
library(data.table)
library(dplyr)
#intframe=read_csv("/home/dpapakos/moderating_variables/mod_var_mbart.csv")
# intframe=data.table::fread(file = "/home/dpapakos/moderating_variables/mod_var_mbart.csv", na.strings = c("", "NA", "#N/A"))
# #First 25350 are B1G1
# intframe_B1G1=intframe[, 2:25351]
# #second 25350 are B1G0
# intframe_B1G0=intframe[, 25352:(25352+25349)]
# intframe_G=intframe[ , 50702:76051]
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
set.seed(12296)
post_list=sample.int(2000, 500, replace=F)
#RR=RR[post_list, ]
EVal=RR+sqrt(RR*(RR-1))

Eval_mean=colMeans(EVal[post_list, ])
mean(Eval_mean)
Eval_mean[post_list]
#get these data from the cik archive on the web
cik_ticker=read_delim('/home/dpapakos/sensitivity_analysis/cik_ticker.csv', 
                      "|", escape_double = FALSE, trim_ws = TRUE)
cik_ticker=cik_ticker%>%dplyr::select(Name, CIK)
data_CIK_fix=read_csv('/home/dpapakos/sensitivity_analysis/data_CIK_fix.csv')
data_CIK_fix

#### updated data ####
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

root_path_post <- '/home/dpapakos/cred_intervals/credinterval_constrained_integral_newvar/'
fill3_B1=data.table::fread(file = paste0(root_path_post, 'mean0sd0.5_credint_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill3_B0=data.table::fread(file = paste0(root_path_post, 'mean0sd0.5_credint_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]
##make this a ratio not a difference
tau_frame=t(fill3_B1/fill3_B0)

tau_frame_2=data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                       bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, tau_frame)

tau_frame_2=right_join(cik_ticker,tau_frame_2)
B0_frame=t(data.frame(fill3_B0))
B1_frame=t(data.frame(fill3_B1))
B0_frame_2=right_join(cik_ticker, data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                                             bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, B0_frame))
B1_frame_2=right_join(cik_ticker, data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                                             bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, B1_frame))
jet_blue_2007_B0=data.frame(B0_frame_2%>%
                              dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504]
jet_blue_2007_B1=data.frame(B1_frame_2%>%
                              dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504]
jet_blue_2009_B0=data.frame(B0_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]
jet_blue_2009_B1=data.frame(B1_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]
jet_blue_2006_B0=data.frame(B0_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
jet_blue_2006_B1=data.frame(B1_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
apple_2001_B0=data.frame(B0_frame_2%>%
                           filter( grepl( "Apple Inc" , Name)))[1, 6:504]
apple_2001_B1=data.frame(B1_frame_2%>%
                           filter( grepl( "Apple Inc" , Name)))[1, 6:504]
build_bear_2010_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]
build_bear_2010_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]


build_bear_2013_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]
build_bear_2013_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]
build_bear_2012_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]
build_bear_2012_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]

build_bear_2014_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]
build_bear_2014_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]

blockbuster_2009_B0=data.frame(B0_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[3, 6:504]
blockbuster_2009_B1=data.frame(B1_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[3, 6:504]

blockbuster_2004_B0=data.frame(B0_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[5, 6:504]
blockbuster_2004_B1=data.frame(B1_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[5, 6:504]

radioshack_2014_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
radioshack_2014_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
radioshack_2013_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]
radioshack_2013_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]

sixflags_2009_B0=data.frame(B0_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[1, 6:504]
sixflags_2009_B1=data.frame(B1_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[1, 6:504]

sixflags_2006_B0=data.frame(B0_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[4, 6:504]

sixflags_2006_B1=data.frame(B1_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[4, 6:504]

sdpt5_frame_B1B0=data.frame(
  firms=c('jetblue_2007', 'jetblue_2009', 'apple_2001', 'build_a_bear_2010', 
          'build_a_bear_2014', 'radioshack_2014', 'blockbuster_2004', 
          'blockbuster_2009', 'sixflags_2006', 'sixflags_2009'),
  B0_mean=c(rowMeans(jet_blue_2007_B0), rowMeans(jet_blue_2009_B0),
            rowMeans(apple_2001_B0), rowMeans(build_bear_2010_B0), rowMeans(build_bear_2014_B0), rowMeans(radioshack_2014_B0), 
            rowMeans(blockbuster_2004_B0), rowMeans(blockbuster_2009_B0), 
            rowMeans(sixflags_2006_B0), rowMeans(sixflags_2009_B0)), 
  B1_mean=c(rowMeans(jet_blue_2007_B1),  rowMeans(jet_blue_2009_B1),
            rowMeans(apple_2001_B1), rowMeans(build_bear_2010_B1), rowMeans(build_bear_2014_B1)
            , rowMeans(radioshack_2014_B1), 
            rowMeans(blockbuster_2004_B1), rowMeans(blockbuster_2009_B1), 
            rowMeans(sixflags_2006_B1), rowMeans(sixflags_2009_B1)))
#write.csv(sdpt5_frame_B1B0,'/home/dpapakos/sensitivity_analysis/sdpt5_frame_B1B0.csv')
write.csv(sdpt5_frame_B1B0,'/home/dpapakos/sensitivity_analysis/individual_results/sdpt5_frame_B1B0.csv')
sdpt5_frame_B1B0



#### observed risk ratio ####
RR_frame=data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                    bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, RR_obs=t(RR[post_list,]))
#View(tau_frame[1:10, 1:10])
RR_frame_2=right_join(cik_ticker,RR_frame)

jet_blue_2007_RR=data.frame(RR_frame_2%>%
                              dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504] 


jet_blue_2009_RR=data.frame(RR_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]
jet_blue_2006_RR=data.frame(RR_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
apple_2001_RR=data.frame(RR_frame_2%>%
                           filter( grepl( "Apple Inc" , Name)))[1, 6:504]
build_bear_2010_RR=data.frame(RR_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]
build_bear_2013_RR=data.frame(RR_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]
build_bear_2012_RR=data.frame(RR_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]
build_bear_2014_RR=data.frame(RR_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]

blockbuster_2009_RR=data.frame(RR_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[3, 6:504]
blockbuster_2004_RR=data.frame(RR_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[5, 6:504]

radioshack_2014_RR=data.frame(RR_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
radioshack_2013_RR=data.frame(RR_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]

sixflags_2009_RR=data.frame(RR_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[1, 6:504]

sixflags_2006_RR=data.frame(RR_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[4, 6:504]


RRobs_frame=data.frame(
  firms=c('jetblue_2007', 'jetblue_2009', 'apple_2001', 'build_a_bear_2010', 
          'build_a_bear_2014','radioshack_2014', 'blockbuster_2004', 
          'blockbuster_2009', 'sixflags_2006', 'sixflags_2009'),
  RR_mean=c(rowMeans(jet_blue_2007_RR), rowMeans(jet_blue_2009_RR),
            rowMeans(apple_2001_RR), rowMeans(build_bear_2010_RR), rowMeans(build_bear_2014_RR), 
            rowMeans(radioshack_2014_RR), 
            rowMeans(blockbuster_2004_RR), rowMeans(blockbuster_2009_RR), 
            rowMeans(sixflags_2006_RR), rowMeans(sixflags_2009_RR)))
RRobs_frame
#write.csv(RRobs_frame,'/home/dpapakos/sensitivity_analysis/individual_results/RRobs.csv')
write.csv(RRobs_frame,'/home/dpapakos/sensitivity_analysis/individual_results/RRobs.csv')
#####


#head(tau_frame_2)
#View(tau_frame_2%>%
#       filter(bank==0&GC==1))
#View(data.frame(tau_frame_2%>%
#             filter( grepl( "Pharm" , Name))))
str(data.frame(tau_frame_2%>%
                 dplyr::filter( grepl( "CVS" , Name)))[, 7:504])
rowMeans(data.frame(tau_frame_2%>%
                      dplyr::filter( grepl( "CVS" , Name)))[, 6:504])
#View(head(data.frame(tau_frame_2%>%
#                      dplyr::filter( grepl( "CVS" , Name)))[, 12:504]))


#### Let's check which of these are NA's now

rowMeans(data.frame(tau_frame_2%>%
                      dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[, 6:504])
jet_blue_2007=data.frame(tau_frame_2%>%
                           dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504]
jet_blue_2009=data.frame(tau_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]

#NOW MISSING!
jet_blue_2006=data.frame(tau_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
apple_2001=data.frame(tau_frame_2%>%
                        filter( grepl( "Apple Inc" , Name)))[1, 6:504]
build_bear_2010=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]
build_bear_2013=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]

build_bear_2012=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]
build_bear_2014=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]

blockbuster_2009=data.frame(tau_frame_2%>%
                              filter( grepl( "Blockbuster" , Name)))[3, 6:504]
blockbuster_2004=data.frame(tau_frame_2%>%
                              filter( grepl( "Blockbuster" , Name)))[5, 6:504]

radioshack_2014=data.frame(tau_frame_2%>%
                             filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
####no longer in!
radioshack_2013=data.frame(tau_frame_2%>%
                             filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]

sixflags_2009=data.frame(tau_frame_2%>%
                           filter( grepl( "Six Flags" , Name)))[1, 6:504]

sixflags_2006=data.frame(tau_frame_2%>%
                           filter( grepl( "Six Flags" , Name)))[4, 6:504]
quantile(radioshack_2014, c(0.025, .5, .975))
quantile(jet_blue_2007, c(0.025, .5, .975))
sdpt5_frame=data.frame(
  firms=c('jetblue_2007', 'jetblue_2009', 'apple_2001', 'build_a_bear_2010', 
          'build_a_bear_2014', 'radioshack_2014', 'blockbuster_2004', 
          'blockbuster_2009', 'sixflags_2006', 'sixflags_2009'),
  post_mean=c(rowMeans(jet_blue_2007),rowMeans(jet_blue_2009),
              rowMeans(apple_2001), rowMeans(build_bear_2010), rowMeans(build_bear_2014), rowMeans(radioshack_2014), 
              rowMeans(blockbuster_2004), rowMeans(blockbuster_2009), 
              rowMeans(sixflags_2006), rowMeans(sixflags_2009)),
  lower=c(
    unlist(quantile(jet_blue_2007, c(0.025, .5, .975))[1]),
    unlist(quantile(jet_blue_2009, c(0.025, .5, .975))[1]),
    unlist(quantile(apple_2001, c(0.025, .5, .975))[1]),
    unlist(quantile(build_bear_2010, c(0.025, .5, .975))[1]),
    unlist(quantile(build_bear_2014, c(0.025, .5, .975))[1]),
    unlist(quantile(radioshack_2014, c(0.025, .5, .975))[1]), 
    # unlist(quantile(radioshack_2013, c(0.025, .5, .975))[1]), 
    unlist(quantile(blockbuster_2004, c(0.025, .5, .975))[1]), 
    unlist(quantile(blockbuster_2009, c(0.025, .5, .975))[1]),
    unlist(quantile(sixflags_2006, c(0.025, .5, .975))[1]), 
    unlist(quantile(sixflags_2009, c(0.025, .5, .975))[1])),
  upper=c(
    unlist(quantile(jet_blue_2007, c(0.025, .5, .975))[3]),
    unlist(quantile(jet_blue_2009, c(0.025, .5, .975))[3]),
    unlist(quantile(apple_2001, c(0.025, .5, .975))[3]),
    unlist(quantile(build_bear_2010, c(0.025, .5, .975))[3]),
    unlist(quantile(build_bear_2014, c(0.025, .5, .975))[3]),
    unlist(quantile(radioshack_2014, c(0.025, .5, .975))[3]), 
    #   unlist(quantile(radioshack_2013, c(0.025, .5, .975))[3]), 
    unlist(quantile(blockbuster_2004, c(0.025, .5, .975))[3]), 
    unlist(quantile(blockbuster_2009, c(0.025, .5, .975))[3]),
    unlist(quantile(sixflags_2006, c(0.025, .5, .975))[3]), 
    unlist(quantile(sixflags_2009, c(0.025, .5, .975))[3]))
)
sdpt5_frame
write.csv(sdpt5_frame,'/home/dpapakos/sensitivity_analysis/individual_results/m0sdpt5individfirm.csv')
library(tidyverse)
dim(apple_2001)
tauframe_apple=t(apple_2001)[1:499,]
tauframe_apple=as.numeric(tauframe_apple)
df_apple=data.frame(tau=tauframe_apple)

par(mfrow=c(2,1))
apple_1=df_apple%>%
  ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),
               color='white',fill='#1d2951', bins=20)+
  ggtitle('Apple 2001')+
  ylab('Density')+xlab('Inducement')+xlim(0,200)+ylim(0,0.32)+theme_minimal(base_size = 16)+
  theme(plot.title = element_text(hjust = 0.5,size=16))
apple_1
tauframe_jetblue=t(jet_blue_2007)[1:499,]
tauframe_jetblue=as.numeric(tauframe_jetblue)
df_jetblue=data.frame(tau=tauframe_jetblue)

jetblue_1=df_jetblue%>%
  ggplot(aes(x=tau))+
  geom_histogram(aes(y=..count../sum(..count..)),color='white',
                 fill='black')+
  ggtitle('Jetblue 2007')+
  ylab('Density')+xlab('Inducement')+
  xlim(0,200)+theme_minimal(base_size = 16)+
  theme(plot.title = element_text(hjust = 0.5,size=16))
library(gridExtra)
grid.arrange(apple_1, jetblue_1,nrow=1)

#write.csv(tau_frame_2,'/home/dpapakos/sensitivity_analysis/constrained_integration/norm0pt5_individual.csv')

#plot(density(as.numeric(radioshack_2014)))
#lines(density(as.numeric(radioshack_2013)))


####Repeat for right bump!, mean=0
library(readr)
library(data.table)

#cik_ticker=read_delim('/home/dpapakos/sensitivity_analysis/cik_ticker.csv', 
#   "|", escape_double = FALSE, trim_ws = TRUE)
#cik_ticker=cik_ticker%>%dplyr::select(Name, CIK)
#data_CIK_fix=read_csv('/home/dpapakos/sensitivity_analysis/data_CIK_fix.csv')
fill3_B1=data.table::fread(file = paste0(root_path_post, 'rightbumpsigma0.48_B1.csv'), na.strings = c("", "NA", "#N/A"))[-1]
fill3_B0=data.table::fread(file = paste0(root_path_post, 'rightbumpsigma0.48_B0.csv'), na.strings = c("", "NA", "#N/A"))[-1]
##make this a ratio not a difference
tau_frame=t(fill3_B1/fill3_B0)

tau_frame_2=data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                       bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, tau_frame)

B0_frame=t(data.frame(fill3_B0))
B1_frame=t(data.frame(fill3_B1))
B0_frame_2=right_join(cik_ticker, data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                                             bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, B0_frame))
B1_frame_2=right_join(cik_ticker, data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                                             bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, B1_frame))
jet_blue_2007_B0=data.frame(B0_frame_2%>%
                              dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504]
jet_blue_2007_B1=data.frame(B1_frame_2%>%
                              dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504]
jet_blue_2009_B0=data.frame(B0_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]
jet_blue_2009_B1=data.frame(B1_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]
jet_blue_2006_B0=data.frame(B0_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
jet_blue_2006_B1=data.frame(B1_frame_2%>%
                              filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
apple_2001_B0=data.frame(B0_frame_2%>%
                           filter( grepl( "Apple Inc" , Name)))[1, 6:504]
apple_2001_B1=data.frame(B1_frame_2%>%
                           filter( grepl( "Apple Inc" , Name)))[1, 6:504]
build_bear_2010_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]
build_bear_2010_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]


build_bear_2013_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]
build_bear_2013_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]
build_bear_2012_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]
build_bear_2012_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]

build_bear_2014_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]
build_bear_2014_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]

blockbuster_2009_B0=data.frame(B0_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[3, 6:504]
blockbuster_2009_B1=data.frame(B1_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[3, 6:504]

blockbuster_2004_B0=data.frame(B0_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[5, 6:504]
blockbuster_2004_B1=data.frame(B1_frame_2%>%
                                 filter( grepl( "Blockbuster" , Name)))[5, 6:504]

radioshack_2014_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
radioshack_2014_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
radioshack_2013_B0=data.frame(B0_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]
radioshack_2013_B1=data.frame(B1_frame_2%>%
                                filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]

sixflags_2009_B0=data.frame(B0_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[1, 6:504]
sixflags_2009_B1=data.frame(B1_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[1, 6:504]

sixflags_2006_B0=data.frame(B0_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[4, 6:504]

sixflags_2006_B1=data.frame(B1_frame_2%>%
                              filter( grepl( "Six Flags" , Name)))[4, 6:504]

sdpt5_frame_B1B0=data.frame(
  firms=c('jetblue_2007', 'jetblue_2009', 'apple_2001', 'build_a_bear_2010', 
          'build_a_bear_2014', 'radioshack_2014', 'blockbuster_2004', 
          'blockbuster_2009', 'sixflags_2006', 'sixflags_2009'),
  B0_mean=c(rowMeans(jet_blue_2007_B0), rowMeans(jet_blue_2009_B0),
            rowMeans(apple_2001_B0), rowMeans(build_bear_2010_B0), rowMeans(build_bear_2014_B0), 
            rowMeans(radioshack_2014_B0), 
            rowMeans(blockbuster_2004_B0), rowMeans(blockbuster_2009_B0), 
            rowMeans(sixflags_2006_B0), rowMeans(sixflags_2009_B0)), 
  B1_mean=c(rowMeans(jet_blue_2007_B1), rowMeans(jet_blue_2009_B1),
            rowMeans(apple_2001_B1), rowMeans(build_bear_2010_B1), rowMeans(build_bear_2014_B1), 
            rowMeans(radioshack_2014_B1), 
            rowMeans(blockbuster_2004_B1), rowMeans(blockbuster_2009_B1), 
            rowMeans(sixflags_2006_B1), rowMeans(sixflags_2009_B1)))

sdpt5_frame_B1B0
#write.csv(sdpt5_frame_B1B0, '/home/dpapakos/sensitivity_analysis/rightbumpindividfirmB1B0.csv')
write.csv(sdpt5_frame_B1B0, '/home/dpapakos/sensitivity_analysis/individual_results/rightbumpindividfirmB1B0.csv')
library(tidyverse)

#head(tau_frame_2)
#View(tau_frame_2%>%
#       filter(bank==0&GC==1))
#View(data.frame(tau_frame_2%>%
#             filter( grepl( "Pharm" , Name))))
#str(data.frame(tau_frame_2%>%
#                 dplyr::filter( grepl( "CVS" , Name)))[, 7:504])
tau_frame=t(fill3_B1/fill3_B0)

tau_frame_2=data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                       bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, tau_frame)

tau_frame_2=right_join(cik_ticker,tau_frame_2)
rowMeans(data.frame(tau_frame_2%>%
                      dplyr::filter( grepl( "CVS" , Name)))[, 6:504])
#View(head(data.frame(tau_frame_2%>%
#                      dplyr::filter( grepl( "CVS" , Name)))[, 12:504]))
rowMeans(data.frame(tau_frame_2%>%
                      dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[, 6:504])
jet_blue_2007=data.frame(tau_frame_2%>%
                           dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504]
jet_blue_2009=data.frame(tau_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]
jet_blue_2006=data.frame(tau_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
apple_2001=data.frame(tau_frame_2%>%
                        filter( grepl( "Apple Inc" , Name)))[1, 6:504]
build_bear_2010=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]
build_bear_2013=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]
build_bear_2012=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]
build_bear_2014=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]

blockbuster_2009=data.frame(tau_frame_2%>%
                              filter( grepl( "Blockbuster" , Name)))[3, 6:504]
blockbuster_2004=data.frame(tau_frame_2%>%
                              filter( grepl( "Blockbuster" , Name)))[5, 6:504]

radioshack_2014=data.frame(tau_frame_2%>%
                             filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
radioshack_2013=data.frame(tau_frame_2%>%
                             filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]

sixflags_2009=data.frame(tau_frame_2%>%
                           filter( grepl( "Six Flags" , Name)))[1, 6:504]

sixflags_2006=data.frame(tau_frame_2%>%
                           filter( grepl( "Six Flags" , Name)))[4, 6:504]
quantile(radioshack_2014, c(0.025, .5, .975))
quantile(jet_blue_2007, c(0.025, .5, .975))
sdpt5_frame=data.frame(
  firms=c('jetblue_2007', 'jetblue_2009', 'apple_2001', 'build_a_bear_2010', 
          'build_a_bear_2014', 'radioshack_2014', 'blockbuster_2004', 
          'blockbuster_2009', 'sixflags_2006', 'sixflags_2009'),
  post_mean=c(rowMeans(jet_blue_2007), rowMeans(jet_blue_2009),
              rowMeans(apple_2001), rowMeans(build_bear_2010), rowMeans(build_bear_2014), 
              rowMeans(radioshack_2014),  
              rowMeans(blockbuster_2004), rowMeans(blockbuster_2009), 
              rowMeans(sixflags_2006), rowMeans(sixflags_2009)),
  lower=c(
    unlist(quantile(jet_blue_2007, c(0.025, .5, .975))[1]),
    unlist(quantile(jet_blue_2009, c(0.025, .5, .975))[1]),
    unlist(quantile(apple_2001, c(0.025, .5, .975))[1]),
    unlist(quantile(build_bear_2010, c(0.025, .5, .975))[1]),
    unlist(quantile(build_bear_2014, c(0.025, .5, .975))[1]),
    #unlist(quantile(radioshack_2013, c(0.025, .5, .975))[1]), 
    unlist(quantile(radioshack_2014, c(0.025, .5, .975))[1]), 
    unlist(quantile(blockbuster_2004, c(0.025, .5, .975))[1]), 
    unlist(quantile(blockbuster_2009, c(0.025, .5, .975))[1]),
    unlist(quantile(sixflags_2006, c(0.025, .5, .975))[1]), 
    unlist(quantile(sixflags_2009, c(0.025, .5, .975))[1])),
  upper=c(
    unlist(quantile(jet_blue_2007, c(0.025, .5, .975))[3]),
    unlist(quantile(jet_blue_2009, c(0.025, .5, .975))[3]),
    unlist(quantile(apple_2001, c(0.025, .5, .975))[3]),
    unlist(quantile(build_bear_2010, c(0.025, .5, .975))[3]),
    unlist(quantile(build_bear_2014, c(0.025, .5, .975))[3]),
    # unlist(quantile(radioshack_2013, c(0.025, .5, .975))[3]), 
    unlist(quantile(radioshack_2014, c(0.025, .5, .975))[3]), 
    unlist(quantile(blockbuster_2004, c(0.025, .5, .975))[3]), 
    unlist(quantile(blockbuster_2009, c(0.025, .5, .975))[3]),
    unlist(quantile(sixflags_2006, c(0.025, .5, .975))[3]), 
    unlist(quantile(sixflags_2009, c(0.025, .5, .975))[3]))
)
sdpt5_frame
#write.csv(sdpt5_frame, '/home/dpapakos/sensitivity_analysis/rightbumpindividfirm.csv')
write.csv(sdpt5_frame, '/home/dpapakos/sensitivity_analysis/individual_results/rightbumpindividfirm.csv')
library(tidyverse)
dim(apple_2001)
tauframe_apple=t(apple_2001)[1:499,]
tauframe_apple=as.numeric(tauframe_apple)
df_apple=data.frame(tau=tauframe_apple)

par(mfrow=c(2,1))
apple_1=df_apple%>%
  ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),
                color='white',fill='#1d2951', bins=25)+
  ggtitle('Apple 2001')+
  ylab('Density')+xlab('Inducement')+xlim(0,100)+ylim(0,0.3)+theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=16))

tauframe_radioshack=t(radioshack_2014)[1:499,]
tauframe_radioshack=as.numeric(tauframe_radioshack)
df_radioshack=data.frame(tau=tauframe_radioshack)

radioshack_1=df_radioshack%>%
  ggplot(aes(x=tau))+
  geom_histogram(aes(y=..count../sum(..count..)),color='white',
                 fill='#1d2951', bins=25)+
  ggtitle('Radioshack 2014')+
  ylab('Density')+xlab('Inducement')+
  xlim(0,100)+theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=16))
library(gridExtra) 
grid.arrange(apple_1, radioshack_1,nrow=1)
ggsave("/home/dpapakos/sensitivity_analysis/Apple_radioshack_right_RR_navy.pdf", 
       grid.arrange(apple_1, radioshack_1,nrow=1), height=5, width=10 )



library(tidyverse)

tau_frame=t(fill3_B1-fill3_B0)

tau_frame_2=data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                       bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, tau_frame)

tau_frame_2=right_join(cik_ticker,tau_frame_2)
rowMeans(data.frame(tau_frame_2%>%
                      dplyr::filter( grepl( "CVS" , Name)))[, 6:504])
#View(head(data.frame(tau_frame_2%>%
#                      dplyr::filter( grepl( "CVS" , Name)))[, 12:504]))
rowMeans(data.frame(tau_frame_2%>%
                      dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[, 6:504])
jet_blue_2007=data.frame(tau_frame_2%>%
                           dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504]
jet_blue_2009=data.frame(tau_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]
jet_blue_2006=data.frame(tau_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
apple_2001=data.frame(tau_frame_2%>%
                        filter( grepl( "Apple Inc" , Name)))[1, 6:504]
build_bear_2010=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]
build_bear_2013=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]
build_bear_2012=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]
build_bear_2014=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]

blockbuster_2009=data.frame(tau_frame_2%>%
                              filter( grepl( "Blockbuster" , Name)))[3, 6:504]
blockbuster_2004=data.frame(tau_frame_2%>%
                              filter( grepl( "Blockbuster" , Name)))[5, 6:504]

radioshack_2014=data.frame(tau_frame_2%>%
                             filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
radioshack_2013=data.frame(tau_frame_2%>%
                             filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]

sixflags_2009=data.frame(tau_frame_2%>%
                           filter( grepl( "Six Flags" , Name)))[1, 6:504]

sixflags_2006=data.frame(tau_frame_2%>%
                           filter( grepl( "Six Flags" , Name)))[4, 6:504]
quantile(radioshack_2014, c(0.025, .5, .975))
quantile(jet_blue_2007, c(0.025, .5, .975))
sdpt5_frame=data.frame(
  firms=c('jetblue_2007', 'jetblue_2009', 'apple_2001', 'build_a_bear_2010', 
          'build_a_bear_2014', 'radioshack_2014', 'blockbuster_2004', 
          'blockbuster_2009', 'sixflags_2006', 'sixflags_2009'),
  post_mean=c(rowMeans(jet_blue_2007), rowMeans(jet_blue_2009),
              rowMeans(apple_2001), rowMeans(build_bear_2010), rowMeans(build_bear_2014), 
              rowMeans(radioshack_2014),  
              rowMeans(blockbuster_2004), rowMeans(blockbuster_2009), 
              rowMeans(sixflags_2006), rowMeans(sixflags_2009)),
  lower=c(
    unlist(quantile(jet_blue_2007, c(0.025, .5, .975))[1]),
    unlist(quantile(jet_blue_2009, c(0.025, .5, .975))[1]),
    unlist(quantile(apple_2001, c(0.025, .5, .975))[1]),
    unlist(quantile(build_bear_2010, c(0.025, .5, .975))[1]),
    unlist(quantile(build_bear_2014, c(0.025, .5, .975))[1]),
    #unlist(quantile(radioshack_2013, c(0.025, .5, .975))[1]), 
    unlist(quantile(radioshack_2014, c(0.025, .5, .975))[1]), 
    unlist(quantile(blockbuster_2004, c(0.025, .5, .975))[1]), 
    unlist(quantile(blockbuster_2009, c(0.025, .5, .975))[1]),
    unlist(quantile(sixflags_2006, c(0.025, .5, .975))[1]), 
    unlist(quantile(sixflags_2009, c(0.025, .5, .975))[1])),
  upper=c(
    unlist(quantile(jet_blue_2007, c(0.025, .5, .975))[3]),
    unlist(quantile(jet_blue_2009, c(0.025, .5, .975))[3]),
    unlist(quantile(apple_2001, c(0.025, .5, .975))[3]),
    unlist(quantile(build_bear_2010, c(0.025, .5, .975))[3]),
    unlist(quantile(build_bear_2014, c(0.025, .5, .975))[3]),
    # unlist(quantile(radioshack_2013, c(0.025, .5, .975))[3]), 
    unlist(quantile(radioshack_2014, c(0.025, .5, .975))[3]), 
    unlist(quantile(blockbuster_2004, c(0.025, .5, .975))[3]), 
    unlist(quantile(blockbuster_2009, c(0.025, .5, .975))[3]),
    unlist(quantile(sixflags_2006, c(0.025, .5, .975))[3]), 
    unlist(quantile(sixflags_2009, c(0.025, .5, .975))[3]))
)

#write.csv(sdpt5_frame, '/home/dpapakos/sensitivity_analysis/rightbumpindividfirm.csv')
#write.csv(sdpt5_frame, '/home/dpapakos/sensitivity_analysis/individual_results/rightbumpindividfirm.csv')
library(tidyverse)
dim(apple_2001)
tauframe_apple=t(apple_2001)[1:499,]
tauframe_apple=as.numeric(tauframe_apple)
df_apple=data.frame(tau=tauframe_apple)

dim(apple_2001)
tauframe_apple=t(apple_2001)[1:499,]
tauframe_apple=as.numeric(tauframe_apple)
df_apple=data.frame(tau=tauframe_apple)

par(mfrow=c(2,1))
apple_1=df_apple%>%
  ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),
                                    color='white',fill='#1d2951', bins=25)+
  ggtitle('Apple 2001')+
  ylab('Density')+xlab('Risk Difference')+scale_x_continuous(limits = 
                              c(0, max(df_apple$tau)-0.1),
                      breaks = round(seq(0.0, max(df_apple$tau),
                                         length.out = 10), 
                    2))+  scale_y_continuous(limits = 
                                               c(0, 0.35,
                                                 breaks = round(seq(0.0, 0.40,
                                                                    length.out = 10), 
                                                                2)))+
  theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=16))
apple_1
tauframe_radioshack=t(radioshack_2014)[1:499,]
tauframe_radioshack=as.numeric(tauframe_radioshack)
df_radioshack=data.frame(tau=tauframe_radioshack)

radioshack_1=df_radioshack%>%
  ggplot(aes(x=tau))+
  geom_histogram(aes(y=..count../sum(..count..)),color='white',
                 fill='#1d2951', bins=25)+
  ggtitle('Radioshack 2014')+
  ylab('Density')+xlab('Risk Difference')+
  scale_x_continuous(limits = 
                       c(0, max(df_apple$tau)-0.1),
                     breaks = round(seq(0.0, max(df_apple$tau),
                                        length.out = 10), 
                                    2))+
  scale_y_continuous(limits = 
                       c(0, 0.40,
                     breaks = round(seq(0.0, 0.40,
                                        length.out = 10), 
                                    2)))+theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=16))
radioshack_1
library(gridExtra) 
grid.arrange(apple_1, radioshack_1,nrow=1)
ggsave("/home/dpapakos/sensitivity_analysis/Apple_radioshack_right_treat_navy.pdf", 
       grid.arrange(apple_1, radioshack_1,nrow=1), height=5, width=10 )

####TREATMENT EFFECTS!####
tau_frame=t(fill3_B1-fill3_B0)

tau_frame_2=data.frame(CIK=data_CIK_fix$CIK,GC=data_CIK_fix$going_concern,
                       bank=data_CIK_fix$bankrptobs, date=data_CIK_fix$sig_date_of_op_s, tau_frame)

tau_frame_2=right_join(cik_ticker,tau_frame_2)
rowMeans(data.frame(tau_frame_2%>%
                      dplyr::filter( grepl( "CVS" , Name)))[, 6:504])
#View(head(data.frame(tau_frame_2%>%
#                      dplyr::filter( grepl( "CVS" , Name)))[, 12:504]))
rowMeans(data.frame(tau_frame_2%>%
                      dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[, 6:504])
jet_blue_2007=data.frame(tau_frame_2%>%
                           dplyr::filter( grepl( "Jetblue Airways Corp" , Name)))[1, 6:504]
jet_blue_2009=data.frame(tau_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[2, 6:504]
jet_blue_2006=data.frame(tau_frame_2%>%
                           filter( grepl( "Jetblue Airways Corp" , Name)))[3, 6:504]
apple_2001=data.frame(tau_frame_2%>%
                        filter( grepl( "Apple Inc" , Name)))[1, 6:504]
build_bear_2010=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[1, 6:504]
build_bear_2013=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[2, 6:504]
build_bear_2012=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[3, 6:504]
build_bear_2014=data.frame(tau_frame_2%>%
                             filter( grepl( "Build A Bear Workshop Inc" , Name)))[4, 6:504]

blockbuster_2009=data.frame(tau_frame_2%>%
                              filter( grepl( "Blockbuster" , Name)))[3, 6:504]
blockbuster_2004=data.frame(tau_frame_2%>%
                              filter( grepl( "Blockbuster" , Name)))[5, 6:504]

radioshack_2014=data.frame(tau_frame_2%>%
                             filter( grepl( "Radioshack Corp" , Name)))[1, 6:504]
radioshack_2013=data.frame(tau_frame_2%>%
                             filter( grepl( "Radioshack Corp" , Name)))[2, 6:504]

sixflags_2009=data.frame(tau_frame_2%>%
                           filter( grepl( "Six Flags" , Name)))[1, 6:504]

sixflags_2006=data.frame(tau_frame_2%>%
                           filter( grepl( "Six Flags" , Name)))[4, 6:504]
quantile(radioshack_2014, c(0.025, .5, .975))
quantile(jet_blue_2007, c(0.025, .5, .975))
sdpt5_frame=data.frame(
  firms=c('jetblue_2007', 'apple_2001', 'build_a_bear_2010', 
          'build_a_bear_2014', 'radioshack_2014', 'blockbuster_2004', 
          'blockbuster_2009', 'sixflags_2006', 'sixflags_2009'),
  post_mean=c(rowMeans(jet_blue_2007), 
              #rowMeans(jet_blue_2009), 
              rowMeans(apple_2001), rowMeans(build_bear_2010), rowMeans(build_bear_2014), 
              rowMeans(radioshack_2014), 
              #rowMeans(radioshack_2013), 
              rowMeans(blockbuster_2004), rowMeans(blockbuster_2009), 
              rowMeans(sixflags_2006), rowMeans(sixflags_2009)),
  lower=c(
    unlist(quantile(jet_blue_2007, c(0.025, .5, .975))[1]),
    # unlist(quantile(jet_blue_2009, c(0.025, .5, .975))[1]),
    unlist(quantile(apple_2001, c(0.025, .5, .975))[1]),
    unlist(quantile(build_bear_2010, c(0.025, .5, .975))[1]),
    unlist(quantile(build_bear_2014, c(0.025, .5, .975))[1]),
    #   unlist(quantile(radioshack_2013, c(0.025, .5, .975))[1]), 
    unlist(quantile(radioshack_2014, c(0.025, .5, .975))[1]), 
    unlist(quantile(blockbuster_2004, c(0.025, .5, .975))[1]), 
    unlist(quantile(blockbuster_2009, c(0.025, .5, .975))[1]),
    unlist(quantile(sixflags_2006, c(0.025, .5, .975))[1]), 
    unlist(quantile(sixflags_2009, c(0.025, .5, .975))[1])),
  upper=c(
    unlist(quantile(jet_blue_2007, c(0.025, .5, .975))[3]),
    # unlist(quantile(jet_blue_2009, c(0.025, .5, .975))[3]),
    unlist(quantile(apple_2001, c(0.025, .5, .975))[3]),
    unlist(quantile(build_bear_2010, c(0.025, .5, .975))[3]),
    unlist(quantile(build_bear_2014, c(0.025, .5, .975))[3]),
    # unlist(quantile(radioshack_2013, c(0.025, .5, .975))[3]), 
    unlist(quantile(radioshack_2014, c(0.025, .5, .975))[3]), 
    unlist(quantile(blockbuster_2004, c(0.025, .5, .975))[3]), 
    unlist(quantile(blockbuster_2009, c(0.025, .5, .975))[3]),
    unlist(quantile(sixflags_2006, c(0.025, .5, .975))[3]), 
    unlist(quantile(sixflags_2009, c(0.025, .5, .975))[3]))
)
sdpt5_frame
#write.csv(sdpt5_frame, '/home/dpapakos/sensitivity_analysis/rightbumpindividfirm_treat.csv')
write.csv(sdpt5_frame, '/home/dpapakos/sensitivity_analysis/individual_results/rightbumpindividfirm_treat.csv')
library(tidyverse)
dim(apple_2001)
tauframe_apple=t(apple_2001)[1:499,]
tauframe_apple=as.numeric(tauframe_apple)
df_apple=data.frame(tau=tauframe_apple)
par(mfrow=c(2,1))
apple_1=df_apple%>%
  ggplot(aes(x=tau))+geom_histogram(aes(y=..count../sum(..count..)),
         color='white',fill='#1d2951', bins=22)+
  ggtitle('Apple 2001')+
  ylab('Density')+xlab('Risk Difference')+xlim(0,0.3)+ylim(0,0.5)+theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=16))

tauframe_jetblue=t(jet_blue_2007)[1:499,]
tauframe_jetblue=as.numeric(tauframe_jetblue)
df_jetblue=data.frame(tau=tauframe_jetblue)

jetblue_1=df_jetblue%>%
  ggplot(aes(x=tau))+
  geom_histogram(aes(y=..count../sum(..count..)),color='white',
                 fill='#1d2951', bins=22)+
  ggtitle('Jetblue 2007')+
  ylab('Density')+xlab('Risk Difference')+
  xlim(0,.3)+ylim(0,0.50)+theme_minimal(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=16))
library(gridExtra)
grid.arrange(apple_1, jetblue_1,nrow=1)

#write.csv(tau_frame_2,'/home/dpapakos/sensitivity_analysis/constrained_integration/norm0pt5_individual.csv')

#plot(density(as.numeric(radioshack_2014)))
#lines(density(as.numeric(radioshack_2013)))

