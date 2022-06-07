library(pROC)
library(dplyr)
library(foreign)
library(tidyverse)
library(foreach)
library(bcf)
library(monbart)
#library(fastbart)
library(dbarts)
library(foreach)
library(doParallel)

rm(list = ls())
# load data

setwd("/home/dpapakos/moderating_variables/")
data      = read.dta( "AuditorFinal20160419.dta" )
#bartstuff=read.csv('RawScores4_monotone_bart_lgc_020160907.csv')




#### updated data ####
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
setwd("/home/dpapakos/moderating_variables/")

##Try to calculate how accurate auditor is, compare vs bart and random forest##
#merging in the going concern=1, no going concern
data_sasub2<-data_sasub %>%
  #  dplyr::select(year,bankrptobs,going_concern,logassets_r,lev_r,investments_r,cash_r,roa_r,logprice_r,Intangible_r,RD_r,
  #                RDmissing,fracnonfees_r,feemissing,NoRate,RateC,numyears,downgrade,
  #                ey,dt,kpmg,pwc,gt,bdo)%>%
  # group_by(year,going_concern, bankrptobs)%>%
  mutate(nogc=year<=year&going_concern==0,
         gc=year<=year&going_concern==1)




logiclist<-rep(0<1,length(data_sasub$year))  #any column works
#want a logical return, all trues, 
ynogc<-data_sasub2$nogc  #for x predictors, our list of logical vectors and no going concern

ygc<-data_sasub2$gc  # to get length full, going concern and year

xcomponent<-data_sasub[ynogc,var0]


xcomponent2<-data_sasub[ygc,var0]

#xcomponent3<-data_sasub[,var0]
xcomponent3<-data_sasub[logiclist,var0] #predicting going concern.  #not sure
#why need the logical passed in for the rows


#y-components, ie the response

ycomponent<-as.factor(data_sasub[ynogc, 'bankrptobs'])  #what were predicting

ycomponent2=data_sasub[, 'bankrptobs']

ycomponent2<-as.factor(data_sasub[ygc, 'bankrptobs'])
class(ycomponent2)
#ycomponent3<-as.factor(data_sasub[, 'going_concern'])
ycomponent3<-as.factor(data_sasub[logiclist, 'going_concern'])


#Need this so we can predict on whole 1000 list
#xpred<-data_sasub[, var0]  #provided the whole data frame with the exogeneous variables selected
xpred<-data_sasub[logiclist,var0]
# monotone_bart = function(y, z, x, xpred, nskip=5000, ndpost=5000, m = 50, n=N) {
#   
#   sort_ix = order(z, y)
#   x = x[sort_ix,]
#   z = z[sort_ix]
#   y = y[sort_ix]
#   
#   n0 = sum(z==0)
#   n00 = sum(z==0 & y==0)
#   yobs = y
#   yobs0 = rbinom(n, 1, 0.5)
#   yobs[1:n00] = rbinom(n00, 1, 0.5) # These are DA variables
#   yobs0[1:n00] = rbinom(n00, 1, 0.5) # These are DA variables
#   yobs0[1:n00][yobs[1:n00]==1] = 0 # To statisfy DA constraints
#   yobs0[z==0 & y==1] = 1
#   yobs0 = yobs0[1:n0] 
#   
#   offset =  0#qnorm(mean(y[z==1]))
#   offset0 = 0# qnorm(mean(flu$grp==0 & flu$fluy2==1)/mean(flu$grp==1 & flu$fluy2==1)) #<- wtf
#   
#   zz = offset + 3*yobs - 3*(1-yobs)
#   z0 = offset0 + 3*yobs0 - 3*(1-yobs0)
#   
#   ################################################################################
#   # MCMC
#   ################################################################################
#   
#   set.seed(1022)
#   
#   xi = lapply(1:ncol(x), function(i) bcf:::.cp_quantile(x[,i]))
#   fit.mono = bartRcppMono(yobs, zz, t(as.matrix(x)), t(xpred),
#                           yobs0, z0, t(as.matrix(x)),t(xpred),
#                           n00,
#                           xi,
#                           nskip, ndpost, m, 3.0,
#                           offset, offset0)
#   
#   #xpred.exp = rbind(data.frame(xpred, z=1), data.frame(xpred, z=0))
#   #fit =  bart(cbind(x, z), y, xpred.exp,
#   #               nskip=5000, ndpost=5000,
#   #               ntree=50, usequants=T)
#   
#   
#   # P(Y|Z=1, X)
#   pr1 = pnorm(fit.mono$postpred)
#   
#   # P(Y|Z=0, X)
#   pr0 = pr1*pnorm(fit.mono$postpred0)
#   
#   return(list(pr0 = pr0, pr1 = pr1))
#   
# }
var0<-vars
N=dim(data_sasub)[1]
mbart_pred_func=function(y_here, x_here, y_here_2, x_here_2){
  fold=5
  
  test_index <- purrr::map(1:fold,function(i, y_1, y_0){c(y_1[[i]], y_0[[i]])},
                           y_1=split(sample(which(y_here==1)), 1:fold),
                           y_0=split(sample(which(y_here==0)), 1:fold))
  
  test_index_2 <- purrr::map(1:fold,function(i, y_1, y_0){c(y_1[[i]], y_0[[i]])},
                             y_1=split(sample(which(y_here_2==1)), 1:fold),
                             y_0=split(sample(which(y_here_2==0)), 1:fold))
  
  
  y_train <- list()
  y_test <- list()
  
  
  mbart_train_pred <- list()
  mbart_test_pred <- list()
  
  mbart_train_pred_prob <- list()
  mbart_test_pred_prob <- list()
  pred1<-list()
  pred2<-list()
  
  #### CV: filter with ROC in parallel way
  
  covtreat=data_sasub%>%filter(going_concern==1)  #xcomponent
  covcontrol=data_sasub%>%filter(going_concern==0) #xcomponent2
  covtreat$treat=as.factor(covtreat$going_concern) #xcomponent3
  y_train0<-c()
  y_train1<-c()
  y_test0<-c()
  y_test1<-c()
  true_test<-c()
  xtest<-c()
  imp_frame<-list()
  for (cv in 1:fold) {
    
    train_x_0 = t(x_here)[var0,-test_index[[cv]] ]
    test_x_0 = t(x_here)[var0,test_index[[cv]] ]
    train_y_0 = y_here[-test_index[[cv]]]
    test_y_0 = y_here[test_index[[cv]]]
    train_x_1 = t(x_here_2)[var0,-test_index_2[[cv]] ]
    test_x_1 = t(x_here_2)[var0,test_index_2[[cv]] ]
    train_y_1 = y_here_2[-test_index_2[[cv]]]
    test_y_1 = y_here_2[test_index_2[[cv]]]
    x_train_0=t(train_x_0)
    x_test_0=t(test_x_0)
    x_train_1=t(train_x_1)
    x_test_1=t(test_x_1)
    y_train0[[cv]] = train_y_0
    y_test0[[cv]] = test_y_0
    y_train1[[cv]] = train_y_1
    y_test1[[cv]] = test_y_1
    ytrain0    = as.factor( covcontrol$bankrptobs )
    ytrain1    = as.factor( covtreat$bankrptobs )
    
    
    xtest[[cv]]=data_sasub[setdiff(test_index[[cv]],test_index_2[[cv]]), var0]
    
    
    bart_mono = monotone_bart(y = as.numeric(c(y_train1[[cv]], y_train0[[cv]])==1),
                              z = 1-c(rep(1, length(y_train1[[cv]])), rep(0, length(y_train0[[cv]]))),
                              x = rbind(x_train_1, x_train_0),
                              xpred = xtest[[cv]], nskip = 2000, ndpost = 2000,m=100)
    pred1[[cv]]=1-colMeans(bart_mono$pr0)
    
    pred2[[cv]]= 1-colMeans(bart_mono$pr1)
    
    true_test[[cv]]=data_sasub[setdiff(test_index[[cv]],test_index_2[[cv]]), c('bankrptobs', 'going_concern')]
    imp_frame[[cv]]=data.frame(true_test[[cv]], BG1=pred1[[cv]],  BG0=pred2[[cv]])
  }
  
  print(cv)
  return(list(true_test = unlist(true_test),
              true_train_0=unlist(y_train0),
              true_test_0=unlist(y_test0),
              true_train_1=unlist(y_train1),
              true_test_1=unlist(y_test1),
              pred1=unlist(pred1), 
              pred2=unlist(pred2), 
              imp_frame=sapply(1:5, function(cv)list(imp_frame[[cv]]))
  ))
}

x_here=xcomponent
x_here_2=xcomponent2
y_here=ycomponent
y_here_2=ycomponent2
MBART_CV=mbart_pred_func(y_here=ycomponent, x_here=xcomponent, y_here_2=ycomponent2, x_here_2=xcomponent2)
imp_frame=rbind(MBART_CV$imp_frame[[1]],MBART_CV$imp_frame[[2]],
      MBART_CV$imp_frame[[2]], MBART_CV$imp_frame[[4]], MBART_CV$imp_frame[[5]])
setwd('/home/dpapakos/ROC_plots/')
#write.csv(imp_frame, 'mBART_roc.csv')
#imp_frame=read.csv('mBART_roc.csv')
BG1_mbart=auc(imp_frame$bankrptobs[imp_frame$going_concern==1],imp_frame$BG1[imp_frame$going_concern==1])
BG0_mbart=auc(imp_frame$bankrptobs[imp_frame$going_concern==0],imp_frame$BG0[imp_frame$going_concern==0])
BG1_mbart


rf_fit_func=function(y_here, x_here){
  fold=5
  test_index <- purrr::map(1:fold,function(i, y_1, y_0){c(y_1[[i]], y_0[[i]])},
                           y_1=split(sample(which(y_here==1)), 1:fold),
                           y_0=split(sample(which(y_here==0)), 1:fold))
  
  
  
  y_train <- list()
  y_test <- list()
  
  
  rf_train_pred <- list()
  rf_test_pred <- list()
  
  rf_train_pred_prob <- list()
  rf_test_pred_prob <- list()
  
  
  #### CV: filter with ROC in parallel way
  
  
  test_index[[1]]
  for (cv in 1:fold) {
    train_x = t(x_here)[,-test_index[[cv]] ]
    test_x = t(x_here)[,test_index[[cv]] ]
    train_y = y_here[-test_index[[cv]]]
    test_y = y_here[test_index[[cv]]]
    x_train=t(train_x)
    x_test=t(test_x)
    y_train[[cv]] = train_y
    y_test[[cv]] = test_y
    cv_train_n = length(train_y)
    rf_fit = randomForest::randomForest(x=x_train[,], y=factor(y_train[[cv]], levels=c(0,1)),
                                        xtest=x_test[, ], keep.forest=T)
    rf_train_pred_prob[[cv]]=predict( rf_fit, x_train, type='prob')[, 2] 
    rf_test_pred_prob[[cv]]=predict( rf_fit, x_test, type='prob')[, 2] 
    rf_train_pred[[cv]] = as.numeric(rf_fit$predicted) - 1
    rf_test_pred[[cv]] = as.numeric(rf_fit$test$predicted) - 1
  }
  return(
    list(true_train = unlist(y_train),
         true_test = unlist(y_test),
         rf_train_prob = unlist(rf_train_pred_prob),
         rf_test_prob = unlist(rf_test_pred_prob),
         rf_train = unlist(rf_train_pred),
         rf_test = unlist(rf_test_pred)
    ))
}
rando1=rf_fit_func(ycomponent, xcomponent)
rando2=rf_fit_func(ycomponent2, xcomponent2)
rando3=rf_fit_func(ycomponent3, xcomponent3)
bart_fit<-function(y_here, x_here){
  fold=5
  test_index <- purrr::map(1:fold,function(i, y_1, y_0){c(y_1[[i]], y_0[[i]])},
                           y_1=split(sample(which(y_here==1)), 1:fold),
                           y_0=split(sample(which(y_here==0)), 1:fold))
  
  
  
  y_train <- list()
  y_test <- list()
  
  
  bart_train_pred <- list()
  bart_test_pred <- list()
  
  bart_train_pred_prob <- list()
  bart_test_pred_prob <- list()
  
  
  #### CV: filter with ROC in parallel way
  

  for (cv in 1:fold) {
    train_x = t(x_here)[,-test_index[[cv]] ]
    test_x = t(x_here)[,test_index[[cv]] ]
    train_y = y_here[-test_index[[cv]]]
    test_y = y_here[test_index[[cv]]]
    x_train=t(train_x)
    x_test=t(test_x)
    y_train[[cv]] = train_y
    y_test[[cv]] = test_y
    cv_train_n = length(train_y)
    
    bart_fit=bart(x_train, as.factor(y_train[[cv]]),x_test,ndpost = 2000, nskip = 2000,ntree=100,
                  verbose=T,usequants = TRUE,numcut = 1000)
    
    
    bart_train_pred_prob[[cv]]=colMeans(pnorm(bart_fit$yhat.train))
    bart_test_pred_prob[[cv]]=colMeans(pnorm(bart_fit$yhat.test)) 
    bart_train_pred[[cv]] = 1#as.numeric(bart_fit$predicted) - 1
    bart_test_pred[[cv]] = 1#as.numeric(bart_fit$test$predicted) - 1
    print(cv)
  }
  
  return(
    list(true_train = unlist(y_train),
         true_test = unlist(y_test),
         bart_train_prob = unlist(bart_train_pred_prob),
         bart_test_prob = unlist(bart_test_pred_prob),
         bart_train = unlist(bart_train_pred),
         bart_test = unlist(bart_test_pred)
    ))
  
}
getwd()
setwd("/home/dpapakos/ROC_plots")

bart3=bart_fit(ycomponent3, xcomponent3)
#write.csv(bart3, 'bart3.csv')
auc(bart3$true_test, bart3$bart_test_prob)
boost_fit<-function(y_here, x_here){
  fold=5
  test_index <- purrr::map(1:fold,function(i, y_1, y_0){c(y_1[[i]], y_0[[i]])},
                           y_1=split(sample(which(y_here==1)), 1:fold),
                           y_0=split(sample(which(y_here==0)), 1:fold))
  
  
  
  y_train <- list()
  y_test <- list()
  
  
  boost_train_pred <- list()
  boost_test_pred <- list()
  
  boost_train_pred_prob <- list()
  boost_test_pred_prob <- list()
  
  
  #### CV: filter with ROC in parallel way
  
  
  test_index[[1]]
  for (cv in 1:fold) {
    train_x = t(x_here)[,-test_index[[cv]] ]
    test_x = t(x_here)[,test_index[[cv]] ]
    train_y = y_here[-test_index[[cv]]]
    test_y = y_here[test_index[[cv]]]
    x_train=t(train_x)
    x_test=t(test_x)
    y_train[[cv]] = train_y
    y_test[[cv]] = test_y
    cv_train_n = length(train_y)
    boost_fit<-xgboost::xgboost(data=as.matrix(x_train),label=as.matrix(y_train[[cv]]),
                                max.depth = 4, eta = 1,
                                nthread = 2, nrounds = 5,
                                objective = "binary:logistic")
    
    boost_train_pred_prob[[cv]]=predict( boost_fit, x_train, type='prob')
    boost_test_pred_prob[[cv]]=predict( boost_fit, x_test, type='prob') 
    boost_train_pred[[cv]] = as.numeric(boost_fit$predicted) - 1
    boost_test_pred[[cv]] = as.numeric(boost_fit$test$predicted) - 1
    print(cv)
  }
  
  return(
    list(true_train = unlist(y_train),
         true_test = unlist(y_test),
         boost_train_prob = unlist(boost_train_pred_prob),
         boost_test_prob = unlist(boost_test_pred_prob),
         boost_train = unlist(boost_train_pred),
         boost_test = unlist(boost_test_pred)
    ))
  
}
dim(xcomponent2)
library(xgboost)

boost1=boost_fit(ycomponent, xcomponent)
boost2=boost_fit(ycomponent2, xcomponent2)
boost3=boost_fit(ycomponent3, xcomponent3)
auc(boost1$true_test, boost1$boost_test_prob)
auc(boost2$true_test, boost2$boost_test_prob)
auc(boost3$true_test, boost3$boost_test_prob)
#rando1=rf_fit_func(ycomponent, xcomponent)
#rando2=rf_fit_func(ycomponent2, xcomponent2)
#rando3=rf_fit_func(ycomponent3, xcomponent3)
auc(rando1$true_test, rando1$rf_test_prob)
auc(rando2$true_test, rando2$rf_test_prob)
auc(rando3$true_test, rando3$rf_test_prob)
library(pROC)


# rando1<-randomForest(x=xcomponent,y=ycomponent )
# rando1
# 
# rando2<-randomForest(x=xcomponent2,y=ycomponent2)
# rando3<-randomForest(x=xcomponent3,y=ycomponent3)
# predrando1
# 
# #predict the values
# predrando1  <- rando1$votes#predict( rando1, xpred, type='vote', norm=T) 
# 
# predrando2<-rando2$votes#predict(rando2,xpred, type='prob')
# predrando3<-rando3$votes#predict(rando3,xpred, type='response', norm=T, predict.all=TRUE)
# rando3$votes
# table(predrando3[,2 ]>0.5, dd_2$going_concern)
# rf1=auc(dd_2$going_concern, rando3$votes[, 2])
# 
# 
# rf2=auc(dd_2$bankrptobs[dd_2$going_concern==1], 
#         rando2$votes[,2])#rf$bprob_gc1[dd_2$going_concern==1])
# 
# plot(roc(dd_2$bankrptobs[dd_2$going_concern==1], 
#          rando2$votes[,2]))#rf$bprob_gc1[dd_2$going_concern==1]))
# 
# 
# 
# rf3=auc(dd_2$bankrptobs[dd_2$going_concern==0], 
#         rando1$votes[,2])
# rf3
#rf$bprob_gc0[dd_2$going_concern==0])
#plot(roc(data_sasub2$going_concern, bartstuff$gprob))

#boost=read_csv('outcomeboost.csv')
# colnames(boost)
# boost1auc=auc(data_sasub2$going_concern, boost1$prob_g)
# boost2auc=auc(data_sasub2$bankrptobs[data_sasub2$going_concern==1], 
#            boost$prob_b1_g1[data_sasub2$going_concern==1])
# boost3auc=auc(data_sasub2$bankrptobs[data_sasub2$going_concern==0], 
#            boost$prob_b1_g0[data_sasub2$going_concern==0])

#verification::roc.plot(data_sasub2$going_concern, bartstuff$gprob)
data.frame(BART=c(bart1, bart2, bart3), RF=c(rf1,rf2, rf3), BOOST=c(boost1, boost2, boost3))



########################################################################
######################## 1. initialization 
#rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(pROC)
library(caret)
library(colorspace)



#' Functions plots multiple 'roc' objects into one plot
#' @param rocs
#'   A list of 'roc' objects. Every list item has a name.
#' @param breaks
#'   A vector of integers representing ticks on the x- and y-axis
#' @param legentTitel
#'   A string which is used as legend titel
ggrocs <- function(rocs, breaks = seq(0,1,0.1), tittle = "Fit Model") {
  if (length(rocs) == 0) {
    stop("No ROC objects available in param rocs.")
  } else {
    require(plyr)
    # Store all sensitivities and specifivities in a data frame
    # which an be used in ggplot
    RocVals <- plyr::ldply(names(rocs), function(rocName) {
      if(class(rocs[[rocName]]) != "roc") {
        stop("Please provide roc object from pROC package")
      }
      data.frame(
        fpr = rev(1-rocs[[rocName]]$specificities),
        tpr = rev(rocs[[rocName]]$sensitivities),
        auc = rep(sprintf("%.3f",rocs[[rocName]]$auc), length(rocs[[rocName]]$sensitivities)),
        trials = rep(rocName, length(rocs[[rocName]]$sensitivities)),
        stringsAsFactors = FALSE
      )
    })
    
    
    rocPlot <- ggplot(RocVals, aes(x = fpr, y = tpr, colour = trials)) +
      scale_color_manual(values=c('dodgerblue4', 'firebrick4', 'darkgreen'))+
      #scale_colour_viridis(discrete = TRUE)+
      geom_line(size = 1.5, alpha = 0.4) + 
      geom_segment(aes(x = 0, y = 0, xend = 1,yend = 1), alpha = 0.4, colour = "gray") + 
      geom_step() +
      scale_x_continuous(name = "False Positive Rate (1 - Specificity)",limits = c(0,1), breaks = breaks) + 
      scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,1), breaks = breaks) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))+
      
      ggtitle(tittle) + 
      guides(fill = guide_legend(title=NULL)) + 
      coord_equal() 
    
    
     auc_table <- unique(RocVals[, c("trials", "auc")])
      df.table <- gridExtra::tableGrob(auc_table, theme = gridExtra::ttheme_default(base_size = 8), rows = NULL)
    # Plot chart and table into one object
    return(grid.arrange(rocPlot,# df.table,
                        ncol=2,
                        as.table=TRUE,
                        widths=c(8,1)))
  }
}
library(viridis)

aa=ggrocs(list(BART=roc(bart3$true_test, bart3$bart_test_prob), 
               BOOST=roc(boost3$true_test, boost3$boost_test_prob), 
               RF=roc(rando3$true_test, rando3$rf_test_prob)), 
          breaks = seq(0,1,0.1), 
          tittle = "ROC curve Pr(G|x)")

bb=ggrocs(list(BART=roc(imp_frame$bankrptobs[imp_frame$going_concern==1],imp_frame$BG1[imp_frame$going_concern==1]), 
               BOOST=roc(boost2$true_test, boost2$boost_test_prob), 
               RF=roc(rando2$true_test, rando2$rf_test_prob)), 
          breaks = seq(0,1,0.1), 
          tittle = "ROC curve Pr(B|G=1, x)")

ggsave("/home/dpapakos/ROC_plots/roc_pbg1_cv_png.png", bb, height=4, width=6)

data_sasub2$bankrptobs[data_sasub2$going_concern==0]
cc=ggrocs(list(BART=roc(data_sasub2$bankrptobs[data_sasub2$going_concern==0], 
                        imp_frame$BG0[data_sasub2$going_concern==0]), 
               BOOST=roc(data_sasub2$bankrptobs[data_sasub2$going_concern==0], 
                         boost1$boost_test_prob[data_sasub2$going_concern==0]), 
               RF=roc(data_sasub2$bankrptobs[data_sasub2$going_concern==0], rando1$rf_test_prob)), 
          breaks = seq(0,1,0.1), 
          tittle = "ROC curve Pr(B|G=0, x)")

#full train AUC
intframe2=data.table::fread("/home/dpapakos/moderating_variables/mbart_audit_newvars.csv")
bart_BG0AUCtrain <- pROC::auc(ycomponent, intframe2$notreatpred[data_sasub2$going_concern==0])
bart_BG1AUCtrain <- pROC::auc(ycomponent2, intframe2$treatpred[data_sasub2$going_concern==1])
bart_BG0AUCtrain

bart_GAUCtrain <- pROC::auc(ycomponent3, intframe2$propensity)
bart_BG1AUCtrain
bart_GAUCtrain
