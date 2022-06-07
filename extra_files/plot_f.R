#setwd("C:/Users/demetri/Documents/ASU/dissertation")
beta=seq(-4, 4, by=0.001)
sharkfin2=function(q,s, beta){
  ifelse(beta>0, (2*(1-q)*(1/s))*dnorm(beta/s),2*q*(dnorm(beta)))
}

sharkfin = function(z){
  val = rep(NA, length(z))
  val[z < 0] = 2*q*dnorm(z[z<0], sd = sig)
  val[z>0] = 2*(1-q)*dnorm(z[z>0],0,sd = sig*(1-q)/q)
  return(val)
  
}
shark2=function(z){
  val=dcauchy(z, 0,1)
}

q=0.25
s=0.34
s=(1-q)/q
sig=0.50
library(dplyr)
library(ggplot2) 
theme_set(theme_gray(base_size = 16))
theme_set(theme_minimal(base_size = 16))
sharkfinq25=data.frame(u=beta, density=sharkfin(beta))%>%
  
  ggplot(aes(x=u, y=density))+geom_line(color='black',lwd=1.5)+ggtitle(expression(paste('Shark Fin q=0.25, ', 's=0.5')))+
  ylab(paste0('f(u)'))+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5, size=18))
sharkfinq25
q=0.75
s=(1-q)/q
sig=1.25
sharkfinq75=data.frame(u=beta, density=sharkfin(beta))%>%
  ggplot(aes(x=u, y=density))+geom_line(color='black',lwd=1.5)+ggtitle('Shark Fin q=0.75, s=1.25')+
  theme_minimal()+theme(plot.title = element_text(hjust = 0.5,size=20))+xlim(-4, 4)
newf=function(u){
  .01*dnorm(u, mean=-2, sd=0.05)+0.94*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
}
symf1=function(u){
  .01*dnorm(u, mean=-2, sd=0.05)+0.98*dnorm(u,mean=0,sd=0.05)+.01*dnorm(u,mean=2,sd=0.05)
}
symf2=function(u){
  .05*dnorm(u, mean=-2, sd=0.05)+0.9*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
}
NN=10000
dists=runif(NN)
datastuff=vector(length=NN)
probs=c(0.05, .95, 1)
for (i in 1:NN){
  if(dists[i]<probs[1]){
    datastuff[i]=rnorm(1, mean=-2, sd=0.05)
  }else if (dists[i]< probs[2]){
    datastuff[i]=rnorm(1, mean=0, sd=0.05)
  }
  else{
    datastuff[i]= rnorm(1, mean=2, sd=0.05)
  }
}
sd(datastuff)
hist(datastuff, 40, freq=F)
dists2=runif(NN)
datastuff2=vector(length=NN)
probs2=c(0.01, .99, 1)
for (i in 1:NN){
  if(dists2[i]<probs2[1]){
    datastuff2[i]=rnorm(1, mean=-2, sd=0.05)
  }else if (dists2[i]< probs2[2]){
    datastuff2[i]=rnorm(1, mean=0, sd=0.05)
  }
  else{
    datastuff2[i]= rnorm(1, mean=2, sd=0.05)
  }
}
sd(datastuff2)
hist(datastuff2, 40, freq=F)
dists3=runif(NN)
datastuff3=vector(length=NN)
probs3=c(0.01, .95, 1)
for (i in 1:NN){
  if(dists3[i]<probs3[1]){
    datastuff3[i]=rnorm(1, mean=-2, sd=0.05)
  }else if (dists3[i]< probs3[2]){
    datastuff3[i]=rnorm(1, mean=0, sd=0.05)
  }
  else{
    datastuff3[i]= rnorm(1, mean=2, sd=0.05)
  }
}
sd(datastuff3)
hist(datastuff3, 40, freq=F)
xstart=-3
xend=3
NN=10000

#create density plots
curve(symf1(x), lwd=2,lty=1,from=xstart, to=xend, n=NN, col='black', ylab='Density') 


curve(symf2(x), lwd=2,lty=3,from=xstart, to=xend, n=NN, col='dodgerblue4', ylab='Density', add=TRUE) 
curve(newf(x), lwd=2, lty=2,from=xstart, to=xend, n=NN,col='firebrick4',  ylab='Density', add=TRUE)
#add legend
legend('topleft', legend=c("98% peak", "90% peak", "right bump"),
       col=c("black", "dodgerblue", "firebrick4"), lty=c(1,2,3), cex=1.)

#curve(symf1(x), lwd=2,lty=1,from=xstart, to=xend, n=NN, col='black', ylab='Density') 


#curve(symf2(x), lwd=2,lty=3,from=xstart, to=xend, n=NN, col='dodgerblue4', ylab='Density', add=TRUE) 
#curve(newf(x), lwd=2, lty=2,from=xstart, to=xend, n=NN,col='firebrick4',  ylab='Density', add=TRUE)
curve(dnorm(x, 0, sd=0.1), lwd=2,lty=1,from=xstart, to=xend, n=NN,col='black', add=T, ylab='Density')
curve(dnorm(x, 0, sd=0.5), lwd=2,lty=2,from=xstart, to=xend, n=NN,col='firebrick4', add=T, ylab='Density')
curve(dnorm(x, 0, sd=1), lwd=2,lty=3,from=xstart, to=xend, n=NN,col='darkgreen', add=T,ylab='Density')
q=0.25
s=0.34
s=(1-q)/q
sig=0.50
curve(sharkfin(x), lwd=2,lty=4,from=xstart, to=xend, n=NN,col='darkorchid', add=T, ylab='Density') 
beta=seq(-4, 4, by=0.001)
q=0.75
s=(1-q)/q
sig=1.25
curve(sharkfin(x), lwd=2,lty=5,from=xstart, to=xend, n=NN, col='dodgerblue4', ylab='Density', add=TRUE) 

#curve(dt(x, df=300), lwd=2,lty=1,from=-4, to=4, col='black', add=TRUE, ylim=c(0, 0.4), ylab='Density')

legend('topleft', legend=c("N(0,0.1)", "N(0,0.5)","N(0,1)", expression(paste("shark q=0.25, s=0.5, ", sigma, '=1.05')), 
                           expression(paste("shark q=0.75, s=1.25, ", sigma, '=0.88'))  ),
       col=c("black", "firebrick4", "darkgreen", "darkorchid", "dodgerblue4"), lty=c(1,2,3,4,5), cex=1.)

sharkfinq75


q=0.1
var=0.5^2  #the true variance we wanna set
a=0.3633802
b=0.7978846
varexp=(1-q)*a*((1-q)/q)^2 + q*a + (b*(1 + (1-q)/q))^2*q*(1-q)
#sig is really s
sig<-sqrt(var/varexp)
beta=seq(-4, 4, by=0.001)
print(sig)

curve(sharkfin(x), lwd=2.2,lty=1,from=xstart, to=xend, n=NN, col='black', ylab='Density') 
q=0.9
var=0.5^2  #the true variance we wanna set
a=0.3633802
b=0.7978846
varexp=(1-q)*a*((1-q)/q)^2 + q*a + (b*(1 + (1-q)/q))^2*q*(1-q)
#sig is really s
sig<-sqrt(var/varexp)
beta=seq(-4, 4, by=0.001)
#q=0.1
print(sig)
#sig=1.25
curve(sharkfin(x), lwd=2.2,lty=2,from=xstart, to=xend, n=NN, col='dodgerblue4', ylab='Density', add=T) 
#curve(dt(x, df=300), lwd=2,lty=1,from=-4, to=4, col='black', add=TRUE, ylim=c(0, 0.4), ylab='Density')

legend('topright', legend=c( expression(paste("shark q=0.1, s=0.09, ", sigma, '=0.25')), 
                             expression(paste("shark q=0.9, s=0.79, ", sigma, '=0.25'))  ),
       col=c("dodgerblue4", "black"), lty=c(2,1),cex=.75)
#setwd("C:/Users/demetri/Documents/ASU/RTG/fixed_plots/")
#ggsave('sharkfin25.pdf', sharkfinq25, width = 5, height = 4)
#ggsave('sharkfin75.pdf', sharkfinq75)
# 
# 
# plot(beta, sharkfin2(q,s,beta))
# #density(sharkfin2(.5,beta))
# 
# 
# 
# cauchymix=function(beta,v,m){
#   0.5*dt(beta, v, -m)+0.5*dt(beta, v, m)
# }
# 
# plot(beta,cauchymix(beta, 1,1.5))
# 
# cauchyplot=data.frame(u=beta, density=cauchymix(beta,1, 1.5))%>%
#   ggplot(aes(x=u, y=density))+geom_line(color='dodgerblue4',lwd=1.5)+ggtitle('Cauchy Mixture M=1.5, v=1')+
#   theme_minimal()+theme(plot.title = element_text(hjust = 0.5))
# 
# 
# #ggsave('cauchyplot.pdf', cauchyplot)
# 
# 
# #theme_set(theme_gray(base_size = 20))
theme_set(theme_minimal(base_size = 20))
i                  = 1
# 
vu                 = seq( -4.99,5, by=0.02)
# 
mf                 = matrix( NA, nrow = 3, ncol = 500 )
#library( 'truncnorm' )
# 
# 
# 
# 
# 
# f2=function(u){
#   .01*dnorm(u, mean=-2, sd=0.02)+0.94*dnorm(u,mean=0,sd=0.02)+.05*dnorm(u,mean=2,sd=0.02)
# }
# 
# mf[1,] = f2( vu )
# 
# 
# binded=cbind(vu, mf[1,])
# 
# binded<-as.data.frame(binded)
# #pdf( "bumpsNormal.pdf" )
# 
# 
# binded<-binded%>%
#   mutate(normed1=V2/sum(V2))
# binded
# 
# 
# rightbump<-binded%>% ggplot()+geom_density(aes(x=vu,y=normed1),color='dodgerblue4',stat='identity', lwd=.95)+xlim(-3,3)+labs(colour = "Cylinders")+
#   xlab('u')+ylab('f(u)')+
#   ggtitle(expression(paste('mean=-2,0,2, sigma=0.02'))) +
#   theme(plot.title = element_text(hjust = 0.5))+annotate('text', x = c(-2,0,2), y =c(.02,.35,.03), label=c("area=0.01 ",'area=0.94','area=0.05'), size=8)#+geom_vline(xintercept=c(c(-2.2,-1.8),c(-0.2,.2),c(1.8,2.2)))
# rightbump
# #setwd("C:/Users/demetri/Documents/ASU/RTG/newplots")
# 
# #ggsave('rightbump.pdf', rightbump)
# 
# 
# 
# 
# 
newf=function(u){
  .01*dnorm(u, mean=-2, sd=0.05)+0.94*dnorm(u,mean=0,sd=0.05)+.05*dnorm(u,mean=2,sd=0.05)
}
# 
mf[1,] = newf( vu )
# 
# 
binded2=cbind(vu, mf[1,])

binded2<-as.data.frame(binded2)
# #pdf( "bumpsNormal.pdf" )
# 
# 
binded2<-binded2%>%
  mutate(normed1=V2/sum(V2))
# 
rightbump2<-binded2%>% ggplot()+geom_density(aes(x=vu,y=normed1),color='#1d2951',stat='identity', lwd=.95)+xlim(-3,3)+labs(colour = "Cylinders")+
  xlab('u')+ylab('f(u)')+
  ggtitle(expression(paste('Asymmetric Mixture(',  sigma, '=0.49)'))) +
  theme(plot.title = element_text(hjust = 0.5,size=16))#+annotate('text', x = c(-2,0,2), y =c(.02,.18,.03), label=c("area=0.01 ",'area=0.94','area=0.05'), size=8)#+geom_vline(xintercept=c(c(-2.2,-1.8),c(-0.2,.2),c(1.8,2.2)))

rightbump2
ggsave("/home/dpapakos/sensitivity_analysis/rightbump2_navy.pdf", 
       rightbump2, height=6, width=6)
# #setwd("C:/Users/demetri/Documents/ASU/RTG/fixed_plots")
# #ggsave('rightbump2.pdf', rightbump2)
# rightbump2
# 
# ###################################################################
# # Plot our densities for visual aid
# ###################################################################
# i                  = 1
# #vsd                = c(0.05, 0.2, 1, 2)
# area=c(.01,.1,.5)
# vu                 = seq( -4.99,5, by=0.02)
# 
# mf                 = matrix( NA, nrow = length(area), ncol = 500 )
# #library( 'truncnorm' )
# 
# for ( j in  area ){
#   
#   #	f = function(u){	
#   #		dtruncnorm(u, a = 0, b = Inf, mean = 1, sd = sd)
#   #	}
#   #f=function(u){
#   # .2*dnorm(u, mean=-2, sd=0.1)+0.7*dnorm(u,mean=0,sd=0.1)+.1*dnorm(u,mean=2,sd=0.1)
#   #}
#   
#   f=function(u){
#     .05*dnorm(u, mean=-2, sd=0.1)+0.05*dnorm(u, mean=-.75, sd=0.1)+0.75*dnorm(u,mean=0,sd=0.1)+0.05*dnorm(u,mean=.75, sd=0.1)+.1*dnorm(u,mean=2,sd=0.1)
#   }
#   
#   mf[i,] = f( vu )
#   i = i + 1
# }
# stu=cbind(vu, mf[1,])
# stu=cbind(stu, mf[2,])
# stu=cbind(stu, mf[3,])
# stu<-as.data.frame(stu)
# #pdf( "bumpsNormal.pdf" )
# 
# 
# stu<-stu%>%
#   mutate(normed1=V2/sum(V2),
#          normed2=V3/sum(V3),
#          normed3=V4/sum(V4))
# 
# #stu<-stu%>%
# # select(vu, normed1, normed2, normed3)
# #all<-stu%>%
# # gather(vu,c(normed1,normed2,normed3))
# #colnames(all) <- c("run","u")
# all<-cbind(all, vu)
# 
# stu%>% ggplot()+geom_density(aes(x=vu,y=normed1,color='sigma=0.01'),stat='identity', lwd=1.1)+xlim(-3,3)+geom_density(aes(x=vu, y=normed2, color='sigma=0.1'), stat='identity', lwd=1)+geom_density(aes(x=vu, y=normed3, color='sigma=0.5'), stat='identity',lwd=1)+labs(colour = "Cylinders")+
#  scale_color_manual( values =c("firebrick4",'dodgerblue2', 'black'),
#                 name="sigma")+
# xlab('u')+ylab('f(u)')+
#  ggtitle(expression(paste('mean=0'))) +
# theme(plot.title = element_text(hjust = 0.5))
# 
# 
# weird2plot=stu%>% ggplot()+geom_density(aes(x=vu,y=normed1),stat='identity', fill='dodgerblue4',lwd=1, alpha=0.29)+xlim(-3,3)+
#   xlab('u')+ylab('f(u)')+
#   ggtitle(expression(paste('multimodal 2'))) +theme_minimal()+
#   theme(plot.title = element_text(hjust = 0.5))+annotate("text", x = -2, y = .016, label = '5 %', size=6)+
#   annotate("text", x = 0, y = .062, label = '75 %', size=6)+
#   annotate("text", x = 2, y = .016, label = '10 %', size=6)+
#   annotate("text", x = -.75, y = .016, label = '5 %', size=6)+
#   annotate("text", x = .75, y = .016, label = '5 %', size=6)
# 
# #ggsave('weirdplot.pdf', weird2plot)
# weird2plot
# 
# 
# #define marginal density of u; must be valid over (-Inf, Inf)
# 
# ###################################################################
# # Plot truncated normal density. 
# ###################################################################
# 
# i                  = 1
# #vsd                = c(0.05, 0.2, 1, 2)
# area=c(.02,.1,.2)
# vu                 = seq( -4.99,5, by=0.02)
# 
# mf                 = matrix( NA, nrow = length(area), ncol = 500 )
# #library( 'truncnorm' )
# 
# 
# 
# 
# for ( j in  area ){
#   
#   
#   f=function(u){
#     (j/2)*dnorm(u, mean=-2, sd=0.05)+(1-j)*dnorm(u,mean=0,sd=0.05)+(j/2)*dnorm(u,mean=2,sd=0.05)
#   }
#   
#   mf[i,] = f( vu )
#   i = i + 1
# }
# 
# stu=cbind(vu, mf[1,])
# stu=cbind(stu, mf[2,])
# stu=cbind(stu, mf[3,])
# stu<-as.data.frame(stu)
# #pdf( "bumpsNormal.pdf" )
# 
# 
# stu<-stu%>%
#   mutate(normed1=V2/sum(V2),
#          normed2=V3/sum(V3),
#          normed3=V4/sum(V4))
# 
# #stu<-stu%>%
# # select(vu, normed1, normed2, normed3)
# #all<-stu%>%
# # gather(vu,c(normed1,normed2,normed3))
# #colnames(all) <- c("run","u")
# all<-cbind(all, vu)
# 
# multimod1<-stu%>% ggplot()+geom_density(aes(x=vu,y=normed1,color='area=0.01'),stat='identity', lwd=.8)+xlim(-3,3)+geom_density(aes(x=vu, y=normed2, color='area=0.05'), stat='identity', lwd=.8)+geom_density(aes(x=vu, y=normed3, color='area=0.1'), stat='identity',lwd=.8)+labs(colour = "Cylinders")+
#   scale_color_manual( values =c("firebrick4",'dodgerblue2', 'black'),
#                       name="area of little bump")+
#   xlab('u')+ylab('f(u)')+
#   ggtitle(expression(paste(mu,'=-2,0,2,',  sigma, '=0.05'))) +
#   theme(plot.title = element_text(hjust = 0.5))
# multimod1
# #ggsave(filename="multimod1.pdf", plot=multimod1)
# #ggtitle("mu=0.2,1,2,sd=.05") +
# #theme(plot.title = element_text(hjust = 0.5))#+geom_vline(xintercept=c(x_0,x_1,x_2))+annotate('text', x = c(.2,1,2), y =c(.01,.01,.01), label=c(".05 ",'.75','.20'))
# 
# #plot( range( vu ), range( c(mf) ), type = "n", main = "Different normals" )
# 
# #for (i in 1:length(area) ){
# #	lines( vu, mf[i,], lty = i, col = i, lwd = 2 )
# #}
# #legend( "topright", legend = paste( "sd=0.1, area in each bump =", area/2), lty = 1:length(area), col = 1:length(area) )
# #dev.off()