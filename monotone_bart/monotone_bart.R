
library(dbarts)
library(bcf)
library(fastbart)

set.seed(18022)

monotone_bart = function(y, z, x, xpred, nskip=5000, ndpost=5000, m = 50) {
  
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

# 
# 
# n = 500
# x1 = runif(n)
# x2 = runif(n)
# x = cbind(x1, x2)
# f1 = function(x) -0.5 + x[1]^2
# f0 = function(x) -.5 + 0.75*x[1]^2 - 0.1*x[2]
# z = rbinom(n, 1, 0.5)
# py = pnorm(apply(x, 1, f1))*z + pnorm(apply(x, 1, f0))*(1-z)
# y = rbinom(n, 1, py)
# xpred = as.matrix(expand.grid(x1=(0:100)/100, x2 = c(0, 0.5)))
# 
# res = monotone_bart(y, z, x, xpred, 4000, 4000, 50)
# 
# plot(xpred[1:101,1], colMeans(res$pr1[,1:101]), type='l', ylim=c(0,1))
# lines(xpred[1:101,1], colMeans(res$pr0[,1:101]), type='l', col='red', ylim=c(0,1))
# 
# plot(xpred[1:101,1], colMeans(res$pr1[,102:202]), type='l', ylim=c(0,1))
# lines(xpred[1:101,1], colMeans(res$pr0[,102:202]), type='l', col='red', ylim=c(0,1))
