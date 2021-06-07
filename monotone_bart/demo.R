
library(dbarts)
library(bcf)
library(fastbart)

set.seed(18022)

setwd('~/Dropbox/Observation model/code/binaryiv')

n = 500
x1 = runif(n)
x2 = runif(n)
x = cbind(x1, x2)
f1 = function(x) -0.5 + x[1]^2
f0 = function(x) -.5 + 0.75*x[1]^2 - 0.1*x[2]
z = rbinom(n, 1, 0.5)
py = pnorm(apply(x, 1, f1))*z + pnorm(apply(x, 1, f0))*(1-z)
y = rbinom(n, 1, py)

sort_ix = order(z, y)
x = x[sort_ix,]
z = z[sort_ix]
y = y[sort_ix]

#subs = !(x[,1] < 0.4 & x[,1] > 0.2)
#x = x[subs,]
#z = z[subs]
#y = y[subs]

n0 = sum(z==0)
n00 = sum(z==0 & y==0)
yobs = y
yobs0 = rbinom(n, 1, 0.5)
yobs[1:n00] = rbinom(n00, 1, 0.5) # These are DA variables
yobs0[1:n00] = rbinom(n00, 1, 0.5) # These are DA variables
yobs0[1:n00][yobs[1:n00]==1] = 0 # To statisfy DA constraints
yobs0[z==0 & y==1] = 1
yobs0 = yobs0[1:n0] # I guess I made this too long? ¯\_(ツ)_/¯

offset =  0#qnorm(mean(y[z==1]))
offset0 = 0# qnorm(mean(flu$grp==0 & flu$fluy2==1)/mean(flu$grp==1 & flu$fluy2==1)) #<- wtf
zz  = rtnorm(length(yobs), mean=0, lower = ifelse(yobs,-offset, -Inf), upper = ifelse(yobs, Inf, -offset))
z0 = rtnorm(length(yobs0), mean=0, lower = ifelse(yobs0,-offset0, -Inf), upper = ifelse(yobs0, Inf, -offset0))

zz = offset + 3*yobs - 3*(1-yobs)
z0 = offset0 + 3*yobs0 - 3*(1-yobs0)

################################################################################
# MCMC
################################################################################

nskip=5000
ndpost=5000

xpred = as.matrix(expand.grid(x1=(0:100)/100, x2 = c(0, 0.5)))

set.seed(1022)

xi = lapply(1:ncol(x), function(i) bcf:::.cp_quantile(x[,i]))
fit.mono = bartRcppMono(yobs, zz, t(as.matrix(x)), t(xpred),
                        yobs0, z0, t(as.matrix(x)),t(xpred),
                        n00,
                        xi,
                        nskip, ndpost, 50, 3.0,
                        offset, offset0)

xpred.exp = rbind(data.frame(xpred, z=1), data.frame(xpred, z=0))
fit =  bart(cbind(x, z), y, xpred.exp,
               nskip=5000, ndpost=5000,
               ntree=50, usequants=T)


# P(Y|Z=1, X)
pr  = pnorm(fit.mono$postpred)

# P(Y|Z=0, X)
pr0 = pr*pnorm(fit.mono$postpred0)

# Z=1, X2=0
plot(xpred[1:101,1], colMeans(pr[,1:101]), type='l', ylim=c(0,1))
lines(xpred[1:101,1], pnorm(apply(xpred[1:101,], 1, f1)), type='l', ylim=c(0,1))
lines(xpred[1:101,1],pnorm(colMeans(fit$yhat.test[,1:101])), lty=2)

# Z=0, X2=0
lines(xpred[1:101,1], colMeans(pr0[,1:101]), col='red')
lines(xpred[1:101,1], pnorm(apply(xpred[1:101,], 1, f0)), type='l', ylim=c(0,1), col='red')
lines(xpred[1:101,1],pnorm(colMeans(fit$yhat.test[,203:303])), lty=2, col='red')

#plot(xpred[1:101,1], colMeans(pnorm(fit.mono$postpred0)[,1:101]), type='l', ylim=c(0,1))
#lines(xpred[102:202,1], colMeans(pnorm(fit.mono$postpred0)[,102:202]), type='l', ylim=c(0,1), col='red')

# Z=1, X2=0.5
plot(xpred[102:202,1], colMeans(pr[,102:202]), type='l', ylim=c(0,1))
lines(xpred[102:202,1], pnorm(apply(xpred[102:202,], 1, f1)), type='l', ylim=c(0,1))
lines(xpred[102:202,1],pnorm(colMeans(fit$yhat.test[,102:202])), lty=2)
# Z=0, X2=0.5
lines(xpred[102:202,1], colMeans(pr0[,102:202]), col='red')
lines(xpred[102:202,1], pnorm(apply(xpred[102:202,], 1, f0)), type='l', ylim=c(0,1), col='red')
lines(xpred[102:202,1], pnorm(colMeans(fit$yhat.test[,304:404])), lty=2, col='red')
