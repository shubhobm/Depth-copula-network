## Clayton copula depth
library(VineCopula)
###### required for parallel computing
library(parallel)
library(doSNOW)

set.seed(006)
N = 100

theta.range = seq(.1,10,0.1)
dep.summary = matrix(rep(0,3*length(theta.range)), ncol=3)


## Simulates max depth copula estimators for all values in theta.range

cl <- makeCluster(detectCores())
registerDoSNOW(cl)

loopfun = function(theta){
  library(VineCopula)
  maxdeps = rep(0,1000); theta.cnt=1
  # Generate 1000 samples of size 100 each
  for(i in 1:1000){
    # Generate sample from Clayton copula
    X = rgamma(N, shape=1/theta, rate=1)
    V = matrix(runif(2*N), ncol=2)
    U = (1-log(V)/X)^(-1/theta)
    
    u = U[,1]; v = U[,2]
    thetavals = seq(1.1,20,.1); j=1
    deps = rep(0,length(thetavals))
    for(theta1 in thetavals){
      del.loglik = BiCopDeriv(u,v,family=3,par=theta1,deriv="par") / BiCopPDF(u,v,family=3,par=theta1)
      deps[j] = min(sum(del.loglik>0), sum(del.loglik<0)) + sum(del.loglik==0)
      j = j+1
    }

    maxdeps[i] = thetavals[which.max(deps)]
  }
  
 return(list(theta, mean(maxdeps), sd(maxdeps)))
}

system.time(dep.summary <- foreach(theta = theta.range) %dopar% loopfun(theta))
dep.summary = matrix(unlist(dep.summary), ncol=3, byrow=T)
stopCluster(cl)

plot(dep.summary[,2]~theta.range, xlim=c(1,10),ylim=c(1,20), type="p", pch=19, cex=1.2,
     main="Simulated maximum Clayton copula depth estimators",
     xlab="Real parameter", ylab="Parameter with maximum depth")
for(i in 1:5){
  segments(dep.summary[i,1],dep.summary[i,2]-dep.summary[i,3],
           dep.summary[i,1],dep.summary[i,2]+dep.summary[i,3], lwd=2)
}
abline(0,1, lty=2)
head(dep.summary)

# Fit linear regression model
th = dep.summary[,1]; bth = dep.summary[,2]
summary(m <- lm(th ~ bth))
m$coef

