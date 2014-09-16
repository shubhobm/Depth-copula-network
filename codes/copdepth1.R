## Normal copula depth
require(MASS)
set.seed(004)

N = 100
dep.summary = matrix(rep(0,27), ncol=3)
rho.range = seq(.1,.9,.1)

## Simulates max depth copula estimators for rho  = 0.1,...,0.9
for(rho in rho.range){
  maxdeps = rep(0,1000)
  # Generate 1000 samples of size 100 each
  for(i in 1:1000){
    # Generate sample
    X = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1),nrow=2))
    U = apply(X,2,pnorm)
    
    x=X[,1]; y = X[,2]
    rhovals = seq(0,.99,.01)
    deps = rep(0,length(rhovals))
    for(rho1 in rhovals){
      # find out likelihood depth at rho1
      del.loglik = (-rho1*y^2+(1+rho1^2)*x*y+rho1-rho1^3-rho1*x^2)*(1-rho1^2)^(-2)/(dnorm(x)*dnorm(y))
      j = rho1*100+1
      deps[j] = min(sum(del.loglik>0), sum(del.loglik<0)) + sum(del.loglik==0)
    }
    maxdeps[i] = rhovals[which.max(deps)]
  }
  
  dep.summary[rho*10,] = c(rho, mean(maxdeps), sd(maxdeps))
}

plot(dep.summary[,2]~rho.range, xlim=c(0,1),ylim=c(0,1), type="p", pch=19, cex=1.2,
     main="Simulated maximum Gaussian copula depth estimators",
     xlab="Real parameter", ylab="Parameter with maximum depth")
for(i in 1:9){
  segments(dep.summary[i,1],dep.summary[i,2]-dep.summary[i,3],
           dep.summary[i,1],dep.summary[i,2]+dep.summary[i,3], lwd=2)
}
abline(0,1, lty=2)
