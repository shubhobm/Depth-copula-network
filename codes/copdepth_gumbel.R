## Gumbel copula depth
library(gumbel)
set.seed(005)

N = 100
dep.summary = matrix(rep(0,15), ncol=3)
theta.range = c(1.5,2,5,8,10); theta.cnt=1

## Simulates max depth copula estimators for rho  = 0.1,...,0.9
for(theta in theta.range){
  maxdeps = rep(0,1000)
  # Generate 1000 samples of size 100 each
  for(i in 1:1000){
    # Generate sample from Gumbel copula
    U = rgumbel(N,theta)
    
    u = U[,1]; v = U[,2]
    thetavals = seq(1.1,20,.1); j=1
    deps = rep(0,length(thetavals))
    for(theta1 in thetavals){
      del.loglik = BiCopDeriv(u,v,family=4,par=theta1,deriv="par") / BiCopPDF(u,v,family=4,par=theta1)
      deps[j] = min(sum(del.loglik>0), sum(del.loglik<0)) + sum(del.loglik==0)
      j = j+1
    }
    
    maxdeps[i] = thetavals[which.max(deps)]
  }
  
  dep.summary[theta.cnt,] = c(theta, mean(maxdeps), sd(maxdeps))
  theta.cnt = theta.cnt+1
}

plot(dep.summary[,2]~theta.range, xlim=c(1,10),ylim=c(1,20), type="p", pch=19, cex=1.2,
     main="Simulated maximum Gumbel copula depth estimators",
     xlab="Real parameter", ylab="Parameter with maximum depth")
for(i in 1:5){
  segments(dep.summary[i,1],dep.summary[i,2]-dep.summary[i,3],
           dep.summary[i,1],dep.summary[i,2]+dep.summary[i,3], lwd=2)
}
abline(0,1, lty=2)
dep.summary