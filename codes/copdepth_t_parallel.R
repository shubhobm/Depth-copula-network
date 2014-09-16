## t-copula depth
library(VineCopula)
###### required for parallel computing
library(parallel)
library(doSNOW)

set.seed(007)
N = 100; df = 10
dep.summary = matrix(rep(0,27), ncol=3)
rho.range = seq(.1,.9,.1)

## Simulates max depth copula estimators for rho  = 0.1,...,0.9

cl <- makeCluster(detectCores())
registerDoSNOW(cl)

for(rho in rho.range){
  maxdeps = rep(0,1000)
  # Generate 1000 samples of size 100 each then caculates depth
  # uses parallel computing  
  loopfun = function(rho){
    # Generate sample
    require(VineCopula)
    U = BiCopSim(N, 2, rho, par2 = df)
    
    u=U[,1]; v = U[,2]
    rhovals = seq(0,.99,.01); j=1
    deps = rep(0,length(rhovals))
    for(rho1 in rhovals){
      del.loglik = BiCopDeriv(u,v,family=2,par=rho1, par2=df, deriv="par") / 
        BiCopPDF(u,v,family=2,par=rho1, par2=df)
      deps[j] = min(sum(del.loglik>0), sum(del.loglik<0)) + sum(del.loglik==0)
      j = j+1
    }
    return(rhovals[which.max(deps)])
  }  
  system.time(maxdeps <- foreach(i=1:1000) %dopar% loopfun(rho))
  maxdeps = matrix(unlist(maxdeps), ncol=1)
  
  dep.summary[rho*10,] = c(rho, mean(maxdeps), sd(maxdeps))
}
stopCluster(cl)

plot(dep.summary[,2]~rho.range, xlim=c(0,1),ylim=c(0,1), type="p", pch=19, cex=1.2,
     main=paste("Simulated maximum t-copula (df=",df,") depth estimators"),
     xlab="Real parameter", ylab="Parameter with maximum depth")
for(i in 1:9){
  segments(dep.summary[i,1],dep.summary[i,2]-dep.summary[i,3],
           dep.summary[i,1],dep.summary[i,2]+dep.summary[i,3], lwd=2)
}
abline(0,1, lty=2)
dep.summary