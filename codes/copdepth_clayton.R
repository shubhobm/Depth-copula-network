## Clayton copula depth
library(VineCopula)
set.seed(006)

N = 100; maxit=50

theta.range = c(1.5,2,5,8,10)
dep.summary = matrix(rep(0,3*length(theta.range)), ncol=3)
theta.cnt=1

## Simulates max depth copula estimators for all values in theta.range
for(theta in theta.range){
  maxdeps = rep(0,1000)
  # Generate 1000 samples of size 100 each
  for(i in 1:1000){
    # Generate sample from Clayton copula
    X = rgamma(N, shape=1/theta, rate=1)
    V = matrix(runif(2*N), ncol=2)
    U = (1-log(V)/X)^(-1/theta)
    
    u = U[,1]; v = U[,2]
    
    est = BiCopEst(u, v, family=3, method='mle')
    interval = c(1,10*est$par)
    
    iterating=TRUE; iter=1
    while(iterating){      
      theta.list = seq(interval[1], interval[2], length.out=10)    
      dep.list = rep(0,10)
      for(j in 1:10){
        del.loglik = BiCopDeriv(u,v,family=3,par=theta.list[j], deriv="par") / 
          BiCopPDF(u,v,family=3,par=theta.list[j])
        
        dep.list[j] =(min(sum(del.loglik>0), sum(del.loglik<0)) + sum(del.loglik==0))/N
      }
      
      # stopping criteria
      maxdep = max(na.omit(dep.list))
      if(maxdep==0.5 || iter==maxit) iterating = FALSE
      # otherwise determine new interval
      else{
        max.index = which(dep.list==maxdep)
        if(length(max.index)==1){
          lm = ifelse(max.index==1, interval[1], theta.list[max.index-1])
          rm = ifelse(max.index==10, interval[2], theta.list[max.index+1])
        }
        else{
          m1 = min(max.index); m2 = max(max.index)
          lm = ifelse(m1==1, interval[1], theta.list[m1-1])
          rm = ifelse(m2==10, interval[2], theta.list[m2+1])
        }
        interval = c(lm,rm)
        iter=iter+1
      }
    }

    maxdeps[i] = theta.list[which.max(dep.list)]
  }
  
  dep.summary[theta.cnt,] = c(theta, mean(maxdeps), sd(maxdeps))
  theta.cnt = theta.cnt+1
}

plot(dep.summary[,2]~theta.range, xlim=c(1,10),ylim=c(1,20), type="p", pch=19, cex=1.2,
     main="Clayton copula",
     xlab="Real parameter", ylab="Parameter with maximum depth")
for(i in 1:5){
  segments(dep.summary[i,1],dep.summary[i,2]-dep.summary[i,3],
           dep.summary[i,1],dep.summary[i,2]+dep.summary[i,3], lwd=2)
}
abline(0,1, lty=2)
dep.summary