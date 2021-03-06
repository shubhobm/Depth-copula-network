## Clayton copula depth
library(VineCopula)
###### required for parallel computing
library(parallel)
library(doSNOW)

set.seed(006); N = 100
theta.range=seq(0.1,.9,.1)
dfs=c(7:9); mat=matrix(rep(0,6*length(dfs)), ncol=6)

for(i in 1:length(dfs)){
  
  df = dfs[i]
  
  ## Simulates max depth copula estimators for all values in theta.range
  loopfun = function(theta){
    library(VineCopula)
    maxdeps = rep(0,1000)
    # Generate 1000 samples of size 100 each
    for(i in 1:1000){
      # Generate sample from t copula
      U = BiCopSim(N, family=2, par=theta, par2=df)
      u = U[,1]; v = U[,2]
      
      est = BiCopEst(u, v, family=2, method='mle')
      interval = c(0,.99)
      
      # interval search algorithm
      iterating=TRUE; iter=1; maxit=10
      while(iterating){      
        theta.list = seq(interval[1], interval[2], length.out=10)    
        dep.list = rep(0,10)
        for(j in 1:10){
          del.loglik = BiCopDeriv(u,v,family=2,par=theta.list[j], par2=df, deriv="par") / 
            BiCopPDF(u,v,family=2,par=theta.list[j], par2=df)
          
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
    
    return(list(theta, mean(maxdeps), sd(maxdeps)))
  }
  
  dep.summary = matrix(rep(0,3*length(theta.range)), ncol=3)
  
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  system.time(dep.summary <- foreach(theta = theta.range) %dopar% loopfun(theta))
  dep.summary = matrix(unlist(dep.summary), ncol=3, byrow=T)
  stopCluster(cl) 
  
  plot(dep.summary[,2]~theta.range, xlim=c(0,1),ylim=c(0,1), type="p", pch=19, cex=1.2,
       main=paste("t-copula (df=",df,")"),
       xlab="Real parameter", ylab="Parameter with maxim um depth")
  for(i in 1:9){
    segments(dep.summary[i,1],dep.summary[i,2]-dep.summary[i,3],
             dep.summary[i,1],dep.summary[i,2]+dep.summary[i,3], lwd=2)
  } 
  abline(0,1, lty=2)
  
  # Fit cubic regression model
  th = dep.summary[,1]; bth = dep.summary[,2]
  summary(m <- lm(th ~ bth+I(bth^2)+I(bth^3)))
  
  mat[i,] = c(df, m$coef, dep.summary[1,2])
}

mat