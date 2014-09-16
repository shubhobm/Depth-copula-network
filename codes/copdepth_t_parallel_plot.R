## t copula depth
library(VineCopula)
###### required for parallel computing
library(parallel)
library(doSNOW)

set.seed(010)
N = 100; df = 5

rho.range = seq(0.1,0.9,0.1)
dep.summary = matrix(rep(0,3*length(rho.range)), ncol=3)


## Simulates max depth copula estimators for all values in theta.range

loopfun = function(rho){
  library(VineCopula)
  maxdeps = rep(0,1000); rho.cnt=1
  # Generate 1000 samples of size 100 each
  for(i in 1:1000){
    # Generate sample from t copula
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

    maxdeps[i] = rhovals[which.max(deps)]
  }
  
 return(list(rho, mean(maxdeps), sd(maxdeps)))
}

cl <- makeCluster(detectCores())
registerDoSNOW(cl)

system.time(dep.summary <- foreach(rho = rho.range) %dopar% loopfun(rho))
(dep.summary1 = matrix(unlist(dep.summary), ncol=3, byrow=T))
stopCluster(cl)

plot(dep.summary1[,2]~rho.range, xlim=c(0,1),ylim=c(0,1), type="p", pch=19, cex=1.2,
     main=paste("t-copula (df=",df,")"),
     xlab="Real parameter", ylab="Parameter with maximum depth")
for(i in 1:9){
  segments(dep.summary1[i,1],dep.summary1[i,2]-dep.summary1[i,3],
           dep.summary1[i,1],dep.summary1[i,2]+dep.summary1[i,3], lwd=2)
}
abline(0,1, lty=2)

# Fit cubic regression model
th = dep.summary1[,1]; bth = dep.summary1[,2]
summary(m <- lm(th ~ bth+I(bth^2)+I(bth^3)))
m$coef
head(dep.summary1)

 save.image(filename="fuck.rda")