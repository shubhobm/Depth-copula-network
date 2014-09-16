## Code to determine copula that fits the data best
MD.bias = function(est, family){
  if(family==1){
    if(est>0.461) corr.est = -1.22217*est^3 + 3.6434*est^2 - 1.42154*est + 0.000396
    else corr.est = 0
  }
  else if(family==3) corr.est = -0.50675 + 0.713898*est
  else corr.est = 0.015 + 0.71*est
  
  return(corr.est)
}

# generate data
library(VineCopula)
library(copula)
N = 1000; df=5
U = BiCopSim(N, family=1, par=.4)
u = U[,1]; v = U[,2]

# Generate empirical copula for the (pseudo-)data
empC = C.n(U, U=U)

# Given data, calculate ML parameters for different copula.. Gaussian, Clayton, Gumbel, BB1
MLparams = matrix(rep(NA,8),ncol=2)
for(i in 1:4){
  est = BiCopEst(u, v, family=i, method='mle')
  MLparams[i,1] = est$par
  if(i==2){
    fl = floor(est$par2)
    dfval = ifelse(est$par2-fl<0.5,fl,fl+1)
    MLparams[i,2] = ifelse(dfval<3,3,dfval)
  }
}
MLparams

# obtain distance of all copulas from empirical copula
D2.ml = rep(0,4)
for(i in 1:4){
  if(i==2)
    D2.ml[i] = sum((pCopula(U, tCopula(MLparams[i,1],df=MLparams[i,2])) - empC)^2)
  else
    D2.ml[i] = sum((BiCopCDF(u,v,family=i,par=MLparams[i,1]) - empC)^2)
}
D2.ml

# Given data, calculate Max depth parameters for different copula
MDparams = matrix(rep(NA,8),ncol=2)
for(i in 1:4){
  if(i==1||i==2) parvals = seq(0,.99,.01)
  else parvals = seq(1.1,20,.01)

  deps = rep(0,length(parvals)); j=1
  for(par1 in parvals){
    if(i==2){
      del.loglik = BiCopDeriv(u,v,family=i,par=par1, par2=MLparams[i,2], deriv="par") / 
        BiCopPDF(u,v,family=i,par=par1, par2=MLparams[i,2])
    }
    else{
      del.loglik = BiCopDeriv(u,v,family=i,par=par1, deriv="par") / 
        BiCopPDF(u,v,family=i,par=par1)
    }      
    
    deps[j] = min(sum(del.loglik>0), sum(del.loglik<0)) + sum(del.loglik==0)
    j = j+1
  }
  
  MDparams[i,1] = MD.bias(parvals[which.max(deps)], family=i)

}
MDparams[2,2] = MLparams[2,2]
MDparams

# obtain distance of all copulas from empirical copula
D2.md = rep(0,4)
for(i in 1:4){
  if(i==2)
    D2.md[i] = sum((pCopula(U, tCopula(MDparams[i,1],df=MDparams[i,2])) - empC)^2)
  else
    D2.md[i] = sum((BiCopCDF(u,v,family=i,par=MDparams[i,1]) - empC)^2)
}
D2.md
 