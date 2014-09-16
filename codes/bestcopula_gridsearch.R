setwd('C:/study/my projects/depth')
source('functions.R')

# generate data
library(VineCopula)
library(copula)
N = 1000; df=13
U = BiCopSim(N, family=1, par=.4)
u = U[,1]; v = U[,2]

# Generate empirical copula for the (pseudo-)data
empC = C.n(U, U=U)

# Given data, calculate ML parameters for different copula.. Gaussian, Clayton, Gumbel, BB1
MLparams = matrix(rep(NA,8),ncol=2); MDparams = MLparams
for(i in 1:4){
  eml = CopEst(u,v, family=i, method='mle')
  emd = CopEst(u,v, family=i, method='mde')
  MLparams[i,] = c(eml$par, eml$par2)
  MDparams[i,] = c(emd$par, emd$par2)
}
MLparams
MDparams

# obtain distance of all copulas from empirical copula. ML estimators
D2.ml = rep(0,4)
for(i in 1:4){
  if(i==2)
    D2.ml[i] = sum((pCopula(U, tCopula(MLparams[i,1],df=MLparams[i,2])) - empC)^2)
  else
    D2.ml[i] = sum((BiCopCDF(u,v,family=i,par=MLparams[i,1]) - empC)^2)
}
D2.ml

# obtain distance of all copulas from empirical copula... MLD estimators
D2.md = rep(0,4)
for(i in 1:4){
  if(i==2)
    D2.md[i] = sum((pCopula(U, tCopula(MDparams[i,1],df=MDparams[i,2])) - empC)^2)
  else
    D2.md[i] = sum((BiCopCDF(u,v,family=i,par=MDparams[i,1]) - empC)^2)
}
D2.md
 