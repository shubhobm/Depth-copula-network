g = function(x,mu){ return (cbind(x-mu, x^3-mu^3-3*mu))}

gprime = function(mu){ return (c(-1,-3*mu^2-3))}

logstar = function(z,del=1e-5){
  if(z>del)
    return(log(z))
  else
    return(log(del)-1.5+2*z/del-z^2/(2*del^2))
}

fun2 <- function(z,x,mu){return(-sum(apply(as.matrix(1+g(x,mu) %*% z),1,logstar)))}

## Generate sample
n=100
x = rnorm(n, mean=5, sd=1)

lt = length(tval); val = matrix(nrow=lt, ncol=1)
# mu=9
# for(i in 1:lt){
#   fval = fun1(t[i],x,mu)
#   val[i] = fval
# }
# plot(val~t, type='l')


require(numDeriv)
muvals = seq(4.9,5.1,.01)
de.mu = rep(0,length(muvals))
i=1

for(mu in muvals){
  ## Solve fun2 for t using Newton's method
  maxit = 1e5; tol = 1e-10; t0 = c(.2,.3)
  for(iter in 1:maxit){
    gp = gprime(mu)
    gvals = colSums(g(x,mu))
    denom = 1+sum(t0*gvals)
    grad = -gp/denom
    Hess = crossprod(t(gp))/denom^2
    d0 = -grad
    t1 = t0 + d0
    err = (t1-t0)/t0
    # cat(iter, "\t", t0, "->", t1, "\t", grad, "\t", Hess, "\t", fun2(t1,x,mu),"\n")
    if(abs(err)<tol)
      break()
    else
      t0 = t1
  }
  t=t1
  
  ## then calculate depth
  lprime = t/(1+t*(x-mu)); nm<-sum(lprime>0)
  de.mu[i] = min(nm,n-nm)/n; i=i+1
  cat(mu,"\t",t,"\t",grad,"\t",nm,"\t",n-nm,"\n")
}

plot(de.mu~muvals,type='l',lwd=2)