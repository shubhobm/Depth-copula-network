fun1 <- function(z,x,mu){return(sum((x-mu)/(1+z*(x-mu))))}

logstar = function(z,del=1e-5){
  if(z>del)
    return(log(z))
  else
    return(log(del)-1.5+2*z/del-z^2/(2*del^2))
}

fun2 <- function(z,x,mu){return(-sum(apply(as.matrix(1+z*(x-mu)),1,logstar)))}

## Generate sample
n=1000
x = rnorm(n, mean=5, sd=2)

x = rt(n,df=100,ncp=5)

x = rcauchy(n,location=5, scale=2)

x = rchisq(n,df=10,ncp=5)
lb = max((1/n-1)/x)
tval = seq(lb,20,.001)

lt = length(tval); val = matrix(nrow=lt, ncol=1)
# mu=9
# for(i in 1:lt){
#   fval = fun1(t[i],x,mu)
#   val[i] = fval
# }
# plot(val~t, type='l')

muvals = seq(2.9,7.1,.01)
de.mu = rep(0,length(muvals))
i=1

for(mu in muvals){
  ## Solve fun2 for t using Newton's method
  maxit = 1e5; tol = 1e-10; t0 = 2
  for(iter in 1:maxit){
    grad = -fun1(t0,x,mu)
    Hess = sum((x-mu)^2/(1+t0*(x-mu))^2)
    d0 = -grad/Hess
    t1 = t0 + d0
    err = (t1-t0)
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