## sigma = mu ; y~N(mu,mu^2)
mu = seq(1,10,by=.01)
dep1.mu = pnorm(mean=3,sd=3,abs(mu)*(1+sqrt(5)/2))-pnorm(mean=3,sd=3,abs(mu)*(1-sqrt(5)/2))
dep2.mu = 1-dep1.mu
dep.mu = ifelse(dep1.mu>dep2.mu,dep2.mu,dep1.mu)
plot(mu,dep.mu,type='l')

d1.mu = pnorm(mean=3,sd=3,mu)
d.mu = ifelse(d1.mu<.5,d1.mu,1-d1.mu)

plot(mu,d.mu,type='l',lwd=2,col='red',xlim=c(1,10))
lines(mu,dep.mu,col='blue',lwd=2)

mu[which.max(dep.mu)]