## function to generate randomly from stable distribition
rstable = function (n, alpha, beta, gamma, delta){
  x.list = rep(0,n)
  
  for(i in 1:n){
    Th = runif(1)*pi - pi/2
    W = rexp(1, rate=1)
    theta0 = atan(beta*tan(pi*alpha/2))/alpha
    pi.by2.plus.betaTh = pi/2+beta*Th
    
    if(alpha==1){
      Z = 2/pi*(pi.by2.plus.betaTh*tan(Th) - 
                  beta*log(pi/2*W*cos(Th)/pi.by2.plus.betaTh))
      x.list[i] = gamma*Z + delta + 2*beta/pi*gamma*log(gamma)
    }
    else{
      Z = sin(alpha*(theta0+Th))/(cos(alpha*theta0)*cos(Th))^(1/alpha)*
        (cos(alpha*theta0+(alpha-1)*theta)/W)^(1/alpha-1)
      x.list[i] = gamma*Z + delta
    }
  }
  
  return(x.list)
}