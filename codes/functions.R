setwd('C:/study/my projects/depth')
load('t.coefs.rda')

## Code to determine copula that fits the data best
MD.bias = function(est, family, tdf=NULL){
  aest = abs(est)
  if(family==1){
     if(aest>0.461) corr.est = sign(est)*(-1.22217*aest^3 + 3.6434*aest^2 - 1.42154*aest + 0.000396)
     else corr.est = 0
  }
  else if(family==2){
    info = t.coefs[t.coefs$df==tdf,]
    if(aest>info[6]) corr.est = as.numeric(sign(est)*(info[5]*aest^3 + info[4]*aest^2 + info[3]*aest + info[2]))
    else corr.est = 0
  }
  else if(family==3) corr.est = max(-0.5301586 + 0.716288*aest, 0.1)
  else if(family==23) corr.est = min(0.5301586 - 0.716288*aest, -0.1)
  else if(family==4) corr.est = max(0.02220367 + 0.70569984*aest, 1)
  else if(family==24) corr.est = min(-0.02220367 - 0.70569984*aest, -1)
  else corr.est = est
  
  return(corr.est)
}

CopEst = function(u, v, family, method, maxit=50){
  est = BiCopEst(u, v, family=family, method='mle')
  if(family==2){
    fl = floor(est$par2)
    dfval = ifelse(est$par2-fl<0.5,fl,fl+1)
    intdf = ifelse(dfval<3,3,dfval)
  }
  else intdf=NA
  
  if(method=='mle')
    return(list(par=est$par, par2=intdf))
  else{
    if(family==1||family==2){
      if(sign(cor(u,v,method='kendall'))>0)
        interval = c(0,.99)
      else
        interval=c(-.99,0)
      }
    else if(family==3) interval = c(0.1,10*est$par)
    else if(family==23) interval = c(10*est$par, -0.1)
    else if(family==4) interval = c(1,10*est$par)
    else interval = c(10*est$par, -1)
    
    iterating=TRUE; iter=1
    while(iterating){
      
      theta.list = seq(interval[1], interval[2], length.out=10)    
      dep.list = rep(0,10)
      for(j in 1:10){
        if(family==2){
          del.loglik = BiCopDeriv(u,v,family=family,par=theta.list[j], , par2=intdf, deriv="par") / 
            BiCopPDF(u,v,family=family,par=theta.list[j], par2=intdf)
        }
        else{
          del.loglik = BiCopDeriv(u,v,family=family,par=theta.list[j], deriv="par") / 
            BiCopPDF(u,v,family=family,par=theta.list[j])
        }      
        
        dep.list[j] =(min(sum(del.loglik>0), sum(del.loglik<0)) + sum(del.loglik==0))/length(u)
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
    
    if(maxdep<0.5 && iter==maxit)
      warning("Maximum number of iter ations reached: estimator does not have maximum possible likelihood depth")
    return(list(par=MD.bias(theta.list[which.max(dep.list)], family=family, tdf=intdf), par2=intdf))
  }
}

getBestCop = function(u,v, fam.list, method){
  
  # Obtain estimates for family = 1,2,(3,4) or (23,24)
  params = matrix(rep(NA,8),ncol=2)
  for(i in 1:4){
    e = CopEst(u,v, family=fam.list[i], method=method)
    params[i,] = c(e$par, e$par2)
  }
  
  # Generate empirical copula for the (pseudo-)data
  U = cbind(u,v)
  empC = C.n(U, U=U)
  
  # Compare them using D2 distance
  D2.list = rep(0,4)
  for(i in 1:4){
    if(i==2)
      D2.list[i] = sum((pCopula(U, tCopula(params[i,1],df=params[i,2])) - empC)^2)
    else
      D2.list[i] = sum((BiCopCDF(u,v,family=fam.list[i],par=params[i,1]) - empC)^2)
  }
  
  bestC = which.min(D2.list)
  # takes care of the case t30 is selected as best copula
  # then selects Gaussian copula instead
  if(bestC==2 && params[bestC,2]==30)
    bestC = 1
  return(list(family=fam.list[bestC], par=params[bestC,1], par2=params[bestC,2]))
}

CopNetwork = function(pdata, is.U=TRUE, method, quiet=TRUE){
  if(!is.U) U = pobs(pdata)
  else U = pdata
  
  n = nrow(U); p = ncol(U)
  npairs = p*(p-1)/2
  out.mat = matrix(rep(NA,5*npairs), ncol=5)
  sum.till.i = 0
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      out.mat[sum.till.i+(j-i), 1:2] = c(i,j)
    }
    sum.till.i = sum.till.i + (p-i)
  }
  
  for(i in 1:npairs){
    u = U[,out.mat[i,1]]; v = U[,out.mat[i,2]]
    if(!quiet)
      cat(out.mat[i,1], out.mat[i,2], "\n")
    if(BiCopIndTest(u,v)$p.value < 0.05){
      if(cor(u,v,method='kendall')>=0)
        fam.list=c(1,2,3,4)
      else
        fam.list=c(1,2,23,24)
      bestC = getBestCop(u,v, fam.list=fam.list, method=method)
      out.mat[i,3:5] = c(bestC$family, bestC$par, bestC$par2)
    }  
  }
  
  return(graph=out.mat)
}

plot.CopNetwork = function(NetworkMat, varnames, titletext, leg=T){
  index.list = which(is.na(NetworkMat[,3])==F)
  
  # setting colors for each edge
  col.list = rep(0,length(index.list))
  for(i in 1:length(index.list)){
    posi = index.list[i]
    if(NetworkMat[posi,3]==1) col.list[i] = 'black'
    else if (NetworkMat[posi,3]==2) col.list[i] = 'green'
    else if (NetworkMat[posi,3]==3) col.list[i] = 'blue'
    else if (NetworkMat[posi,3]==23) col.list[i] = 'cyan'
    else if (NetworkMat[posi,3]==4) col.list[i] = 'red'
    else if (NetworkMat[posi,3]==24) col.list[i] = 'coral'
  }
  graphdata = data.frame(
    from = NetworkMat[index.list,1],
    to = NetworkMat[index.list,2]
  )
  
  #par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(2,2,2,2))
  plot(mlgraph <- qgraph(graphdata, edge.color = col.list, directed=F, labels = varnames))
  title(paste(titletext),line = 0)
  if(leg){
  legend(-1.1,-.5, cex=.9, title="Copula type",
         legend=c("Gaussian", "t", "Clayton", "90° rotated Clayton", "Gumbel", "90° rotated Gumbel"),
         lty=1, lwd=2, col=c('black','green','blue','cyan','red','coral'))}
}

summary.CopNetwork = function(graph){
  ## within-group and between-group graphs
  npairs = nrow(graph)
  group.conns = cbind(matrix(rep(0,npairs*2),ncol=2),graph[,3])
  for(i in 1:npairs){
    group.conns[i,1] = ifelse(graph[i,1]<=6, 1, ifelse(graph[i,1]<=10, 2, 3))
    group.conns[i,2] = ifelse(graph[i,2]<=6, 1, ifelse(graph[i,2]<=10, 2, 3))
  }
  
  # intra-group connections
  type1 = group.conns[which(group.conns[,1]==group.conns[,2]),]; type1[is.na(type1[,3]),3] = 0
  num1 = cbind(c(11,22,33),matrix(rep(0,21),ncol=7))
  cop.list=c(0,1,2,3,23,4,24)
  for(i in 1:3){
    for(j in 1:7)
      num1[i,j+1] = sum(type1[type1[,3]==cop.list[j],1]==i)
  }
  colnames(num1) = c("Type",as.character(cop.list))
  
  # inter-group connections
  type2 = group.conns[which(group.conns[,1]!=group.conns[,2]),]; type2[is.na(type2[,3]),3] = 0
  type2 = cbind(as.numeric(paste(type2[,1],type2[,2],sep="")),type2[,3])
  num2 = cbind(c(12,13,23),matrix(rep(0,21),ncol=7))
  for(i in 1:3){
    for(j in 1:7)
      num2[i,j+1] = sum(type2[type2[,2]==cop.list[j],1]==num2[i,1])
  }
  return(rbind(num1,num2))
}