tmpS <- function(Y,add_x,U,V,thetahat,bHat1,bHat2,alphahat,tau,gamma,h)
{
  bHat <- c(bHat1,bHat2)
  n <- length(Y)
  loss <- 0
  pen <- 0
  for (ii in 1:n)
  {
    tmp <- Y[ii] -t(add_tmp[ii,])%*%thetahat-t(U[ii,])%*%bHat-alphahat
    loss <- loss + srho(tmp,tau,h)
  }
  for (jj in 1:2){
    for(j in 1:M + (jj-1)*(M+d))
    {
      index = j:(j+d)
      betaNormj[j] = t(bHat[j:(j+d)])%*%W[,,j-(jj-1)*(M+d)]%*%bHat[j:(j+d)]
      cjk = Dpfunc(betaNormj[j]*L2NNer,lambda[jj],a)
      scjk[index] = scjk[index] + cjk*(L2NNer/(epsilon + betaNormj[j]))*t(bHat[index])%*%W[,,j-(jj-1)*(M+d)] 
      lqaWk[index,index] = lqaWk[index,index] + cjk*(L2NNer/(epsilon + betaNormj[j]))*W[,,j-(jj-1)*(M+d)]
    }
  }
  
  lqaW = lqaWk
  lqaW = lqaW/2
  
  out1 <- t(bHat)%*%(n*lqaW)%*%bHat + t(bHat1)%*%V%*%bHat1*gamma[1] + t(bHat2)%*%V%*%bHat2*gamma[2] + loss
  out2 <- loss
  return(c(out1,out2))
}

if (sum(bZero)!=length(bHat)){
  bHat[bZero] = 0
  bNonZero = !bZero
  UU = as.matrix(U[,bNonZero])
  VV = bdiag(gamma[1]*V,gamma[2]*V)[bNonZero,bNonZero]

  # print(dim(VV))

  lambda <- c(0,0)
  full_x <- cbind(add_x, UU)
  tmp <- rq(Y~ full_x,tau)
  bHat[bNonZero] <- tmp$coefficients[-c(1:(nonfunc+1))]
  bHat1 <- bHat[1:ncol(U1)]
  bHat2 <- bHat[1:ncol(U2)+ncol(U1)]
  thetahat <-  tmp$coefficients[2:(nonfunc+1)]
  alphahat <- tmp$coefficients[1]
  change <- 1
  while(change>=tol)
  {
    betaold <- c(thetahat,bHat)
    tmpchange <- 1
    tmpobj2 <- tmpS(Y,add_x,U,V,thetahat,bHat1,bHat2,alphahat,tau,gamma,h)
    while(tmpchange>=tol)
    {
      eta <- 2^0
      tmpobj1 <- tmpobj2
      ds = 0
      s = 0
      n <- length(Y)
      for (ii in 1:n){
        ds <- ds + dGh(Y[ii]-sum(add_x[ii,]*thetahat)-sum(UU[ii,]*bHat[bNonZero])-alphahat,h)*as.matrix(c(1,add_x[ii,],UU[ii,]))%*%c(1,add_x[ii,],UU[ii,])
        s <- s-(Gh(Y[ii]-sum(add_x[ii,]*thetahat)-sum(UU[ii,]*bHat[bNonZero])-alphahat,h)+tau-1)*c(1,add_x[ii,],UU[ii,])}
      s <- as.matrix(s)
      tmpmatrix <- bdiag(matrix(0,nrow=nonfunc+1, ncol=nonfunc+1),VV)
      tmpmatrix <- tmpmatrix + ds
      tmp <- svd(tmpmatrix)
      tmpinv <- (tmp$u)%*%diag((tmp$d)^-1)%*%t(tmp$v)
      scjk <- t(bHat[bNonZero])%*%VV
      sdiff <- -1
      while(sdiff < 0)
      {
        scjk <- as.matrix(scjk)
        tmpstep <- eta*tmpinv%*%(s+as.matrix(c(rep(0,nonfunc+1),scjk)))
        tmpbHat <- rep(0,length(bHat))
        tmpbHat[bNonZero] <- bHat[bNonZero] - tmpstep[-c(1:(nonfunc+1))]
        tmpbHat1 <-  tmpbHat[1:ncol(U1)]
        tmpbHat2 <-  tmpbHat[1:ncol(U2)+ncol(U1)]
        tmptheta <- thetahat - tmpstep[2:(nonfunc+1)]
        tmpalphahat <- alphahat - tmpstep[1]
        tmpchange <- mean(tmpstep^2)
        tmpobj2 <- tmpS(Y,add_x,U,V,tmptheta,tmpbHat1,tmpbHat2,tmpalphahat,tau,gamma,h)
        sdiff <- ((tmpobj1[1]-tmpobj2[1])/n)[1]
        eta <- eta/2
      }
      bHat1 <- tmpbHat1
      bHat2 <- tmpbHat2
      alphahat <- tmpalphahat
      thetahat <- tmptheta
      bHat <- c(bHat1,bHat2)
  }
  change <- mean((c(thetahat,bHat)-betaold)^2)
  }
}else{
  bHat <- rep(0,length(bHat))
}

