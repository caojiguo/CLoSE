if (sum(bZero)!=length(bHat)){
  bHat[bZero] = 0
  bNonZero = !bZero
  U1 = as.matrix(U[,bNonZero])
  V1 = V[bNonZero,bNonZero]

  # print(dim(V1))

  lambda <- 0
  # cand2 <- c(exp(-8),exp(-9),exp(-10),exp(-11),exp(-12))
  # gamma = cand2[cvout]
  full_x <- cbind(add_x, U1)
  tmp <- rq(Y~ full_x,tau)
  bHat[bNonZero] <- tmp$coefficients[-c(1:(nonfunc+1))]
  thetahat <-  tmp$coefficients[2:(nonfunc+1)]
  alphahat <- tmp$coefficients[1]
  change <- 1
  n <- nrow(full_x)
  while(change>=tol)
  {
    betaold <- c(thetahat,bHat)
    tmpchange <- 1
    tmpobj2 <- tmpS(Y,add_x,U,V,thetahat,bHat,alphahat,L2NNer,lambda,tau,M,gamma,h)
    while(tmpchange>=tol)
    {
      eta <- 2^0
      tmpobj1 <- tmpobj2
      ds = 0
      s = 0
      for (ii in 1:n){
        ds <- ds + dGh(Y[ii]-sum(add_x[ii,]*thetahat)-sum(U1[ii,]*bHat[bNonZero])-alphahat,h)*as.matrix(c(1,add_x[ii,],U1[ii,]))%*%c(1,add_x[ii,],U1[ii,])
        s <- s-(Gh(Y[ii]-sum(add_x[ii,]*thetahat)-sum(U1[ii,]*bHat[bNonZero])-alphahat,h)+tau-1)*c(1,add_x[ii,],U1[ii,])}
      s <- as.matrix(s)
      tmpmatrix <- bdiag(matrix(0,nrow=nonfunc+1, ncol=nonfunc+1),gamma*V1)
      tmpmatrix <- tmpmatrix + ds
      tmp <- svd(tmpmatrix)
      tmpinv <- (tmp$u)%*%diag((tmp$d)^-1)%*%t(tmp$v)
      scjk <- t(bHat[bNonZero])%*%V1*gamma
      sdiff <- -1
      while(sdiff < 0)
      {
        tmpstep <- eta*tmpinv%*%(s+as.matrix(c(rep(0,nonfunc+1),scjk)))
        tmpbHat <- rep(0,length(bHat))
        tmpbHat[bNonZero] <- bHat[bNonZero] - tmpstep[-c(1:(nonfunc+1))]
        tmptheta <- thetahat - tmpstep[2:(nonfunc+1)]
        tmpalphahat <- alphahat - tmpstep[1]
        tmpchange <- mean(tmpstep^2)
        tmpobj2 <- tmpS(Y,add_x,U,V,tmptheta,tmpbHat,tmpalphahat,L2NNer,lambda,tau,M,gamma,h)
        sdiff <- (tmpobj1-tmpobj2)/n
        eta <- eta/2
      }
      bHat <- tmpbHat
      alphahat <- tmpalphahat
      thetahat <- tmptheta
  }
  change <- mean((c(thetahat,bHat)-betaold)^2)
  }
}else{
  bHat <- rep(0,length(bHat))
}

