library(Matrix)

# define the BIC function
bic1 <- function(Y,add_x,U,thetahat,bHat,alphahat,Cutoff,tau){
  output <- 0
  n <- length(Y)
  res <- Y- add_x%*%thetahat - U%*%bHat - rep(alphahat,n)
  p <- ncol(U)
  loss1 <- sum(rho(res,tau))
  loss2 <- sum(srho(res,tau,h))
  bZero = (abs(bHat) < Cutoff)
  bNonZero = !bZero
  output <- c(log(loss1) + sum(bNonZero)*log(n)/(2*n)*(2^-1),
              log(loss2) + sum(bNonZero)*log(n)/(2*n)*(2^-1))
  return(output)
}


# candidate values for the two penalties
cand1 <- c(exp(-8),exp(-9),exp(-10),exp(-11),exp(-12))
cand2 <- c(exp(-8),exp(-9),exp(-10),exp(-11),exp(-12))

# quantiles of interest
quans <- c(2:8)/10
est1 <- c()
est2 <- c()
c1 <- 1e-3
ncv <- 10
# random seed for generating the data
tmpid <- 1
for (t in quans){
  tau <- t
  id <- tmpid
  set.seed(id)
  # sample size
  nobs <- 500
  cvresult1 <- c()
  source('simple_DataGen.R')
  add_xx <- add_x
  # initial fit
  gammacv <- c()
  for (gamma in cand2){
    cvloss <- 0
    for (cv in 1:ncv){
      lc <- length(daty)/ncv
      tind <- (cv-1)*lc + c(1:lc)
      trind <- c(1:length(daty))[-tind]
      X <- mtx[trind,]
      Y <- daty[trind]
      add_tmp <- add_x[trind,]
      source('initial_fit.R')
      U = GetDesignMatrix2(X=mtx[tind,],beta.basis=beta.basis)
      cvres <- daty[tind]- add_x[tind,]%*%thetahat - U%*%bHat - rep(alphahat,lc)
      cvloss <- cvloss + sum(rho(cvres,tau))
    }
    gammacv <- c(gammacv, cvloss)
  }
  optid <- which(gammacv==min(gammacv))[1]
  X = mtx
  Y = daty
  add_tmp = add_x
  gamma <- cand2[optid]
  n <- nrow(mtx)
  source('initial_fit.R')
  initial_est <- c(alphahat, thetahat, bHat)
  
  for (l in 1:length(cand1))
  {
    for (g in 1:length(cand2))
    {
      lambda = cand1[l]
      gamma = cand2[g]
      X = mtx
      Y = daty
      source('CLoSE_fit.R')
      gammacv <- c()
      UU <- U
      for (gamma in cand2){
        cvloss <- 0
        for (cv in 1:ncv){
          lc <- length(daty)/ncv
          tind <- (cv-1)*lc + c(1:lc)
          trind <- c(1:length(daty))[-tind]
          X <- mtx[trind,]
          Y <- daty[trind]
          add_x <- add_xx[trind,]
          U <- UU[trind,]
          source('step2_fit.R')
          U = GetDesignMatrix2(X=mtx[tind,],beta.basis=beta.basis)
          cvres <- daty[tind]- add_xx[tind,]%*%thetahat - U%*%bHat - rep(alphahat,lc)
          cvloss <- cvloss + sum(srho(cvres,tau,h))
        }
        gammacv <- c(gammacv, cvloss)
      }
      optid <- which(gammacv==min(gammacv))[1]
      print(optid)
      X = mtx
      Y = daty
      add_x <- add_xx
      U <- UU
      gamma = cand2[optid]
      source('step2_fit.R')
      gammacv1 <- bic1(Y,add_x,U,thetahat,bHat,alphahat,Cutoff,tau)
      cvresult1 <- rbind(cvresult1, c(l,g,gammacv1,sum(bZero)))
      saveRDS(object=cvresult1,file=paste0(id,'cvresult1',tau,'.RDS'))
    }
  }
  tmpcv <- cvresult1[cvresult1[,3]==min(cvresult1[,3]),]
  print(tmpcv)
  if (is.null(nrow(tmpcv))){
    cvout <- tmpcv[1:3]
  }else{
    cvout <- tmpcv[1,1:3]
  }
  lambda <- cand1[cvout[1]]
  gamma <- cand2[cvout[2]]
  source('CLoSE_fit.R')
  est1 <- as.matrix(bHat)
  
  saveRDS(object=est1,file=paste0(id,'est1',tau,'.RDS'))
}
saveRDS(object=beta.basismat,file='basismat.RDS')
