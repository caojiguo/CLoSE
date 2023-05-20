# define BIC function
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


mtx <- readRDS('soy_x.RDS')
daty <- readRDS('soy_y.RDS')
add_var <- readRDS('soy_add.RDS')
add_var[,3] <- add_var[,3]/add_var[,2]
add_xx <- add_x <- add_var[,-2]
names(add_x) <- c('prcp','irr_ratio')
mtx2 = mtx[,30:335] # minimum temperature
mtx1 = mtx[,395:700] # maximum temperature
t_grid = c(0,1)
Y = daty


cand2 <- cand1 <- c(exp(-5),exp(-6))
cand4 <- cand3 <- c(exp(-6),exp(-7),exp(-8),exp(-9))
cand6 <- cand5 <- c(exp(-6),exp(-7),exp(-8),exp(-9))
cvresult1 <- c()
bhat1s <- c()
bhat2s <- c()
thetas <- c()
alphas <- c()
bind <- c()
tau <- 0.5
ncv <- 10

gammacv <- c()
for (l3 in 1:length(cand3)){
  for (l4 in 1:length(cand4)){
    cvloss <- 0
    for (cv in 1:ncv){
      lc <- length(daty)/ncv
      lc <- round(lc)
      tind <- (cv-1)*lc + c(1:lc)
      if (max(tind)>nrow(mtx1)){
        tind <- c(min(tind):nrow(mtx1)) 
      }
      trind <- c(1:length(daty))[-tind]
      X1 <- mtx1[trind,]
      X2 <- mtx2[trind,]
      Y <- daty[trind]
      add_tmp <- add_x[trind,]
      gamma <- c(cand3[l3], cand4[l4])
      source('initial_fit_2fvars.R')
      U1 = GetDesignMatrix2(X=mtx1[tind,],beta.basis=beta.basis)
      U2 = GetDesignMatrix2(X=mtx2[tind,],beta.basis=beta.basis)
      U <- cbind(U1,U2)
      cvres <- daty[tind]- add_x[tind,]%*%thetahat - U%*%bHat - rep(alphahat,length(tind))
      cvloss <- cvloss + sum(rho(cvres,tau))
    }
    gammacv <- rbind(gammacv, c(l3,l4,cvloss))
  }
}
optid <- which(gammacv[,3]==min(gammacv[,3]))[1]
print(optid)

X1 = mtx1
X2 = mtx2
Y = daty
add_tmp <- add_xx <- add_x
gamma <- c(cand3[gammacv[optid,1]],cand4[gammacv[optid,2]])
n <- nrow(mtx1)
source('initial_fit_2fvars.R')
initial_est <- c(alphahat, thetahat, bHat)

for (l1 in 1:length(cand1))
{
  for (l2 in 1:length(cand2))
  {
    for (l3 in 1:length(cand3))
    {
      for (l4 in 1:length(cand4))
      {
        gammacv <- c()
        lambda = c(cand1[l1],cand2[l2])
        gamma = c(cand3[l3],cand4[l4])
        X1 = mtx1
        X2 = mtx2
        Y = daty
        add_x <- add_xx
        source('CLoSE_fit_2fvars.R')
        if (l1+l2+l3+l4>4&sum(bZero)==110){
          cvresult1 <- rbind(cvresult1, cvresult1[1,])
          saveRDS(object = cvresult1, file = 'cvresult1.RDS')
          bhat1s <- rbind(bhat1s, bhat1s[1,])
          bhat2s <- rbind(bhat2s, bhat2s[1,])
          saveRDS(object = bhat1s, file = 'bhat1s.RDS')
          saveRDS(object = bhat2s, file = 'bhat2s.RDS')
        }else{
          U_full <- U
          for (l5 in 1:length(cand5)){
            for (l6 in 1:length(cand6)){
              cvloss <- 0
              for (cv in 1:ncv){
                lc <- length(daty)/ncv
                lc <- round(lc)
                tind <- (cv-1)*lc + c(1:lc)
                if (max(tind)>nrow(mtx1)){
                  tind <- c(min(tind):nrow(mtx1)) 
                }
                trind <- c(1:length(daty))[-tind]
                U <- U_full[trind,]
                Y <- daty[trind]
                add_x <- add_xx[trind,]
                gamma <- c(cand5[l5], cand6[l6])
                source('step2_fit_2fvars.R')
                U1 = GetDesignMatrix2(X=mtx1[tind,],beta.basis=beta.basis)
                U2 = GetDesignMatrix2(X=mtx2[tind,],beta.basis=beta.basis)
                U <- cbind(U1,U2)
                cvres <- daty[tind]- add_xx[tind,]%*%thetahat - U%*%bHat - rep(alphahat,length(tind))
                cvloss <- cvloss + sum(rho(cvres,tau))
              }
              gammacv <- rbind(gammacv, c(l5,l6,cvloss))
            }
          }
          optid <- which(gammacv[,3]==min(gammacv[,3]))[1]
          gamma <- c(cand5[gammacv[optid,1]],cand6[gammacv[optid,2]])
          U <- U_full
          Y <- daty
          add_x <- add_xx
          source('step2_fit_2fvars.R')
          gammacv1 <- bic1(Y,add_x,U,thetahat,bHat,alphahat,Cutoff,tau)
          cvresult1 <- rbind(cvresult1, c(l1,l2,l3,l4,gammacv[optid,1],gammacv[optid,2],gammacv1,sum(bZero)))
          saveRDS(object = cvresult1, file = 'cvresult1.RDS')
          bhat1s <- rbind(bhat1s, bHat1)
          bhat2s <- rbind(bhat2s, bHat2)
          saveRDS(object = bhat1s, file = 'bhat1s.RDS')
          saveRDS(object = bhat2s, file = 'bhat2s.RDS')
          thetas <- rbind(thetas,thetahat)
          alphas <- rbind(alphas, alphahat)
          print(c(l1,l2,l3,l4))
        }
      }
    }
  }
}
