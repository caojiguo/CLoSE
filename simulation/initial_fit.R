library(psych)
library(pls)
library(plsdof)
library(fda)
library(MASS)
library(quantreg)
# quantile loss
rho <- function(u,tau)
{
  u*(tau - (u<=0))
}

le <- function(u)
{
  (3*u^2/4-u^4/8+3/8)*(abs(u)<=1) + abs(u)*(abs(u)>1)
}

lg <- function(u)
{
  (2/pi)^0.5*exp(-u^2/2)+ u*(1-2*pnorm(-u))
}


# fscad
fscad <- function(u,lambda,a)
{
  out <- 0
  if (u<=lambda){
    out <- lambda*u
  }else if (u < a*lambda){
    -(u^2 - 2*a*lambda*u+lambda^2)/(2*(a-1))
  }else{
    (a+1)*lambda^2/2
  }
}

# objective function
tmpS <- function(Y,add_tmp,U,V,thetahat,bHat,alphahat,L2NNer,lambda,tau,M,gamma,h)
{
  n <- length(Y)
  m <- length(betaNormj)
  loss <- 0
  pen <- 0
  for (ii in 1:n)
  {
    tmp <- Y[ii] -t(add_tmp[ii,])%*%thetahat-t(U[ii,])%*%bHat-alphahat
    loss <- loss + srho(tmp,tau,h)
  }
  lqaWk = matrix(0,L,L)
  for(j in 1:M)
  {
    index = j:(j+d)
    betaNormj[j] = t(bHat[j:(j+d)])%*%W[,,j]%*%bHat[j:(j+d)]
    cjk = Dpfunc(betaNormj[j]*L2NNer,lambda,a)
    scjk[index] = scjk[index] + cjk*(L2NNer/(epsilon + betaNormj[j]))*t(bHat[index])%*%W[,,j] 
    lqaWk[index,index] = lqaWk[index,index] + cjk*(L2NNer/(epsilon + betaNormj[j]))*W[,,j]
  }
  lqaW = lqaWk
  lqaW = lqaW/2
  out <- as.numeric(t(bHat)%*%(n*lqaW)%*%bHat + loss + t(bHat)%*%V%*%bHat*gamma)
  return(out)
}

# calculate the design matrix after B-spline expansion
GetDesignMatrix2 = function(X,beta.basis)
{
  rng     = getbasisrange(beta.basis)
  breaks  = c(rng[1],beta.basis$params,rng[2])
  Mobs    = dim(X)[2]
  tobs    = seq(rng[1],rng[2],length.out = Mobs)
  
  beta.basismat = eval.basis(tobs,beta.basis)
  hDM           = (rng[2]-rng[1])/Mobs
  cefDM         = c(1, rep(c(4,2), (Mobs-2)/2), 1)
  U             = hDM/3*X%*%diag(cefDM)%*%beta.basismat
  return(U=U)
}

# weights matrix of LQA approximation
slos.compute.weights = function(basis)
{
  L       = basis$nbasis
  rng     = getbasisrange(basis)
  breaks  = c(rng[1],basis$params,rng[2])
  M       = length(breaks) - 1
  norder  = L-M+1
  W       = array(0,dim=c(norder,norder,M))
  
  for (j in 1:M)
  {
    temp = inprod(basis,basis,rng=c(breaks[j],breaks[j+1]))
    W[,,j] = temp[j:(j+norder-1),j:(j+norder-1)]
  }
  W
}

# derivative of fSCAD
Dpfunc = function(u,lambda,a)
{
  if (u<=lambda){Dpval = lambda
  }else if(u<a*lambda){
    Dpval = -(u-a*lambda)/(a-1)
  }else{Dpval = 0}
  Dpval
}


Gh = function(x,h){
  out <- pnorm(x/h)
}

# Gaussian kernel with bandwidth h
dGh = function(x,h){
  out <- dnorm(x/h)
  return(out)
}

# smoothed quantile loss
srho <- function(u,tau,h)
{
  h/2*lg(u/h) + (tau-0.5)*u
}

a <- 3.7
Cutoff <- 10^-3
tol <- 10^-5

# define B-spline basis functions
domain <- range(t_grid)
n = length(Y)
m = dim(X)[2]
M = 50
d = 5
norder   = d+1
knots    = seq(domain[1],domain[2], length.out=M+1)
nknots   = length(knots)
nbasis   = length(knots) + norder - 2
beta.basis  = create.bspline.basis(knots,nbasis,norder)
rng     = getbasisrange(beta.basis)
breaks  = c(rng[1],beta.basis$params,rng[2])
Mobs    = dim(X)[2]
tobs    = seq(rng[1],rng[2],length.out = Mobs)
beta.basismat = eval.basis(tobs,beta.basis)
L2NNer = sqrt(M/(rng[2]-rng[1]))  

U = GetDesignMatrix2(X=X,beta.basis=beta.basis)
W = slos.compute.weights(beta.basis)
L = beta.basis$nbasis
V = eval.penalty(beta.basis,int2Lfd(2))
n = length(Y)
h <- (M/n)^0.4

full_x <- cbind(add_tmp, U)
nonfunc <- ncol(add_tmp)
if (is.null(nonfunc)){
  nonfunc <- 1
  add_tmp <- as.matrix(add_tmp)
}
# no fSCAD when calculating the initial estimates
lambda = 0

tmp <- rq(Y~ full_x,tau)
bHat <- tmp$coefficients[-c(1:(nonfunc+1))]
thetahat <-  tmp$coefficients[2:(nonfunc+1)]
alphahat <- tmp$coefficients[1]
minabs <- sort(abs(c(thetahat,bHat)))[which(sort(abs(c(thetahat,bHat)))>0)[1]]
epsilon <- max(tol/(2*n*Dpfunc(0,lambda,a))*minabs,1e-5)

betaNormj = rep(0,M)
obj1 <- rep(0,3)
obj2 <- rep(1,3)
change <- 1
while(change>=tol)
{
  betaold <- c(thetahat,bHat)
  lqaW = NULL
  lqaWk = matrix(0,L,L)
  scjk = rep(0,L)
  tmpchange <- 1
  tmpobj2 <- tmpS(Y,add_tmp,U,V,thetahat,bHat,alphahat,L2NNer,lambda,tau,M,gamma,h)
  while(tmpchange>=tol)
  {
    eta <- 2^0
    tmpobj1 <- tmpobj2
    ds = 0
    s = 0
    for (ii in 1:n){
      ds <- ds + dGh(Y[ii]-sum(add_tmp[ii,]*thetahat)-sum(U[ii,]*bHat)-alphahat,h)*as.matrix(c(1,add_tmp[ii,],U[ii,]))%*%c(1,add_tmp[ii,],U[ii,])
      s <- s-(Gh(Y[ii]-sum(add_tmp[ii,]*thetahat)-sum(U[ii,]*bHat)-alphahat,h)+tau-1)*c(1,add_tmp[ii,],U[ii,])}
    s <- as.matrix(s)
    tmpmatrix <- bdiag(matrix(0,nrow=nonfunc+1, ncol=nonfunc+1),gamma*V)
    tmpmatrix <- tmpmatrix + ds
    tmp <- svd(tmpmatrix)
    tmpinv <- (tmp$u)%*%diag((tmp$d)^-1)%*%t(tmp$v)
    scjk <- t(bHat)%*%V*gamma
    sdiff <- -1
    # tune the step size of each iteration to make sure that the value of the objective function decreases after each iteration
    while(sdiff < 0)
    {
      tmpstep <- eta*tmpinv%*%(s+as.matrix(c(rep(0,nonfunc+1),scjk)))
      dt <- -tmpinv%*%(s+as.matrix(c(rep(0,nonfunc+1),scjk)))
      dt <- as.matrix(dt)
      tmpbHat <- rep(0,length(bHat))
      tmpbHat <- bHat - tmpstep[-c(1:(nonfunc+1))]
      tmptheta <- thetahat - tmpstep[2:(nonfunc+1)]
      tmpalphahat <- alphahat - tmpstep[1]
      tmpchange <- mean(tmpstep^2)
      tmpobj2 <- tmpS(Y,add_tmp,U,V,tmptheta,tmpbHat,tmpalphahat,L2NNer,lambda,tau,M,gamma,h)
      sdiff <- tmpobj1-tmpobj2 + c1*eta*t(dt)%*%(s+as.matrix(c(rep(0,nonfunc+1),scjk)))
      eta <- eta/2
    }
    bHat <- tmpbHat
    alphahat <- tmpalphahat
    thetahat <- tmptheta
  }
  change <- mean((c(thetahat,bHat)-betaold)^2)
}
