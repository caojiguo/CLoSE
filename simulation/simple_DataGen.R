library(e1071)

# time points
t_grid <- newgrids <- seq(-1,1,length.out = 400)
tmp2 <- 'sin(2*pi*t)*(t>=-0.5)*(t<=0.5)'
text <- paste0('function(t)',tmp2)
# the true slope function
rhos <- eval(parse(text=text))

mtx <- matrix(0,nobs,length(t_grid))
add_x <- matrix(rnorm(nobs*2,0,0.1),nrow=nobs)
daty <- c()
Mobs    = ncol(mtx)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM     = (tmax-tmin)/Mobs
cefDM   = c(1, rep(c(4,2), (Mobs-2)/2), 1)
set.seed(1)
for (j in 1:nobs)
{
  X <- rwiener(end=2,frequency = 200)
  mtx[j,] <- c(X)
}
set.seed(id+1)
err <- rnorm(nobs,0,0.02)
daty <- hDM/3*mtx%*%diag(cefDM)%*%rhos(t_grid)
daty <- daty + add_x%*%c(1,1) + err
