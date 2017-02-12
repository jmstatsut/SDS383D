library(mlbench)

ozone = data(Ozone, package='mlbench')

# ?Ozone

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

# compute the estimator
betahat = solve(t(x) %*% x) %*% t(x) %*% y

# Fill in the blank
# betacov = ?

# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as an extra column)
# lm1 = lm(y~x-1)

# summary(lm1)
# betacovlm = vcov(lm1)
# sqrt(diag(betacovlm))

#X: design matrix
#y: observations
#B: how many bootstrap iterations to use
compute_bootstrap_cov <- function(x,y,B) {
  n <- nrow(x)
  p <- ncol(x)
  betahat_boot <- matrix(nrow=B,ncol=p)
  betahat = solve(t(x) %*% x) %*% t(x) %*% y
  yhat <- x %*% betahat
  residuals <- y-yhat
  #create bootstrap samples for beta estimation
  for(i in 1:B) {
    boot_sample_indexes <- sample(1:n,n,replace=TRUE)
    xtemp <- x[boot_sample_indexes,]
    ytemp <- y[boot_sample_indexes]
    betahat_boot[i,] <- solve(t(xtemp) %*% xtemp) %*% t(xtemp) %*% ytemp
  }
  #find covariance matrix of bootstrap samples
  return(cov(betahat_boot))
}

cov_boot <- compute_bootstrap_cov(x,y,10000)