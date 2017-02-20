#must run saratoga.r first

library(mvtnorm)

y <- SaratogaHouses[,1]
n <- nrow(X)
p <- ncol(X) #includes column of 1s
it <- 10000 #number of iterations

#create data structures for the parameters
beta <- matrix(ncol=it,nrow=p)#betas stored as columns
omega <- vector(mode="numeric",length=it)
Lambda <- matrix(ncol=it,nrow=n)

#set initial values for parameters and hyperparameters
beta[,1] <- rep(0,p)
omega[1] <- 1
Lambda[,1] <- rep(1,n)
h <- 2
m <- rep(0,p)
K <- diag(1,p,p)
d <- 2
eta <- 2

d <- n+d#only update once

for (i in 2:it) {
  #update beta
  tXLX <- t(X) %*% diag(Lambda[,i-1]) %*% X
  betahat <- solve(tXLX) %*% t(X) %*% diag(Lambda[,i-1]) %*% y
  Kstar <- K + tXLX
  m <- solve(Kstar) %*% (Kstar %*% m + tXLX %*% betahat)
  beta[,i] <- rmvnorm(1,mean=m,sigma=solve(omega[i-1]*Kstar))
  #update omega
  yhat <- X %*% betahat
  Sy <- t(y - yhat) %*% diag(Lambda[,i-1]) %*% (y - yhat)
  etastar <- eta + Sy + t(m-betahat) %*% solve(solve(K) + solve(tXLX)) %*% (m-betahat)
  omega[i] <- rgamma(1,d/2,etastar/2)
  #update Lambda
  Lambda[,i] <- rgamma(1,(h+1)/2,(omega[i]/2)*t(y-X %*% beta[,i]) %*% (y-X %*% beta[,i])+h/2)
}