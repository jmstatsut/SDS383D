rm(list=ls())
setwd("C:/Users/Jesse/Box Sync/!Jesse/Stats/Courses/2017_1/SDS 383D Statistical Modeling II/Exercise Solutions/04/")
library(mvtnorm)

cheese_data <- read.csv("04_02_cheese.csv",header=TRUE)
stores <- unique(cheese_data[,1])
m <- length(stores)
n <- nrow(cheese_data)
log_cheese_data <- cheese_data
log_cheese_data[,2:3] <- log(log_cheese_data[,2:3])


# Matrix of the fixed effects
X <- model.matrix(~ price + disp + disp:price, data = log_cheese_data)
# Matrix of the fixed effects
Z <- model.matrix(~ price + disp + disp:price, data = log_cheese_data)
p <- ncol(X)
q <- ncol(X)

# Store data and design matrices in a list
yi <- list()
Xi <- list()
Zi <- list()
for (i in 1:m){	
  yi[[i]] <- log_cheese_data[log_cheese_data[,1] == stores[i], 3]
  Xi[[i]] <- X[log_cheese_data[,1] == stores[i],]
  Zi[[i]] <- Z[log_cheese_data[,1] == stores[i],]
}


Gibbs_size <- 11000
burnin <- 1000
thin <- 2

# Create structures for parameters and initiaize their values
sigmasq <- rep(1,Gibbs_size)
B <- matrix(1, nrow = p, ncol = Gibbs_size)
g <- array(0, dim = c(m, q, Gibbs_size))
Tinv <- array(0,dim=c(q,q,Gibbs_size))
Tinv[,,1] <- diag(rep(1, q))


# Hyperparameters
nu <- 2
V <- diag(rep(1, q))

RSS <- rep(0, m) # used for sigmasq update
prod_for_T_update <- rep(0, m) # used for T update
sumXTX <- matrix(0, ncol=p, nrow=p) # used for B update
sumXTymXg <- rep(0, p) # used for B update

# Precache the computation of sumXTX
for(i in 1:m){
  sumXTX <- sumXTX + crossprod(Xi[[i]])
}

for (iter in 2:Gibbs_size){
  if(iter %% 100 == 0){print(iter)}
  #Update sigmasq
  for (i in 1:m){
    temp <- yi[[i]] - as.numeric(Zi[[i]] %*% g[i,,iter-1]) - as.numeric(Xi[[i]] %*% B[,iter-1])
    RSS[i] <- crossprod(temp)
  }
  sigmasq[iter] <- 1/rgamma(1, n/2, sum(RSS)/2)
  #Update Tinv
  prod_for_T_update <- tcrossprod(g[i,,iter-1])
  Tinv[,,iter] <- rWishart(1,nu + m, solve(V + prod_for_T_update))
  # Update the gamma_i's:
  for(i in 1:m){
    Tstar <- solve(Tinv[,,iter] + (1/sigmasq[iter]) * crossprod(Zi[[i]]))
    kappastar <- (1/sigmasq[iter]) * Tstar %*% (crossprod(Zi[[i]], yi[[i]] - Xi[[i]] %*% B[,iter-1]))
    g[i,,iter] <- rmvnorm(1, kappastar, Tstar)
  }
  #Update B:
  for(i in 1:m){
    sumXTymXg <- sumXTymXg + crossprod(Xi[[i]], yi[[i]] - Xi[[i]] %*% g[i,,iter])
  }
  Sigmastar <- solve(diag(rep(1, p)) + (1/sigmasq[iter]) * sumXTX)
  mustar <- Sigmastar %*% ((1/sigmasq[iter])*sumXTymXg)
  B[,iter] <- rmvnorm(1, mustar, Sigmastar)
  # Reset the auxiliary quantity sumXTymXg
  sumXTymXg <- rep(0, p)
}

#Remove burnin and thin the chain
sigmasq <- sigmasq[seq(burnin+1, Gibbs_size, by = thin)]
B <- B[,seq(burnin+1, Gibbs_size, by = thin)]
g <- g[,,seq(burnin+1, Gibbs_size, by = thin)]
Tinv <- Tinv[,,seq(burnin+1, Gibbs_size, by = thin)]


# ======================================
# Analyse the convergence via traceplots
# ======================================

#sigmasq
pdf(file='./04_02_plot_Gibbs_sigmasq.pdf',width=7,height=5)
par(mar=c(4,4,2,2), mfrow = c(1,1), family = 'Palatino',pin=c(5,3))
plot(sigmasq, type = 'l', col="blue", xlab = 'Iterations', ylab = bquote(sigma^2))
dev.off()

#beta
pdf(file='./04_02_plot_Gibbs_beta.pdf',width=7,height=5)
par(mar=c(4,4,2,2), mfrow = c(2,2), family = 'Palatino')
for (i in 1:4)
  plot(B[i,], col="blue",type = 'l', xlab = 'Iterations', ylab = bquote(beta[.(i)]))
dev.off()

# The 88 random effects gamma_1, ..., gamma_m for the intercept
pdf(file='./04_02_plot_Gibbs_gamma.pdf',width=7,height=6)
par(mar=c(2,4,1,1), mfrow = c(2,2), family = 'Palatino')
for (i in 1:4)
  plot(g[i,1,], type = 'l', xlab = 'Iterations', ylab = bquote(gamma[.(1)~.(i)]))
dev.off()


# ===============================
# Analyse the posterior estimates
# ===============================

# Compute the posterior point estimates
betas.post <- rowMeans(B)
gammas.post <- apply(g, c(1,2), mean)

# Draw the demand curves
cols <- c('blue','red')
pdf(file='./04_02_plot_demand_curves.pdf',width=7,height=4)
par(mar=c(2,2,1,1), mfrow = c(2,3), family = 'Palatino')
for (i in c(2,4,11,18,36,68)){
  xgrid <- seq(min(Xi[[i]][,2]), max(Xi[[i]][,2]), length.out = 100)
  plot(exp(Xi[[i]][,2]), exp(yi[[i]]), main = as.character(cheese_data$store[i]),type = 'p', col = cols[Xi[[i]][,3]+1], pch = 16, cex = 0.8)
  lines(exp(xgrid), exp((betas.post[1] + gammas.post[i,1]) + xgrid * (betas.post[2] + gammas.post[i,2])), col = cols[1], lwd = 2)
  lines(exp(xgrid), exp((betas.post[1] + betas.post[3] + gammas.post[i,1]+ gammas.post[i,3]) + xgrid * (betas.post[2] + betas.post[4] + gammas.post[i,2] + gammas.post[i,4])), col = cols[2], lwd = 2)
}
dev.off()


# Draw the fitted lines
pdf(file='./04_02_plot_fitted_lines.pdf',width=7,height=7)
par(mar=c(2,2,1,1), mfrow = c(5,5), family = 'Palatino')
for (i in 1:25){
  xgrid <- seq(min(Xi[[i]][,2]), max(Xi[[i]][,2]), length.out = 100)
  plot(Xi[[i]][,2], yi[[i]], type = 'p', col = cols[Xi[[i]][,3]+1], pch = 16, cex = 0.8)
  lines(xgrid, (betas.post[1] + gammas.post[i,1]) + xgrid * (betas.post[2] + gammas.post[i,2]), col = cols[1], lwd = 2)
  lines(xgrid, (betas.post[1] + betas.post[3] + gammas.post[i,1]+ gammas.post[i,3]) + xgrid * (betas.post[2] + betas.post[4] + gammas.post[i,2] + gammas.post[i,4]), col = cols[2], lwd = 2)              
}
dev.off()

# Fit four stores that are peculiar (few data points for one of the categories)
pdf(file='./04_02_plot_fitted_lines_peculiar.pdf',width=7,height=4)
xgrid <- seq(min(X[,2]), max(X[,2]), length.out = 100)
par(mar=c(2,2,2,2), mfrow = c(2,2), cex = 0.8, family = 'Palatino')
for (i in c(1,9,18,34)){
  plot(Xi[[i]][,2], yi[[i]], type = 'p', main = as.character(cheese_data$store[i]), col = cols[Xi[[i]][,3]+1], pch = 16, cex = 0.8, xlim = c(min(X[,2]), max(X[,2])), ylim = c(min(log_cheese_data$vol), max(log_cheese_data$vol)))
  lines(xgrid, (betas.post[1] + gammas.post[i,1]) + xgrid * (betas.post[2] + gammas.post[i,2]), col = cols[1], lwd = 2)
  lines(xgrid, (betas.post[1] + betas.post[3] + gammas.post[i,1]+ gammas.post[i,3]) + xgrid * (betas.post[2] + betas.post[4] + gammas.post[i,2] + gammas.post[i,4]), col = cols[2], lwd = 2)              
}
dev.off()


# Same four stores via OLS
xgrid <- seq(min(X[,2]), max(X[,2]), length.out = 100)
pdf(file='./04_02_plot_fitted_lines_peculiar_OLS.pdf',width=7,height=4)
par(mar=c(2,2,2,2), mfrow = c(2,2), cex = 0.8, family = 'Palatino')
for (i in c(1,9,18,34)){
  plot(Xi[[i]][,2], yi[[i]], type = 'p', main = as.character(cheese_data$store[i]), col = cols[Xi[[i]][,3]+1], pch = 16, cex = 0.8, xlim = c(min(X[,2]), max(X[,2])), ylim = c(min(log_cheese_data$vol), max(log_cheese_data$vol)))
  fit = lm(yi[[i]] ~ -1 + Xi[[i]])
  betas.lm <- fit$coefficients
  lines(xgrid, (betas.lm[1]) + xgrid * (betas.lm[2]), col = cols[1], lwd = 2)
  lines(xgrid, (betas.lm[1]+ betas.lm[3]) + xgrid * (betas.lm[2] + betas.lm[4]), col = cols[2], lwd = 2)              
}
dev.off()
# The shrinkage effect is evident for the stores whose OLS estimates were distorted.



