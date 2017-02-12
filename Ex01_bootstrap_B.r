simulate_mvn <- function(m,v,n) {
  d <- nrow(v)
  cd <- t(chol(v))
  mvn_samp <- matrix(nrow=n,ncol=d)#store samples as rows
  for(i in 1:n){
    std_norm_samp <- rnorm(d)
    mvn_samp[i,] <- cd %*% std_norm_samp + m
  }
  return(mvn_samp)
}

m <- c(4,5)
v <- matrix(c(1,2,2,3),nrow=2,ncol=2)#make this symmetric
v <- t(v) %*% v#make sure the matrix is positive definite

N <- 10000

mvn_samp <- simulate_mvn(m,v,N)

mu_hat <- colSums(mvn_samp)/N

sigma_hat <- t(mvn_samp - mu_hat) %*% (mvn_samp - mu_hat)/N

#now use bootstrap to estimate sample distribution variance
B <- 10000
n <- nrow(mvn_samp)
p <- ncol(mvn_samp)
mvn_samp_boot <- matrix(nrow=B,ncol=p)
sigma_hat_boot <- matrix(0,nrow=p,ncol=p)
for(i in 1:B) {
  boot_sample_indexes <- sample(1:n,n,replace=TRUE)
  mvn_samp_temp <- mvn_samp[boot_sample_indexes,]
  mu_hat_temp <- colSums(mvn_samp_temp)/N
  sigma_hat_boot <- sigma_hat_boot + (t(mvn_samp_temp - mu_hat_temp) %*% (mvn_samp_temp - mu_hat_temp)/N)
}
sigma_hat_boot <- sigma_hat_boot/B