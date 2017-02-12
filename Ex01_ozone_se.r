# Load the library
# you might have to install this the first time
library(mlbench)

# Load the data
ozone = data(Ozone, package='mlbench')

# Look at the help file for details
?Ozone

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

# compute the estimator
betahat = solve(t(x) %*% x) %*% t(x) %*% y
yhat <- x %*% betahat

sigmahat <- (1/(nrow(x)-ncol(x)-1))*sum((y-yhat)^2)
# Fill in the blank
betacov = sigmahat*solve(t(x) %*% x)

sigmahat1 <- (1/(nrow(x)-ncol(x)))*sum((y-yhat)^2)
# Fill in the blank
betacov1 = sigmahat1*solve(t(x) %*% x)

# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as an extra column)
lm1 = lm(y~x-1)

summary(lm1)
betacovlm = vcov(lm1)
sqrt(diag(betacovlm))

