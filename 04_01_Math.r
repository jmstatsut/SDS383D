setwd("C:/Users/Jesse/Box Sync/!Jesse/Stats/Courses/2017_1/SDS 383D Statistical Modeling II/Exercise Solutions/04/")
test_data <- read.csv("04_01_mathtest.csv",header=TRUE)


grand_mean <- mean(test_data$mathscore)

M <- nrow(test_data)#total number of observations
n <- length(unique(test_data[,1]))#number of groups (that is, schools)

AggData <- matrix(0,nrow=n,ncol=2)

AggData[,1] <- tapply(test_data$mathscore,test_data$school,length)
AggData[,2] <- tapply(test_data$mathscore,test_data$school,mean)

plot(AggData,xlab="sample size",ylab="school sample mean",col="blue")
abline(h=grand_mean,col="green",lty=2)
# plot(AggData[,1],abs(AggData[,2]-grand_mean),xlab="sample size",ylab="mean deviation from grand mean",col="blue")

iterations <- 10000 #how many iterations for Gibbs

#create structures for (hyper)parameters and initialize values
sigmasq <- vector(length=iterations+1)
tausq <- vector(length=iterations+1)
mu <- vector(length=iterations+1)
theta <- matrix(nrow=n,ncol=iterations+1)#row i is Gibbs sample for theta_i
sigmasq[1] <- 1
tausq[1] <- 1
mu[1] <- grand_mean
theta[,1] <- AggData[,2]

for(it in 1:iterations){
  thetaSD <- sigmasq[it]*tausq[it]/(n*tausq[it]+1)
  thetaMean <- (n*AggData[,2]/sigmasq[it]+mu[it]/(sigmasq[it]*tausq[it]))*thetaSD
  theta[,it+1] <- rnorm(n,thetaMean,thetaSD)
  mu[it+1] <- rnorm(1,mean(theta[,it+1]),sigmasq[it]*tausq[it]/n)
  SSgrouped <- sum((test_data-theta[test_data[,1],it+1])^2)
  SStheta <- sum((theta[,it+1]-mu[it+1])^2)
  sigmasq[it+1] <- 1/rgamma(1,(n+M)/2,0.5*(SSgrouped+SStheta/tausq[it]))
  tausq[it+1] <- 1/rgamma(1,n/2,SStheta/(2*sigmasq[it+1]))
}

kappa <- (AggData[,2]-rowMeans(theta))/AggData[,2]

plot(AggData[,1],abs(kappa),xlab="sample size",ylab="shrinkage factor",col="blue")

boxplot(mathscore~school,data=test_data)
points(rowMeans(theta),col="green")
