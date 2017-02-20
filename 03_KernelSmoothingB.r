#x is the new input value
#vecx and vecy are the data vectors, assumed to be centered
#k is the kernel function
#h is the bandwidth
KernelSmoother <- function(x,vecx,vecy,k,h) {
  x <- as.vector(x)
  estvalues <- vector(mode="numeric",length=length(x))
  for (i in 1:length(x)){
    #divide by sum to normalize
    estvalues[i] <- sum((1/h)*k((x[i]-vecx)/h)*vecy)/sum((1/h)*k((x[i]-vecx)/h))
  }
  return(as.vector(estvalues))
}

K <- function(x){
  exp(-x^2/2)/sqrt(2*pi)
}
# h <- 0.5
# vecx <- 1:40/5
# vecy <- sin(vecx) + rnorm(length(vecx),0,0.4^2)
#center the data
# vecx <- vecx - sum(vecx)/length(vecx)
# vecy <- vecy - sum(vecy)/length(vecy)

#vecz contains the predicted values on vecx+0.1
# vecz <- vector(length=length(vecx))
# vecz <- KernelSmoother(vecx+.1,vecx,vecy,K,h)
# plot(vecx,vecz)