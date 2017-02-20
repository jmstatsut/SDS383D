ComputeTestErrors <- function(h,vecx,vecy) {
  h <- as.vector(h)
  #center the data
  vecx <- vecx - sum(vecx)/length(vecx)
  vecy <- vecy - sum(vecy)/length(vecy)
  #split x and y into training and testing sets
  foldassignments = cut(sample(1:length(vecx), length(vecx)), breaks = 3, labels = FALSE)
  x_test <- vecx[foldassignments == 1]
  x_train <- vecx[foldassignments != 1]
  y_test <- vecy[foldassignments == 1]
  y_train <- vecy[foldassignments != 1]
  pred_error <- vector(mode="numeric",length=length(h))
  for (i in 1:length(h)){
    y_pred <- KernelSmoother(x_test,x_train,y_train,K,h[i])
    pred_error[i] <- sum((y_pred-y_test)^2)/length(y_pred)
  }
  return(as.vector(pred_error))
}

#create data
# vecx <- 1:60/5
# vecy <- sin(vecx) + rnorm(length(vecx),0,0.4^2)

# vecerrors <- ComputeTestErrors(c(.1,.5,1,2,5,10,20,30))