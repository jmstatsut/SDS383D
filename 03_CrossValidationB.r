#create data
x <- 1:500/500
y_wiggleln <- sin(100*x) + rnorm(length(x),0,0.4^2)
y_wigglehn <- sin(100*x) + rnorm(length(x),0,0.7^2)
y_smoothln <- x*(1-x) + rnorm(length(x),0,0.1^2)
y_smoothhn <- x*(1-x) + rnorm(length(x),0,0.3^2)

h <- seq(from=0.01,to=5,by=0.01)

vecerrors <- ComputeTestErrors(h,x,y_wiggleln)