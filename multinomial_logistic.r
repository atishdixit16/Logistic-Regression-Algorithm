#input <- scan('Skin_NonSkin.txt',fill=TRUE)
#input.data <- matrix(input,ncol = 4, byrow = T)

linear.regression <- function(input.data,alpha)
{
	input.data <- cbind( c(rep(1,nrow(input.data)))  ,input.data)
	iterations <- 0
	iter.max <- 1000
	epsilon <- .Machine$double.eps
	outcomes <- length(unique(input.data[,ncol(input.data)]))
	theta <- matrix(0,nrow(input.data),outcomes)
	repeat
	{
		theta_old <- theta
		for (i in 1:outcomes )
		{
			sum <- 0
			for (j in 1:nrow(input.data))
				sum <- sum + ( input.data[j,ncol(input.data)] - h_logistic.reg(theta,input.data[j,1:(ncol(input.data) - 1)]) ) * input.data[j,i]
			theta[i] <- theta[i] + alpha*sum
		}
		iterations <- iterations + 1
		cat('iter. no.',iterations,' : ',theta,'\n')
		if ( max(abs(theta - theta_old)) < epsilon || iterations > iter.max )
			return(theta)
	}
}

h_lin.reg <- function(theta,x)
{
	return(sum(theta*x))
}

h_logistic.reg <- function(theta,x)
{
	return(  1/( 1+exp(-(h_lin.reg(theta,x))) )  )
}
