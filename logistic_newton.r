
log.regression <- function(input.data)
{
	input.data <- cbind( c(rep(1,nrow(input.data)))  ,input.data)
	iterations <- 0
	iter.max <- 1000
	epsilon <- .Machine$double.eps
	theta <- as.matrix(c(rep(0,ncol(input.data) - 1)))
	J_vect <- NULL #for analysis
	repeat
	{
		theta_old <- theta
		H <- matrix(0, (ncol(input.data) - 1), (ncol(input.data) - 1) ) 
		dell_theta <- matrix(0, (ncol(input.data) - 1)  ,1)
		m <- nrow(input.data)
		J <- 0		#for analysis
		for (i in 1:m)
		{
			y <- input.data[i,ncol(input.data)]
			h <- h_logistic.reg( theta, input.data[i,1:(ncol(input.data)-1)] )
			x <- as.matrix( input.data[i,1:(ncol(input.data)-1)]  )
			dell_theta <- dell_theta + ( ( h - y ) * x ) * ( 1 / m )
			H <- H + ( h ) * ( 1 - h ) * (1/ m )  * ( x %*% t(x) )
			J <- J + ( (-1/m) * (y*log(h) + (1-y)*log(1-h)) )	#for analysis
		}
		J_vect <- c(J_vect,J)	#for analysis
		theta <- theta - (solve(H) %*% dell_theta)
		iterations <- iterations + 1
		cat('iter. no.',iterations,' : ',theta,'\n')
		if ( max(abs(theta - theta_old)) < epsilon || iterations > iter.max )
		{
			plot(J_vect,type = 'b', pch = 20) 	#for analysis
			return(theta)
		}
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
