## test the pelt and op algorithms with example data
## test results from the changepoint package unless stated
## computed using th BIC not the MBIC
## changepoint returns the last timestep before a change
## this function alters the anomolous output to match this indexing
fidx <- function(tmp){
    out <- sapply(tmp@collective,function(x){x@start})
    out <- out[-1] - 1
    return(out) 
}
set.seed(10)
x <- c(rnorm(100, 0, 1), rnorm(100, 1, 1), rnorm(100, 0, 1), rnorm(100, 0.2, 1))
pen <- 2*log(length(x))
mu <- rep(0,length(x))
sigma <- rep(1,length(x))
expect_equal( fidx(op(x,mu,sigma,gaussMean,pen)), c(97,192) )
expect_equal( fidx(pelt(x,mu,sigma,gaussMean,pen)), c(97,192) )

pen <- 1.5*log(length(x))
expect_equal( fidx(op(x,mu,sigma,gaussMean,pen)), c(97,192,273) )
expect_equal( fidx(pelt(x,mu,sigma,gaussMean,pen)), c(97,192,273) )

## check change in mu has no impact
pen <- 2*log(length(x))
x[50:150] <- x[50:150] + 5
mu[50:150] <- mu[50:150] + 5
expect_equal( fidx(pelt(x,mu,sigma,gaussMean,pen)), c(97,192) )


data("Lai2005fig4")#, package = "changepoint")
x <- Lai2005fig4[, 5]
pen <- 2*log(length(x))
mu <- rep(0,length(x))
sigma <- rep(1,length(x))
expect_equal( fidx(op(x,mu,sigma,gaussMean,pen)), c(81,85,89,96,123,133) )
expect_equal( fidx(pelt(x,mu,sigma,gaussMean,pen)), c(81,85,89,96,123,133) )


data("wind") #, package = "gstat")
x <- diff(wind[, 11])
mu <- rep(0,length(x))
sigma <- rep(1,length(x))
pen <- 2*log(length(x))
##expect_equal( fidx(pelt(x,mu,sigma,gaussVar, pen)),
##             c(3409, 3496, 5054, 5184, 5203, 5373, 5583, 5678, 5728, 6235, 6241, 6542)  ) # currently takes to long to run


data("discoveries", package = "datasets")
x <- discoveries
pen <- 2*log(length(x))
mu <- rep(NA,length(x))
sigma <- rep(NA,length(x))
expect_equal( fidx(op(x,mu,sigma,poisRate,pen)), c(24, 29,73) )
expect_equal( fidx(pelt(x,mu,sigma,poisRate,pen)), c(24, 29,73) )

