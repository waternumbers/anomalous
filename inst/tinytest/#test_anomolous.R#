## test the pelt and op algorithms with example data
## test results from the changepoint package unless stated
## computed using th BIC not the MBIC
set.seed(10)
x <- c(rnorm(100, 0, 1), rnorm(100, 1, 1), rnorm(100, 0, 1), rnorm(100, 0.2, 1))
pen <- 2*log(length(x))
expect_equal(op(x,fCmn,pen), c(97,192))
expect_equal(pelt(x,fCmn,pen), c(97,192))

pen <- 1.5*log(length(x))
expect_equal(op(x,fCmn,pen), c(97,192,273))
expect_equal(pelt(x,fCmn,pen), c(97,192,273))

data("Lai2005fig4")#, package = "changepoint")
x <- Lai2005fig4[, 5]
pen <- 2*log(length(x))
expect_equal(op(x,fCmn,pen), c(81,85,89,96,123,133))
expect_equal(pelt(x,fCmn,pen), c(81,85,89,96,123,133))

data("wind") #, package = "gstat")
x <- diff(wind[, 11])
##changepoint::cpt.var(x-mean(x), penalty="BIC",method = "PELT")
##plot(wind.pelt, xlab = "Index")
##logLik(wind.pelt)
##pelt(x-mean(x),fCv, Beta=2*log(length(x)))
expect_equal( pelt(x,fCv, Beta=2*log(length(x))),
             c(3409, 3496, 5054, 5184, 5203, 5373, 5583, 5678, 5728, 6235, 6241, 6542) )


data("discoveries", package = "datasets")
x <- discoveries
pen <- 2*log(length(x))
expect_equal(op(x,fCpois,pen), c(24, 29,73))
expect_equal(pelt(x,fCpois,pen), c(24, 29,73))

