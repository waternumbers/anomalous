## test the pelt and op algorithms with example data
## test results from the changepoint package unless stated
set.seed(10)
x <- c(rnorm(100, 0, 1), rnorm(100, 1, 1), rnorm(100, 0, 1), rnorm(100, 0.2, 1))
pen <- 2*log(length(x))
expect_equal(op(x,fCmn,pen), c(97,192))
expect_equal(pelt(x,fCmn,pen), c(97,192))

pen <- 1.5*log(length(x))
expect_equal(op(x,fCmn,pen), c(97,192,273))
expect_equal(pelt(x,fCmn,pen), c(97,192,273))

data("Lai2005fig4", package = "changepoint")
x <- Lai2005fig4[, 5]
pen <- 2*log(length(x))
expect_equal(op(x,fCmn,pen), c(81,85,89,96,123,133))
expect_equal(pelt(x,fCmn,pen), c(81,85,89,96,123,133))

data("wind", package = "gstat")
x <- diff(wind[, 11])
#wind.pelt <- cpt.var(diff(wind[, 11]), method = "PELT")
#plot(wind.pelt, xlab = "Index")
#logLik(wind.pelt)
pelt(x,fCv, b=2*log(length(x)))

data("discoveries", package = "datasets")
x <- discoveries
pen <- 2*log(length(x))
expect_equal(op(x,fCpois,pen), c(24, 29,73))
expect_equal(pelt(x,fCpois,pen), c(24, 29,73))

# Placeholder with simple test
expect_equal(1 + 1, 2)

