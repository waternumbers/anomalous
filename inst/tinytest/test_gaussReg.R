## test the pelt and op algorithms with example data
## test results from the changepoint package unless stated
## computed using th BIC not the MBIC
## changepoint returns the last timestep before a change
##rm(list=ls())
##devtools::load_all()

## simple test
set.seed(10)
x <- list()
n <- 120
for(ii in 1:48){
    
    if(ii < 10){ theta = c(1,0); sigma <- 0.1 }
    if(ii >= 10 & ii <12){ theta <- c(10,0); sigma <- 2}
    if(ii >= 12 & ii < 44){ theta <- c(5,1); sigma <- 2}
    if(ii >= 44 ){ theta <- c(1,0); sigma <- 0.1}
    
    X <- cbind(rep(1,n),runif(n,ii-1,ii))
    y <- rnorm(n, X%*%theta, sigma)
    x[[ii]] <- list(y=y,X=X)
}

expect_silent({
    fCost <- gaussRegMeanVar$new(x)
    p <- partition(4*log(sum(sapply(y,length))),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(9,11,43,48) )


## test with singular design matrix

xx <- x
for(ii in c(1:10,44:48)){
    xx[[ii]]$X[,2] <- 0
}
expect_silent({
    fCost <- gaussRegMeanVar$new(xx)
    p <- partition(4*log(sum(sapply(y,length))),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(9,11,43,48) )
