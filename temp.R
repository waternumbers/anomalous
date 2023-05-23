rm(list=ls())
devtools::load_all()

##p <- capaR(rnorm(100),rep(0,100),rep(1,100),gaussKnown,gaussMean,gaussKnown,3*log(100),3*log(100))

## Adapt data to contain a clear point anomaly
data(simulated)

X <- sim.data
X[100,1] <- 45

## adapt X to match the robust scaling applied in v4.0.2
X <- apply(X,2,function(x){ (x-median(x))/mad(x) })

## read in the results from v4.0.2
out <- readRDS("./inst/tinytest/capa_results402_v2.rds") ## TODO need changing to remove path

## #################################
## univariate tests
x <- X[,1]
mu <- rep(0,length(x))
sigma <- rep(1,length(x))
beta <- 4*log(length(x))
betaP <- 3*log(length(x))

tmp <- capa(x,mu,sigma,gaussKnown,gaussMean,gaussPoint,beta,betaP)

collective_anomalies(tmp,"gaussMean")
##anm <- t(sapply(tmp@collective,function(x){c(x@start, x@start+x@n-1, x@n)}))
##cls <- sapply(tmp@collective,class)
##anm[cls=="gaussMean",]
out$single_mean$collective[,c("start","end")]
out$single_mean$point
