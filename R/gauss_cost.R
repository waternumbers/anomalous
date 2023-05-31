## set up the cost class and methods
setGeneric("collectiveCost",function(obj,...){standardGeneric("collectiveCost")})
setGeneric("baseCost",function(obj,...){standardGeneric("baseCost")})
setGeneric("pointCost",function(obj,...){standardGeneric("pointCost")})

#' Gaussian Cost Functions for a univariate series
#'
#' An anomaly is a change in mean
#' @slot summaryStats a record of the summary statistics required, one row per time step
setClass("gaussMean", representation(summaryStats = "matrix"))

#' @export
gaussMean <- function(x,m=0,s=1){
    new("gaussMean",summaryStats = cbind(cumsum((x-m)/s),cumsum(x^2)))
}

setMethod("collectiveCost","gaussMean",function(obj,a,b,pen){
    a <- a-1
    n <- b-a
    if(a<1){
        Syy <- obj@summaryStats[b,2]
        Sy <- obj@summaryStats[b,1]
    }else{
        Syy <- obj@summaryStats[b,2] - obj@summaryStats[a,2]
        Sy <- obj@summaryStats[b,1] - obj@summaryStats[a,1]
    }
    ( n*log(2*pi) + Syy - (Sy^2)/n ) + pen
})

setMethod("baseCost","gaussMean",function(obj,a,b,pen=0){
    a <- a-1
    n <- b-a
    if(a<1){
        Syy <- obj@summaryStats[b,2]
    }else{
        Syy <- obj@summaryStats[b,2] - obj@summaryStats[a,2]
    }
    ( n*log(2*pi) + Syy ) + pen
})

setMethod("pointCost","gaussMean",function(obj,a,pen){
    if(a<2){
        Syy <- obj@summaryStats[a,2]
    }else{
        Syy <- obj@summaryStats[a,2] - obj@summaryStats[a-1,2]
    }
    gamma <- exp(-(1+pen))
    log(2*pi) + log(gamma + Syy) + 1 + pen
})
