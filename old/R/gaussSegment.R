setClass("gaussKnown",
         slots = c(n = "numeric", cost = "numeric", beta = "numeric",
                   param = "numeric", summaryStats = "numeric"0.0),
         prototype = list( n = 0, cost = 0, beta = 0,
                          param = c(0,1), summaryStats = rep(0,4))
         )


setMethod("update","gaussKnown",
          function(obj,x,mu,sigma){
              obj@n <- obj@n + 1.0
              obj@summaryStats <- obj@summaryStats + c(1.0/sigma,log(sigma),(x-mu)/sigma; ((x-mu)^2.0)/sigma)
              kappa = (obj@summaryStats[4] - 2.0*obj@param[1]*obj@summaryStats[3] +
                       (obj@param[1]^2.0)*obj@summaryStats[1] ) / n
              obj@cost = obj@n * log( 2.0*pi*obj@param[2]) + obj@summaryStats[2] +
                  (obj@n*kappa/obj@param[2]) + obj@beta
              return(obj)
          })


#' create a gaussKnown segment
#' @param beta penalty term
#' @param start time step at start
#' 
#' @export
gaussKnown <- function(beta,start){ new("partition",beta=beta,start=start,cost=beta) }

## ##############################
## gauss Unknown mean
setClass("gaussMean",
         slots = c(n = "numeric", cost = "numeric", beta = "numeric",
                   param = "numeric", summaryStats = "numeric"0.0),
         prototype = list( n = 0, cost = 0, beta = 0,
                          param = c(0,1), summaryStats = rep(0,4))
         )


setMethod("update","gaussMean",
          function(obj,x,mu,sigma){
              obj@n <- obj@n + 1.0
              obj@summaryStats <- obj@summaryStats + c(1.0/sigma,log(sigma),(x-mu)/sigma; ((x-mu)^2.0)/sigma)
              obj@param[1] <- obj@summaryStats[3] / obj@summaryStats[1]
              kappa = (obj@summaryStats[4] - 2.0*obj@param[1]*obj@summaryStats[3] +
                       (obj@param[1]^2.0)*obj@summaryStats[1] ) / n
              obj@cost = obj@n * log( 2.0*pi*obj@param[2]) + obj@summaryStats[2] +
                  (obj@n*kappa/obj@param[2]) + obj@beta;
              return(obj)
          })


#' create a gaussMean segment
#' @param beta penalty term
#' @param start time step at start
#' 
#' @export
gaussMean <- function(beta,start){ new("partition",beta=beta,start=start,cost=beta) }






## ##Rcpp::loadModule("gaussSeg",TRUE)

## #' meanSegment
## #' @param beta Penalisation term for the segment
## #' @param t start time of segment
## #' 
## #' @export
## gaussMean <- function(beta,t){
##     out <- new(gaussMeanCpp, start=as.integer(t), penalty=beta)
##     return( out )
## }


## ####################################################################################
## gauss variance
#' gaussVar
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussVar <- function(beta,t){
    out <- new(gaussVarCpp, start=as.integer(t), penalty=beta)
    return( out )
}


## ####################################################################################
## gauss mean & variance
#' gaussMeanVar
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussMeanVar <- function(beta,t){
    out <- new(gaussMeanVarCpp, start=as.integer(t), penalty=beta)
    return( out )
}

## ####################################################################################
## gauss fixed
#' meanSegment
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussFixed <- function(beta,t){
    out <- new(gaussFixedCpp, start=as.integer(t), penalty=beta)
    return( out )
}


## ###################################################################################
#' gaussPoint
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussPoint <- function(beta,t){
    out <- new(gaussPointCpp, start=as.integer(t), penalty=beta)
    return( out )
}
