#' @include generics.R
setClass("gaussKnown",
         slots = c(n = "numeric", cost = "numeric", beta = "numeric", start="integer",
                   param = "numeric", summaryStats = "numeric"),
         prototype = list( n = 0, cost = 0, beta = 0,
                          param = c(0,1), summaryStats = rep(0,4))
         )


setMethod("update","gaussKnown",
          function(obj,x,mu,sigma){
              obj@n <- obj@n + 1.0
              obj@summaryStats <- obj@summaryStats + c(1.0/sigma, log(sigma), (x-mu)/sigma, ((x-mu)^2.0)/sigma)
              kappa = (obj@summaryStats[4] - 2.0*obj@param[1]*obj@summaryStats[3] +
                       (obj@param[1]^2.0)*obj@summaryStats[1] ) / obj@n
              obj@cost = obj@n * log( 2.0*pi*obj@param[2]) + obj@summaryStats[2] +
                  (obj@n*kappa/obj@param[2]) + obj@beta
              return(obj)
          })


#' create a gaussKnown segment
#' @param beta penalty term
#' @param start time step at start
#' 
#' @export
gaussKnown <- function(beta,start){
    beta <- as.numeric(beta)
    start <- as.integer(start)
    new("gaussKnown",beta=beta,start=start,cost=beta)
}

## ##############################
## gauss Unknown mean
setClass("gaussMean",
         slots = c(n = "numeric", cost = "numeric", beta = "numeric", start="integer",
                   param = "numeric", summaryStats = "numeric"),
         prototype = list( n = 0, cost = 0, beta = 0,
                          param = c(0,1), summaryStats = rep(0,4))
         )


setMethod("update","gaussMean",
          function(obj,x,mu,sigma){
              obj@n <- obj@n + 1.0
              obj@summaryStats <- obj@summaryStats + c(1.0/sigma, log(sigma), (x-mu)/sigma, ((x-mu)^2.0)/sigma)
              obj@param[1] <- obj@summaryStats[3] / obj@summaryStats[1]
              kappa = (obj@summaryStats[4] - 2.0*obj@param[1]*obj@summaryStats[3] +
                       (obj@param[1]^2.0)*obj@summaryStats[1] ) / obj@n
              obj@cost = obj@n * log( 2.0*pi*obj@param[2]) + obj@summaryStats[2] +
                  (obj@n*kappa/obj@param[2]) + obj@beta;
              return(obj)
          })


#' create a gaussMean segment
#' @param beta penalty term
#' @param start time step at start
#' 
#' @export
gaussMean <- function(beta,start){
    beta <- as.numeric(beta)
    start <- as.integer(start)
    new("gaussMean",beta=beta,start=start,cost=beta)
}

## ##############################
## gauss Unknown mean
setClass("gaussVar",
         slots = c(n = "numeric", cost = "numeric", beta = "numeric", start="integer",
                   param = "numeric", summaryStats = "numeric"),
         prototype = list( n = 0, cost = 0, beta = 0,
                          param = c(0,1), summaryStats = rep(0,4))
         )

setMethod("update","gaussVar",
          function(obj,x,mu,sigma){
              obj@n <- obj@n + 1.0
              obj@summaryStats <- obj@summaryStats + c(1.0/sigma, log(sigma), (x-mu)/sigma, ((x-mu)^2.0)/sigma)
              kappa = (obj@summaryStats[4] - 2.0*obj@param[1]*obj@summaryStats[3] +
                       (obj@param[1]^2.0)*obj@summaryStats[1] ) / obj@n
              kappa <- max(kappa,.Machine$double.xmin)
              obj@param[2] <- kappa
              obj@cost = obj@n * log( 2.0*pi*obj@param[2]) + obj@summaryStats[2] +
                  (obj@n*kappa/obj@param[2]) + obj@beta;
              return(obj)
          })


#' create a gaussVar segment
#' @param beta penalty term
#' @param start time step at start
#' 
#' @export
gaussVar <- function(beta,start){
    beta <- as.numeric(beta)
    start <- as.integer(start)
    new("gaussVar",beta=beta,start=start,cost=beta)
}

## ##############################
## gauss Unknown mean
setClass("gaussMeanVar",
         slots = c(n = "numeric", cost = "numeric", beta = "numeric", start="integer",
                   param = "numeric", summaryStats = "numeric"),
         prototype = list( n = 0, cost = 0, beta = 0,
                          param = c(0,1), summaryStats = rep(0,4))
         )


setMethod("update","gaussMeanVar",
          function(obj,x,mu,sigma){
              obj@n <- obj@n + 1.0
              obj@summaryStats <- obj@summaryStats + c(1.0/sigma, log(sigma), (x-mu)/sigma, ((x-mu)^2.0)/sigma)
              obj@param[1] <- obj@summaryStats[3] / obj@summaryStats[1]
              kappa = (obj@summaryStats[4] - 2.0*obj@param[1]*obj@summaryStats[3] +
                       (obj@param[1]^2.0)*obj@summaryStats[1] ) / obj@n
              kappa <- max(kappa,.Machine$double.xmin)
              obj@param[2] <- kappa           
              obj@cost = obj@n * log( 2.0*pi*obj@param[2]) + obj@summaryStats[2] +
                  (obj@n*kappa/obj@param[2]) + obj@beta;
              return(obj)
          })


#' create a gaussMeanVar segment
#' @param beta penalty term
#' @param start time step at start
#' 
#' @export
gaussMeanVar <- function(beta,start){
    beta <- as.numeric(beta)
    start <- as.integer(start)
    new("gaussMeanVar",beta=beta,start=start,cost=beta)
}









## ##############################
## gauss Point
setClass("gaussPoint",
         slots = c(n = "numeric", cost = "numeric", beta = "numeric", start="integer",
                   param = "numeric", summaryStats = "numeric"),
         prototype = list( n = 0, cost = 0, beta = 0,
                          param = c(0,1), summaryStats = rep(0,4))
         )


setMethod("update","gaussPoint",
          function(obj,x,mu,sigma){
              obj@n <- obj@n + 1.0
              obj@summaryStats <- obj@summaryStats +
                  c(1.0/sigma, log(sigma), (x-mu)/sigma, ((x-mu)^2.0)/sigma)
              kappa = (obj@summaryStats[4] - 2.0*obj@param[1]*obj@summaryStats[3] +
                       (obj@param[1]^2.0)*obj@summaryStats[1] ) / obj@n
              kappa <- max(kappa,2*.Machine$double.xmin)
              obj@param[2] <- kappa

              obj@cost <- log(2*pi) + log(sigma) + log( exp(-obj@beta) + kappa ) + 1 + obj@beta
##              obj@cost = obj@n * log( 2.0*pi*obj@param[2]) + obj@summaryStats[2] +
##                  (obj@n*kappa/obj@param[2]) + obj@beta;
              return(obj)
          })


#' create a gaussPoint segment
#' @param beta penalty term
#' @param start time step at start
#' 
#' @export
gaussPoint <- function(beta,start){
    beta <- as.numeric(beta)
    start <- as.integer(start)
    new("gaussPoint",beta=beta,start=start,cost=beta)
}


