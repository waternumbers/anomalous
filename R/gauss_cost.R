#' @export
gaussCost <- R6Class("gaussCost",
                     public=list(
                         summaryStats=NULL,
                         initialize = function(x,m=0,s=1){
                             self$summaryStats <- apply(cbind(1/s,log(s),(x-m)/s,((x-m)^2)/s),2,cumsum)
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             n*log(2*pi) + sumStat[2] + sumStat[4] + pen
                         },
                         pointCost = function(a,pen){
                             if(a<2){
                                 sumStat <- self$summaryStats[a,]
                             }else{
                                 sumStat <- self$summaryStats[a,] - self$summaryStats[a-1,]
                             }
                             gamma <- exp(-(1+pen))
                             log(2*pi) + sumStat[2] + log(gamma + sumStat[4]) + 1 + pen
                         }
                     ),
                     private=list(
                         meanChange = function(a,b,pen){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             mhat <- sumStat[3] / sumStat[1]
                             n*log(2*pi) + sumStat[2] + sumStat[4] - (mhat^2)*sumStat[1] + pen
                         },
                         varChange = function(a,b,pen){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             shat <- sumStat[4] / n ## TODO add catch for close to zero
                             shat <- max(shat,.Machine$double.xmin)
                             n*log(2*pi*shat) + sumStat[2] + n + pen
                         },
                         meanVarChange = function(a,b,pen){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             mhat <- sumStat[3] / sumStat[1]
                             shat <- (sumStat[4] - (mhat^2)*sumStat[1])/n
                             shat <- max(shat,.Machine$double.xmin)
                             n*log(2*pi*shat) + sumStat[2] + n + pen
                         }
                     )
                     )
                     
#' @export
gaussMean <- R6Class("gaussMean",
                     inherit = gaussCost,
                     public = list(
                         collectiveCost = function(a,b,pen){ private$meanChange(a,b,pen) }
                     )
                     )

#' @export
gaussVar <- R6Class("gaussVar",
                    inherit = gaussCost,
                    public = list(
                        collectiveCost = function(a,b,pen){ private$varChange(a,b,pen) }
                    )
                    )

#' @export
gaussMeanVar <- R6Class("gaussMeanVar",
                        inherit = gaussCost,
                        public = list(
                            collectiveCost = function(a,b,pen){ private$meanVarChange(a,b,pen) }
                        )
                        )





## ## set up the cost class and methods
## setGeneric("collectiveCost",function(obj,...){standardGeneric("collectiveCost")})
## setGeneric("baseCost",function(obj,...){standardGeneric("baseCost")})
## setGeneric("pointCost",function(obj,...){standardGeneric("pointCost")})

## #' Gaussian Cost Functions for a univariate series
## #'
## #' An anomaly is a change in mean
## #' @slot summaryStats a record of the summary statistics required, one row per time step
## setClass("gaussMean", representation(summaryStats = "matrix"))

## #' @export
## gaussMean <- function(x,m=0,s=1){
##     new("gaussMean",summaryStats = apply(cbind(1/s,log(s),(x-m)/s,((x-m)^2)/s),2,cumsum))
## }

## setMethod("collectiveCost","gaussMean",function(obj,a,b,pen){
##     a <- a-1
##     n <- b-a
##     if(a<1){
##         sumStat <- obj@summaryStats[b,]
##     }else{
##         sumStat <- obj@summaryStats[b,] - obj@summaryStats[a,]
##     }
##     mhat <- sumStat[3] / sumStat[1]
##     n*log(2*pi) + sumStat[2] + sumStat[4] - (mhat^2)*sumStat[1] + pen
## })

## setMethod("baseCost","gaussMean",function(obj,a,b,pen=0){
##     a <- a-1
##     n <- b-a
##     if(a<1){
##         sumStat <- obj@summaryStats[b,]
##     }else{
##         sumStat <- obj@summaryStats[b,] - obj@summaryStats[a,]
##     }
##     n*log(2*pi) + sumStat[2] + sumStat[4] + pen
## })

## setMethod("pointCost","gaussMean",function(obj,a,pen){
##     if(a<2){
##         sumStat <- obj@summaryStats[a,]
##     }else{
##         sumStat <- obj@summaryStats[a,] - obj@summaryStats[a-1,]
##     }
##     gamma <- exp(-(1+pen))
##     log(2*pi) + sumStat[2] + log(gamma + sumStat[4]) + 1 + pen
## })
