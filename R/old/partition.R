## generics

## TODO - change data inputs to matrix so can do multivariate distributions :-)
## TODO - try insertion into all segments to get a classifier?

## basic segment class
#' @include generics.R
setClass("partition", 
         slots = c(collective = "list",
                   point = "list",
                   last_n = "integer",
                   cost = "numeric",
                   isAnom = "logical"),
         prototype = list(
             collective = list(),
             point = list(),
             last_n = integer(1), #NA_integer_,
             cost = 0,#NA_real_,
             isAnom = FALSE
         )
         )



#' Cost of the partioning
#' @param obj partition object
setMethod("cost","partition",function(obj){return(obj@cost)})

#' Add a new Point anomaly
#' @param obj partition object
#' @param mu mean of background distribution
#' @param sigma variance of background distribution
#' @param f function to generate a new class segment
#' @param f function to generate a new class segment
#' @param b penalty for introducing new segment
#' @param t time at start of segment
setMethod("addPoint","partition",
          function(obj,x,mu,sigma,f,b,t){
              tmp <- update(f(b,t),x,mu,sigma)
              obj@cost <- obj@cost + tmp@cost
              obj@point <- c(obj@point,tmp)
              return(obj)
          })

#' Add a Collective Anomaly segment
#' @param obj partition object
#' @param f function to generate a new class segment
#' @param b penalty for introducing new segment
#' @param t time at start of segment
setMethod("addCollective","partition",
          function(obj,f,b,t){
              tmp <- f(b,t)
              obj@collective <- c(obj@collective, tmp)
              obj@cost <- obj@cost + tmp@cost
              obj@isAnom <- TRUE
              return(obj)
          })

#' Add a new baseline period
#' @param obj partition object
#' @param f function to generate a new class segment
#' @param b penalty for introducing new segment
#' @param t time at start of segment
setMethod("addBase","partition",
          function(obj,f,b,t){
              tmp <- f(b,t)
              obj@collective <- c(obj@collective, tmp)
              obj@cost <- obj@cost + tmp@cost
              obj@isAnom <- FALSE
              return(obj)
          })


#' create a new partition
## #' @param min_length minimum length of a non-point segment
## #' @param max_length maximum length of a non-point segment
#' 
#' @export
partition <- function(){
    out <- new("partition")
    return( out )
}


#' Update a partion object
#' @param obj partition object
#' @param x observed value
#' @param mu mean of background distribution
#' @param sigma variance of background distribution
#' @param isPoint add data to the most recent point
setMethod("update","partition",
          function(obj,x,mu,sigma){
              nc <- length(obj@collective)
              tmp <- obj@collective[[nc]]@cost
              obj@collective[[nc]] <- update(obj@collective[[nc]],x,mu,sigma)
              obj@last_n <- as.integer(obj@collective[[nc]]@n)
              obj@cost <- obj@cost - tmp + obj@collective[[nc]]@cost
              ##obj@cost <- sum( as.numeric(sapply(obj@collective,function(x){x$cost})) ) +
              ##    sum( as.numeric(sapply(obj@point,function(x){x$cost})) )              
              return(obj)
          })


#' return summary of Collective Changes
#' @param p a partition object
#' @export
collective_anomalies <- function(p,type=NULL){
    stopifnot(
        "p should be of class partition" = inherits(p,"partition")
    )
    fDF <- function(x){
        data.frame(start=x@start, end = x@start + x@n - 1, type= class(x), t(x@param))
    }
    
    out <- do.call( rbind, lapply(p@collective,fDF) )

    if(!is.null(type)){
        out <- out[out$type==type,]
    }
    
    return(out)
}


#' return summary of point anomalies
#' @param p a partition object
#' @export
point_anomalies <- function(p){
    stopifnot(
        "p should be of class partition" = inherits(p,"partition")
    )

    fDF <- function(x){data.frame(location = x@start, type=class(x), t(x@param))}
    
    out <- do.call(rbind,lapply(p@point,fDF))
##    out <- Reduce(function(x,y){merge(x,y,all=T)},tmp)

    return(out)
}
