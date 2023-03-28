## generics

## TODO - change data inputs to matrix so can do multivariate distributions :-)
## TODO - try insertion into all segments to get a classifier?

## basic segment class
#' @include generics.R
setClass("partition", 
         slots = c(collective = "list",
                   point = "list",
                   last_n = "integer",
                   cost = "numeric"),
         prototype = list(
             collective = list(),
             point = list(),
             last_n = NA_integer_,
             cost = NA_real_
         )
         )



#' Cost of the partioning
#' @param obj partition object
setMethod("cost","partition",function(obj){return(obj@cost)})

#' Add a new collective parrtion
#' @param obj partition object
#' @param f function to generate a new class segment
#' @param b penalty for introducing new segment
#' @param t time at start of segment
setMethod("addCollective","partition",
          function(obj,f,b,t){
              obj@collective <- c(obj@collective, f(b,t))
              return(obj)
          })

#' Add a new point partition
#' @param obj partition object
#' @param f function to generate a new class segment
#' @param b penalty for introducing new segment
#' @param t time at start of segment
setMethod("addPoint","partition",
          function(obj,f,b,t){
              force(obj)
              obj@point <- c(obj@point, f(b,t))
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
          function(obj,x,mu,sigma,isPoint=FALSE){
              np <- length(obj@point)
              nc <- length(obj@collective)
              if( isPoint ){
                  obj@point[[np]]$update(x,mu,sigma)
              }else{
                  obj@collective[[nc]]$update(x,mu,sigma)
              }
              obj@last_n <- as.integer(obj@collective[[nc]]$n)
              obj@cost <- sum( as.numeric(sapply(obj@collective,function(x){x$cost})) ) +
                  sum( as.numeric(sapply(obj@point,function(x){x$cost})) )
              
              return(obj)
          })


#' return summary of Collective Changes
#' @param p a partition object
#' @export
collective_change <- function(p){
    stopifnot(
        "p should be of class partition" = inherits(p,"partition")
    )
    
    fDF <- function(x){data.frame(start=x$start, segment=class(x), t(x$param))}
    
    tmp <- lapply(p@collective,fDF)
    out <- Reduce(function(x,y){merge(x,y,all=T)},tmp)
    out$end <- c(out$start[-1]-1,Inf)
    
    idx <- c("start","end","segment")
    idx <- c(idx,setdiff(names(out),idx))
    out <- out[,idx]
    
    return(out)
}

#' return summary of point anomalies
#' @param p a partition object
#' @export
point_change <- function(p){
    stopifnot(
        "p should be of class partition" = inherits(p,"partition")
    )

    fDF <- function(x){data.frame(index=x$start, segment=class(x), t(x$param))}
    
    tmp <- lapply(p@point,fDF)
    out <- Reduce(function(x,y){merge(x,y,all=T)},tmp)

    return(out)
}
