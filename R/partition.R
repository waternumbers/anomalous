## generics

## TODO - work through gaussian cases to fix where variance terms are assumed to be either 1 or constant
## TODO - Adapt create function to return a blank example then update outside in pert + op. Thi is to make capa easier!
## TODO - work through capa, particularly check all correct updates applied - three moves for everything in ctlg
## TODO - integrate capa examples
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

setMethod("cost","partition",function(obj){return(obj@cost)})

#' meanSegment
#' @param pen Penalisation term for the segment
#' 
#' @export
partition <- function(min_length,max_length){
    out <- new("partition")
    return( out )
}

setMethod("update","partition",
          function(obj,x,mu,sigma,seg=NULL,isPoint=FALSE){

              nc <-  length(obj@collective)
              ncp <- length(obj@point)
              if( is.null(seg) ){
                  obj@collective[[nc]] <- update(obj@collective[[nc]],x,mu,sigma)
              }else{
                  if( isPoint ){
                      ncp <- ncp+1
                      obj@point[[ncp]] <- update(seg,x,mu,sigma)
                  }else{
                      nc <- nc+1
                      obj@collective[[nc]] <- update(seg,x,mu,sigma)
                  }
              }

              if(nc>0){
                  collective_cost <- sapply(obj@collective,function(x){x@cost})
                  obj@last_n <- obj@collective[[nc]]@n
              }else{
                  collective_cost <- 0
                  obj@last_n <- 0
              }

              if(ncp>0){
                  point_cost <- sapply(obj@point,function(x){x@cost})
              }else{
                  point_cost <- 0
              }

              obj@cost <- sum(collective_cost) + sum(point_cost)

              return(obj)
          })
