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
                   is_valid = "logical",
                   cost = "numeric",
                   min_length = "integer",
                   max_length = "integer"),
         prototype = list(
             collective = list(),
             point = list(),
             is_valid = FALSE,
             cost = NA_real_,
             min_length = NA_integer_,
             max_length = NA_integer_
         )
         )

setMethod("cost","partition",function(obj){return(obj@cost)})

#' meanSegment
#' @param pen Penalisation term for the segment
#' 
#' @export
partition <- function(min_length,max_length){
    out <- new("partition",
               min_length = as.integer(min_length)[1],
               max_length = as.integer(max_length)[1])
    return( out )
}


setMethod("update","partition",
          function(obj,x,mu,sigma,seg=NULL){

              nc <-  length(obj@collective)
              if( is.null(seg) ){
                  obj@collective[[nc]] <- update(obj@collective[[nc]],x,mu,sigma)
              }else{
                  nc <- nc+1
                  obj@collective[[nc]] <- update(seg,x,mu,sigma)
              }

              costVec <- rep(NA,nc)
              nVec <- rep(NA,nc)
              for(ii in 1:nc){
                  costVec[ii] <- obj@collective[[ii]]@cost
                  nVec[ii] <- obj@collective[[ii]]@n
              }
              costVec[ nVec<obj@min_length | nVec>obj@max_length ] <- Inf
              
              obj@cost <- sum(costVec)
              
              obj@is_valid <- TRUE
              if( all(is.finite(costVec[-nc])) & (nVec[nc]<obj@min_length) ){
                  obj@is_valid <- FALSE
              }

              return(obj)
          })
