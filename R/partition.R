## generics

## TODO - work through gaussian cases to fix where variance terms are assumed to be either 1 or constant
## TODO - Adapt create function to return a blank example then update outside in pert + op. Thi is to make capa easier!
## TODO - work through capa, particularly check all correct updates applied - three moves for everything in ctlg
## TODO - integrate capa examples
## TODO - change data inputs to matrix so can do multivariate distributions :-)
## TODO - try insertion into all segments to get a classifier?

## basic segment class
#' @include generics.R
setClass("Partition", 
         slots = c(idx = "integer",
                   beta = "numeric",
                   cost = "numeric"),
         prototype = list(
             idx = integer(0),
             beta = rep(NA_real_,3),
             cost = NA_real_)
         )

setMethod("cost","Partition",function(obj){return(obj@cost)})
