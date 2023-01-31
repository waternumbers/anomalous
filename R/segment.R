## generics

## TODO - work through gaussian cases to fix where variance terms are assumed to be either 1 or constant
## TODO - Adapt create function to return a blank example then update outside in pert + op. Thi is to make capa easier!
## TODO - work through capa, particularly check all correct updates applied - three moves for everything in ctlg
## TODO - integrate capa examples
## TODO - change data inputs to matrix so can do multivariate distributions :-)
## TODO - try insertion into all segments to get a classifier?



## basic segment class
setClass("Segment", 
         slots = c(n = "integer",
                   min_length = "integer",
                   max_length = "integer",
                   ss = "numeric",
                   beta = "numeric",
                   cost = "numeric"),
         prototype = list(
             n = NA_integer_,
             min_length = NA_integer_,
             max_length = NA_integer_,
             ss = NA_real_,
             beta = NA_real_,
             cost = NA_real_)
         )

setMethod("cost","Segment",function(obj){return(obj@cost)})




## gauss mean
#' @include generics.R
setClass("meanSegment",contains = "Segment")

#' meanSegment
#' @param x initial data
#' @param beta Penalisation term for the segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#' @export
meanSegment <- function(x,beta,min_length,max_length){
    out <- new("meanSegment",
               n = integer(1),
               ss = numeric(2),
               min_length = as.integer(min_length)[1],
               max_length = as.integer(max_length)[1],
               beta = as.numeric(beta)[1])
    return( update(out,x) )
}

#' @rdname update-methods
setMethod("update","meanSegment",
          function(obj,x){
              obj@n <- obj@n + length(x)
              obj@ss <- obj@ss + c(sum(x), sum(x^2))
              ## compute cost
              if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ obj@cost <- Inf }
              else{ obj@cost <- obj@beta + obj@ss[2] - ( (obj@ss[1]^2)/obj@n ) }
              
              return(obj)
          })

## #' @rdname cost-methods
## setMethod("cost","meanSegment",function(obj){obj@cost})
##           ## function(obj){
##           ##     if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ return(Inf) }
##           ##     return( obj@beta + obj@ss[2] - ( (obj@ss[1]^2)/obj@n ) )
##           ## })

## gauss variance
setClass("varSegment",contains = "Segment")

#' varSegment
#' @param x initial data
#' @param beta Penalisation term for the segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#' @export
varSegment <- function(x,beta,min_length,max_length){
    out <- new("varSegment",
               n = integer(1),
               ss = numeric(1),
               min_length = as.integer(min_length)[1],
               max_length = as.integer(max_length)[1],
               beta = as.numeric(beta)[1])
    return( update(out,x) )
}

#' @rdname update-methods
setMethod("update","varSegment",
          function(obj,x){
              obj@n <- obj@n + length(x)
              obj@ss <- obj@ss + sum(x^2)

              ## compute cost
              if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ obj@cost <- Inf }
              else{ obj@cost <- obj@beta + obj@n*(log(2*pi) + log( obj@ss/obj@n ) +1 ) }
              
              return(obj)
          })

## #' @rdname cost-methods
## setMethod("cost","varSegment",function(obj){obj@cost})
##           ## function(obj){
##           ##     if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ return(Inf) }
##           ##     return( obj@beta + obj@n*(log(2*pi) + log( obj@ss/obj@n ) +1 ) )
##           ## })

## #########################################################################################
## Guass var and mean
setClass("meanvarSegment",contains = "Segment")

#' meanvarSegment
#' @param x initial data
#' @param beta Penalisation term for the segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#' @export
meanvarSegment <- function(x,beta,min_length,max_length){
    out <- new("meanvarSegment",
               n = integer(1),
               ss = numeric(2),
               min_length = as.integer(min_length)[1],
               max_length = as.integer(max_length)[1],
               beta = as.numeric(beta)[1])
    return( update(out,x) )
}

#' @rdname update-methods
setMethod("update","meanvarSegment",
          function(obj,x){
              obj@n <- obj@n + length(x)
              obj@ss <- obj@ss + c(sum(x), sum(x^2))
              
              ## compute cost
              if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ obj@cost <- Inf }
              else{
                  out <- obj@ss[2] - ( (obj@ss[1]^2)/obj@n )
                  out <- obj@n*(log(2*pi) + log( out / obj@n ) +1 )
                  obj@cost <- out + obj@beta
              }
              return(obj)
          })

## #' @rdname cost-methods
## setMethod("cost","meanvarSegment",function(obj){obj@cost})
##           ## function(obj){
##           ##     if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ return(Inf) }
##           ##     out <- obj@ss[2] - ( (obj@ss[1]^2)/obj@n )
##           ##     out <- obj@n*(log(2*pi) + log( out / obj@n ) +1 )
##           ##     return(out + obj@beta)
##           ## })

## ###############################################################################################
## poission segment
setClass("poisSegment",contains = "Segment")

#' poisSegment
#' @param x initial data
#' @param beta Penalisation term for the segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#' @export
poisSegment <- function(x,beta,min_length,max_length){
    out <- new("poisSegment",
               n = integer(1),
               ss = numeric(1),
               min_length = as.integer(min_length)[1],
               max_length = as.integer(max_length)[1],
               beta = as.numeric(beta)[1])
    return( update(out,x) )
}


#' @rdname update-methods
setMethod("update","poisSegment",
          function(obj,x){
              obj@n <- obj@n + length(x)
              obj@ss <- obj@ss + sum(x)

              ## compute cost
              if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ obj@cost <- Inf }
              else{
                  if( obj@ss == 0){ out <- 0 }
                  else{ out <- 2*obj@ss*(log(obj@n)-log(obj@ss)) }
                  obj@cost <- out + obj@beta
              }
              
              return(obj)
          })
## #' @rdname cost-methods
## setMethod("cost","poisSegment",function(obj){obj@cost})
          ## function(obj){
          ##     if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ return(Inf) }
          ##     if( obj@ss == 0){ out <- 0 }
          ##     else{ out <- 2*obj@ss*(log(obj@n)-log(obj@ss)) }
          ##     return( out + obj@beta )
          ## })

## ########################################################################################
## zerooneSegment
setClass("zerooneSegment",contains = "Segment")

#' zerooneSegment
#' @param x initial data
#' @param beta Penalisation term for the segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#' @export
zerooneSegment <- function(x,beta,min_length,max_length){
    out <- new("zerooneSegment",
               n = integer(1),
               ss = numeric(1),
               min_length = as.integer(min_length)[1],
               max_length = as.integer(max_length)[1],
               beta = as.numeric(beta)[1])
    return( update(out,x) )
}


#' @rdname update-methods
setMethod("update","zerooneSegment",
          function(obj,x){
              obj@n <- obj@n + length(x)
              obj@ss <- obj@ss + sum(x)

              ## compute cost
              if( (obj@n < obj@min_length) | (obj@n > obj@max_length) ){ obj@cost <- Inf }
              else{

                  obj@cost <- obj@beta +
                  if( obj@ss == 0){ out <- 0 }
                  else{ out <- 2*obj@ss*(log(obj@n)-log(obj@ss)) }
                  obj@cost <- out + obj@beta
              }
              
              return(obj)
          })


## #' A wrapper function to create segments
## #' @param segType type of segment
## #' @param beta Penalisation term for the segment
## #' @param min_length minimum length of a segment
## #' @param max_length maximum length of a segment
## #' 
## #' @details Basic R implimentation - not efficent
## #' @export
## createSegment <- function(segType,beta,min_length,max_length){
##     ##browser()
##     ## TODO write as proper tests??
##     beta <- as.numeric(beta)[1]
##     min_length <- as.integer(min_length)[1]
##     max_length <- as.integer(max_length)[1]
##     ss <- switch(segType,
##                  "meanSegment" = numeric(2),
##                  "varSegment" = numeric(1),
##                  "meanVarSegment" = numeric(2),
##                  "poisSegment" = numeric(1),
##                  stop("Unknown Segment type")
##                  )
##     return( new(segType,n=integer(1),ss=ss,min_length=min_length,max_length=max_length,beta=beta) )
## }

## #' @export
## segment <-
##     R6::R6Class("segment",
##                 public = list(
##                     initialize = function(x,beta){
##                         self$beta <- beta
##                         self$update(x)
##                     },
##                     n = 0,
##                     ss = NA,
##                     beta = NA,
##                     update = function(x){
##                         self$n <- self$n + 1
##                         self$updateSS(x)
##                         invisible(self)
##                     },
##                     cost = function(min_length=2, max_length=Inf){
##                         if( (self$n < min_length) | (self$n > max_length) ){ return(Inf) }
##                         out <- self$computeCost() + self$beta
##                         return(out)
##                     },
##                     computeCost = function(){Inf},
##                     updateSS = function(){}
##                  )
##                 )

## #' @export
## seg_mean <-
##     R6::R6Class("seg_mean",
##                 inherit = segment,
##                 public = list(
##                     ss = rep(0,2),
##                     updateSS = function(x){
##                         self$ss <- self$ss + c(x, x*x)
##                     },
##                     computeCost = function(){
##                         return( self$ss[2] - ( (self$ss[1]^2)/self$n ) )
##                     }
##                 )
##                 )

## #' @export
## seg_var <-
##     R6::R6Class("seg_var",
##                 inherit = segment,
##                 public = list(
##                     ss = 0,
##                     updateSS = function(x){
##                         self$ss <- self$ss + x*x
##                     },
##                     computeCost = function(){
##                         return( self$n*(log(2*pi) + log( self$ss/self$n ) +1 ) )
##                     }
##                 )
##                 )
## #' @export
## seg_mean_var <- 
##     R6::R6Class("seg_mean_var",
##                 inherit = segment,
##                 public = list(
##                     ss = rep(0,2),
##                     updateSS = function(x){
##                         self$ss <- self$ss + c(x,x*x)
##                     },
##                     computeCost = function(){
##                         out <- self$ss[2] - ( (self$ss[1]^2)/self$n )
##                         out <- self$n*(log(2*pi) + log( out / self$n ) +1 )
##                         return(out)
##                     }
##                 )
##                 )
