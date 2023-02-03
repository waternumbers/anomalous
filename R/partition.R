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
         slots = c(collective = "data.frame",
                   point = "data.frame",
                   is_valid = "logical",
                   beta = "numeric",
                   cost = "numeric",
                   min_length = "integer",
                   max_length = "integer"),
         prototype = list(
             collective = as.data.frame(NULL),
             point = as.data.frame(NULL),
             is_valid = FALSE,
             beta = NA_real_,
             cost = NA_real_,
             min_length = NA_integer_,
             max_length = NA_integer_
         )
         )

setMethod("cost","Partition",function(obj){return(obj@cost)})


## gauss partition
#' @include generics.R
setClass("gaussPartition",contains = "Partition")


#' gaussPartition
#' @param beta Named vector of penalisation terms for starting a new segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#' @export
gaussPartition <- function(beta,min_length,max_length){           
    out <- new("gaussPartition",
               beta = beta,
               min_length = as.integer(min_length)[1],
               max_length = as.integer(max_length)[1])
    return( out )
}

setMethod("update","gaussPartition",
          function(obj,x,mu,sigma,t,type=c("current","background","mean","meanvar","var","point")){

              type <- match.arg(type)
              
              if( type=="point" ){ stop("Not doing points yet...") }

              if( !(type %in% c("current","point")) ){
                  ## start a new time segment
                  obj@collective <- rbind(
                      obj@collective,
                      data.frame(type = type,
                                 strt  = t,
                                 n = integer(1),
                                 sum_w = numeric(1),
                                 sum_log_sig = numeric(1),
                                 sum_eta = numeric(1),
                                 sum_eta2 = numeric(1),
                                 mhat = NA_real_,
                                 shat = NA_real_,
                                 cost = numeric(1),
                                 pen = obj@beta[type]
                                 )                
                  )
              }

              ## group to add data to
              ii <- nrow( obj@collective )

              obj@collective$n[ii] <- obj@collective$n[ii] + 1
              obj@collective$sum_w[ii] <- obj@collective$sum_w[ii] + 1/sigma
              obj@collective$sum_log_sig[ii] <- obj@collective$sum_log_sig[ii] + log(sigma)
              obj@collective$sum_eta[ii] <- obj@collective$sum_eta[ii] + (x-mu)/sigma
              obj@collective$sum_eta2[ii] <- obj@collective$sum_eta2[ii] + ((x-mu)^2)/sigma
              obj@collective$mhat[ii] <- switch(obj@collective$type[ii],
                                                "background" = 0,
                                                "var" = 0,
                                                obj@collective$sum_eta[ii] / obj@collective$sum_w[ii])
              kappa <- ( obj@collective$sum_eta2[ii] - 2*obj@collective$mhat[ii]*obj@collective$sum_eta[ii] +
                         (obj@collective$mhat[ii]^2)*obj@collective$sum_w[ii]) / obj@collective$n[ii]
              kappa <- max(kappa,2*.Machine$double.xmin) ## to catch when s is zero...
              obj@collective$shat[ii] <- switch(obj@collective$type[ii],
                                                "background" = 1,
                                                "mean" = 1,
                                                kappa)
              obj@collective$cost[ii] <- obj@collective$n[ii]*log(2*pi*obj@collective$shat[ii]) +
                  obj@collective$sum_log_sig[ii] + (1/obj@collective$shat[ii])*obj@collective$n[ii]*kappa
 
              ## correct cost for length
              if( (obj@collective$n[ii] < obj@min_length ) | (obj@collective$n[ii] > obj@max_length ) ){
                  obj@collective$cost[ii] <- Inf
              }

              obj@is_valid <- TRUE
              if( (obj@collective$n[ii] < obj@min_length ) & is.finite(sum( obj@collective$cost[-ii] )) ){
                  obj@is_valid <- FALSE
              }
                            
              ## update total cost
              obj@cost <- sum( obj@collective$cost ) + sum( obj@collective$pen )

              return(obj)
          })

