## ####################################################################################
## poission
#' @include generics.R
setClass("poisRate",contains = "Segment",
         prototype = list(
             summaryStats = c(sum_x=0),
             param = c(lambda = NA_real_)
         ))


#' gaussMeanVar
#' @param pen Penalisation term for the segment
#' 
#' @export
poisRate <- function(beta,t){
    out <- new("poisRate", start=as.integer(t), penalty=beta)
    return( out )
}

#' @rdname update-methods
setMethod("update","poisRate",
          function(obj,x,mu,sigma){
              obj@n <- obj@n + as.integer(1)
              ss <- obj@summaryStats
              ss <- ss + x
              obj@summaryStats <- ss
              ## TODO obj@param["lambda"] <- ?? ss["sum_eta"] / ss["sum_w"]
              if( ss == 0){ out <- 0 }
              else{ out <- 2*ss*(log(obj@n)-log(ss)) }

              obj@cost <- out + obj@penalty

              return(obj)
          }
          )


