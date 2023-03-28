##Rcpp::loadModule("gaussSeg",TRUE)

#' meanSegment
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussMean <- function(beta,t){
    out <- new(gaussMeanCpp, start=as.integer(t), penalty=beta)
    return( out )
}


## ####################################################################################
## gauss variance
#' gaussVar
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussVar <- function(beta,t){
    out <- new(gaussVarCpp, start=as.integer(t), penalty=beta)
    return( out )
}


## ####################################################################################
## gauss mean & variance
#' gaussMeanVar
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussMeanVar <- function(beta,t){
    out <- new(gaussMeanVarCpp, start=as.integer(t), penalty=beta)
    return( out )
}

## ####################################################################################
## gauss fixed
#' meanSegment
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussFixed <- function(beta,t){
    out <- new(gaussFixedCpp, start=as.integer(t), penalty=beta)
    return( out )
}


## ###################################################################################
#' gaussPoint
#' @param beta Penalisation term for the segment
#' @param t start time of segment
#' 
#' @export
gaussPoint <- function(beta,t){
    out <- new(gaussPointCpp, start=as.integer(t), penalty=beta)
    return( out )
}
