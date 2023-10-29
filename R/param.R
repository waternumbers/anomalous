#' Get the parameters for a partitioning result
#' @param res the result of a partitioning algorithm
#' @param fCost the cost function
#'
#' @return list of parameters
#' 
#' @details Not yet implimented for all cost functions
#' @export
param <- function(res,fCost){

    if( !("param" %in%  names(fCost)) ){ stop("Parameter method not available for cost function") }
    list(ca = lapply(res$ca, function(v){ fCost$param(v["start"],v["end"]) }),
         pa = lapply(res$pa, function(v){ fCost$param(v["location"],v["location"]) })
         )
}
