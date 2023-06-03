#' Test computations of the cost
#'
#' Note that the point anomaly is not correct...
#' @export
actual_mean_cost <- function(x,ca,pa,beta,betaP){
    ## handle the base distribution
    bdx <- 1:length(x)
    bdx <- setdiff(bdx,pa$location)
    tmp <- unlist(lapply(1:nrow(ca),function(ii){ca$start[ii]:ca$end[ii]}))
    bdx <- setdiff(bdx,tmp)
    cost_base <- -2*sum(dnorm(x[bdx],0,1,log=TRUE))
    ## collective anomaly
    cost_col <- 0
    for(ii in 1:nrow(ca)){
        idx <- ca$start[ii]:ca$end[ii]
        cost_col <- cost_col + -2*sum(dnorm(x[idx],mean(x[idx]),1,log=TRUE)) + beta
    }
    ## point anomalies
    cost_point <- -2*sum(dnorm(x[pa$location],0,abs(pa$location),log=TRUE) + betaP)
    return( cost_base + cost_point + cost_col )
}
