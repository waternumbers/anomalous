#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param seg cost function
#' @param beta Penalisation term
#'
#' @details Basic R implimentation - not efficent
#'@export
op <- function(y,seg,beta){
    n <- length(y)
    F <- rep(NA,n) ## don't actually need this
    cp <- rep(list(NA),n) ## catalog of partitions to try

    
    fC <- function(p){
        sum( sapply(p,function(x){x$cost()}) )
    }

    ## handle tt ==1
    cp[[1]] <- list( seg$new(y[1],0) ) ## no penalty for the first segment
    F[1] <- fC(cp[[1]])
    
    for(tt in 2:length(y)){
        ## first case is to keep the same
        p <- cp[[tt-1]]
        p <- c(p,seg$new(y[tt],beta))
        ## p[[length(p)]]$update(y[tt])
        F[tt] <- fC(p)
        cp[[tt]] <- p

        ## search impact of extending all previous groups
        for(tau in 1:(tt-1)){
            p <- cp[[tau]]
            p[[length(p)]]$update(y[tt])
            ## p <- c(p,seg$new(y[tt],beta))
            tmp <- fC(p)

            if( tmp < F[tt] ){
                #browser()
                F[tt] <- tmp
                cp[[tt]] <- p
            }
        }
    }
    browser()
    cp <- unlist(tail(cp,1)) ## final classification
    return(cp)
    ##return( which(diff(cp)>0) )
}



## op <- function(y,fC,Beta){
##     n <- length(y)
##     F <- rep(NA,n+1)
##     cp <- rep(list(NULL),n+1)
##     F[1] <- -Beta
##     for(tt in 1:n){
##         Fvec <- rep(NA,tt)
##         for(tau in 0:(tt-1)){
##             jj <- tau+1
##             Fvec[jj] <- F[jj] + fC(y[jj:tt]) + Beta
##         }
##         tauhat <- which.min(Fvec)-1
##         cp[[tt+1]] <- c(cp[[tauhat+1]],tauhat)
##         F[tt+1] <- Fvec[tauhat+1]
##     }
##     cp <- unlist(tail(cp,1))[-1] ## final value and trim intial 0
##     return(cp)
## }


