#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param fC cost function
#' @param Beta Penalisation term
#'
#' @details Basic R implimentation - not efficent
#'@export
op <- function(y,fC,beta){
    n <- length(y)
    F <- rep(NA,n+1)
    cp <- rep(list(rep(NA,n)),n)

    ## handle tt ==1
    cp[[1]][1] <- 1
    F[1] <- fC(y,cp[[1]],beta)
    
    
    
    for(tt in 2:n){
        ## first case is to keep the same
        cl <- cp[[tt-1]]
        cl[tt] <- cl[tt-1]
        F[tt] <- fC(y,cl,beta)
        cp[[tt]] <- cl
        
        ## search impact of extending all previous groups
        for(tau in 1:(tt-1)){
            cl <- cp[[tau]]
            cl[(tau+1):tt] <- cl[tau]+1
            tmp <- fC(y,cl,beta)
            if( tmp < F[tt] ){
                #browser()
                F[tt] <- tmp
                cp[[tt]] <- cl
            }
        }
    }
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


