#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
#' @param mu if univariate series y
#' @param sigma variance of univariate series y
#' @param segType function to generate an object of class Segment
#' @param beta Penalisation term for a new segment
#' @param min_length minimum length of a segment
#' @param max_length maximum length of a segment
#' 
#' @details Basic R implimentation - not efficent
#'@export
op <- function(y,mu,sigma,segType,beta,min_length=2,max_length=.Machine$integer.max){
    ##browser()
    n <- length(y)
    opt <- rep(list(NULL),n) ## optimal partions

    ctlg <- rep(list(NULL),n) ## catalog of partitions to try
    last_n <- rep(NA,n)
    cst <- rep(NA,n)
    
    ## initialise
    p <- addCollective(partition(),collectiveSeg,beta,1)
    for(tt in 1:min_length){
        p <- update(p,y[tt],mu[tt],sigma[tt])
    }

    ctlg[[min_length]] <- p
    opt[[min_length]] <- p
        
    ## loop time
    last_n <- rep(NA,n)
    cst <- rep(NA,n)
    for(tt in (min_length+1):n){
        p <- force( opt[[tt-1]] )
        p <- addCollective(p,segType,beta,tt)
        ctlg[[tt]] <- force(p)
        
        
        for(tau in min_length:tt){
            ctlg[[tau]] <- update( ctlg[[tau]], y[tt], mu[tt],sigma[tt] )
            last_n[tau] <- ctlg[[tau]]@last_n
            cst[tau] <- ctlg[[tau]]@cost
        }

        idx <- min_length > last_n | max_length < last_n
        cst[idx] <- Inf

        opt[[tt]] <- force( ctlg[[which.min(cst)]] )
    }
    
    ##     cst <- ctlg[[tt]]@cst
    ##     if( min_length > ctlg[[tt]]@last_n | max_length < ctlg[[tt]]@last_n ){cst <- Inf}
        
    ##     opt[[tt]] <- ctlg[[tt]]
    ##     cst <- opt[[tt]]@cst
    ##     if( opt[
        
    ##     for(tau in 1:(tt-1)){

    ##         ctlg[[tau]] <- update( ctlg[[tau]], y[tt], mu[tt],sigma[tt] )

    ##         tmp <- ctlg[[tau]]@cost
    ##         if( tmp < cst[tt] ){
    ##             opt[[tt]] <- ctlg[[tau]]
    ##             cst[tt] <- tmp
    ##         }
    ##     }
    ## }
    browser()
    return(opt[[tt]])
}
