#' An R implimentation of the segmented search algorithm
#' @param y univariate data series
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
    p <- partition() ##min_length,max_length)
    for(tt in 1:min_length){
        if(tt==1){ seg <- segType(beta,1) }else{seg <- NULL}
       
        p <- update(p,y[1],mu[1],sigma[1],seg)
    }

    ctlg[[min_length]] <- p
    opt[[min_length]] <- p
        
    ## loop time
    last_n <- rep(NA,n)
    cst <- rep(NA,n)
    for(tt in (min_length+1):n){

        ctlg[[tt]] <- update(opt[[tt-1]], y[tt],mu[tt],sigma[tt], segType(beta,tt))
        last_n[tt] <- ctlg[[tt]]@last_n
        cst[tt] <- ctlg[[tt]]@cost
        
        for(tau in min_length:(tt-1)){
            ctlg[[tau]] <- update( ctlg[[tau]], y[tt], mu[tt],sigma[tt] )
            last_n[tau] <- ctlg[[tau]]@last_n
            cst[tau] <- ctlg[[tau]]@cost
        }

        idx <- min_length > last_n | max_length < last_n
        cst[idx] <- Inf

        opt[[tt]] <- ctlg[[which.min(cst)]]
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
    
    return(opt[[tt]])
}
